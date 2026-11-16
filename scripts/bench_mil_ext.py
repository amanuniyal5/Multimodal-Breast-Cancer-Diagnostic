"""ABMIL variants for X3 / X4b / X4c / X4d / X6 (docs/DESIGN_SPEC.md §5.3, 5.4, 5.6). Same runner settings as
bench_mil.py (fixed site-grouped folds, Adam lr 2e-4 wd 1e-4, 15 epochs, 2,048 random patches per step, all <= cap at
test, fp16 autocast, class-weighted BCE), 3 seeds.

--mode regress    X3: ABMIL with a linear output and MSE on the z-scored target (training-fold mean / SD) for
                  ESR1_log2, PGR_log2, ERBB2_log2. Rows: exp "reg", p = prediction on the original log2 scale.
--mode multitask  X4b: one shared ABMIL, heads ER_IHC / PR_IHC / HER2_STRICT (BCE), PAM50 (5-way CE), ESR1 / PGR /
                  ERBB2 (MSE on z-scores); missing targets masked, equal weights. Patients: all with >= 1 target.
--mode single     X4b reference: single-task ABMIL for ER_IHC / PR_IHC / HER2_STRICT on the same patients as multitask
                  (every patient with that label), so the two differ only in the extra heads.
--mode dann       X4c(ii): ABMIL + gradient-reversal site classifier (lam ramped 0 -> 0.1) vs the same network without
                  the site loss; ER group (ER_IHC, ESR1_RNApos) and HER2_STRICT. Test-patient embeddings (seed 0) are
                  saved to <out>/dann_embeddings.npz for the hospital-ID check.
--mode noisy      X4d: ER group (801 patients); targets for training: bce (IHC), consensus (train only on IHC = RNA),
                  soft (latent-class posterior from IHC / RNA / RPPA, fitted on the training fold only), gce (q = 0.7,
                  IHC), smooth (label smoothing 0.1, IHC). Scored later against ER_IHC, ESR1_RNApos and ER_RPPApos.
--mode reliab     X6: ER_IHC, PR_IHC, HER2_STRICT; at test also 20 MC-dropout passes and 20 random half-bags per patient
                  (mean and SD of the probability).
Output: <out>/<mode>_<name>_oof.csv.gz (mode, variant, exp, task, seed, patient, y, p, + extra columns).
Usage: python scripts/bench_mil_ext.py --mode noisy [--name uni_v1_20x] [--smoke]
"""
import argparse
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bench import tasks  # noqa: E402
from bench.latent import dawid_skene  # noqa: E402
from common import site_folds  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402

warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser()
ap.add_argument("--mode", required=True, choices=["regress", "multitask", "single", "dann", "noisy", "reliab"])
ap.add_argument("--name", default="uni_v1_20x")
ap.add_argument("--feats", default="data/features")
ap.add_argument("--cap", type=int, default=4096)
ap.add_argument("--seeds", type=int, default=3)
ap.add_argument("--epochs", type=int, default=15)
ap.add_argument("--lr", type=float, default=2e-4)
ap.add_argument("--train-patches", type=int, default=2048)
ap.add_argument("--passes", type=int, default=20, help="MC-dropout passes and half-bags (reliab)")
ap.add_argument("--smoke", action="store_true", help="150 patients, 1 seed, 1 epoch")
ap.add_argument("--out", default="results/bench/x4")
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
t0 = time.time()
SEEDS = [0] if a.smoke else list(range(a.seeds))
EPOCHS = 1 if a.smoke else a.epochs


def log(*x):
    print(f"[{time.time() - t0:6.0f}s]", *x, flush=True)


import torch  # noqa: E402
import torch.nn as nn  # noqa: E402
import torch.nn.functional as F  # noqa: E402

from mil_models import ABMIL, ABMILDANN, MultiTaskABMIL  # noqa: E402

dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
if dev.type != "cuda" and os.environ.get("SLURM_JOB_ID"):
    sys.exit(f"no usable GPU on {os.environ.get('SLURMD_NODENAME')}: refusing to train on CPU")
AMP = dev.type == "cuda"
idx = pd.read_csv(os.path.join(a.feats, f"{a.name}_index.csv")).set_index("patient")
L = tasks.labels()
L = L[L.index.isin(idx.index)]
if a.smoke:
    L = L[L.index.isin(np.random.default_rng(0).choice(L.index, 150, replace=False))]
t = time.time()
BAGS = np.load(os.path.join(a.feats, f"{a.name}_cap{a.cap}.npy"), mmap_mode="r" if a.smoke else None)
D = BAGS.shape[2]
log(f"mode {a.mode}; {len(L)} patients; bags {BAGS.shape} loaded in {time.time() - t:.0f}s; device {dev}")


def bag(p, train, rng, half=False):
    r, k = idx.at[p, "row"], idx.at[p, "n_kept"]
    x = BAGS[r, :k]
    if train and k > a.train_patches:
        x = x[rng.choice(k, a.train_patches, replace=False)]
    elif half and k > 1:
        x = x[rng.choice(k, k // 2, replace=False)]
    return torch.from_numpy(np.ascontiguousarray(x)).to(dev, non_blocking=True).float()


def train(net, tr, step_loss, seed, on_step=None):
    """generic loop: step_loss(net, i, x) -> loss for training patient position i"""
    rng = np.random.default_rng(seed)
    opt = torch.optim.Adam(net.parameters(), lr=a.lr, weight_decay=1e-4)
    scaler = torch.cuda.amp.GradScaler(enabled=AMP)
    total, step = EPOCHS * len(tr), 0
    for _ in range(EPOCHS):
        net.train()
        for i in rng.permutation(len(tr)):
            if on_step:
                on_step(net, i, step / total)
            with torch.autocast(device_type=dev.type, dtype=torch.float16, enabled=AMP):
                loss = step_loss(net, i, bag(tr[i], True, rng))
            opt.zero_grad(set_to_none=True)
            scaler.scale(loss).backward()
            scaler.step(opt)
            scaler.update()
            step += 1
    net.eval()
    return net


def infer(net, pts, seed, fn=lambda o: torch.sigmoid(o[0].float()).item(), mode="eval"):
    rng = np.random.default_rng(seed + 777)
    out = []
    if mode == "mc":
        net.train()
    with torch.no_grad(), torch.autocast(device_type=dev.type, dtype=torch.float16, enabled=AMP):
        for p in pts:
            out.append(fn(net(bag(p, False, rng, half=(mode == "half")))))
    net.eval()
    return out


def bce_weight(y):
    pos = max(float(np.sum(y)), 1.0)
    return torch.tensor((len(y) - pos) / pos, device=dev)


ROWS = []
fo = os.path.join(a.out, f"{a.mode}_{a.name}_oof.csv.gz")


def save():
    pd.concat(ROWS).to_csv(fo, index=False)


def frame(variant, task, seed, pts, y, p, **extra):
    return pd.DataFrame({"mode": a.mode, "variant": variant, "exp": "1a", "task": task, "seed": seed, "patient": pts,
                         "y": y, "p": p, **extra})


# ---------------- regress (X3) ----------------
if a.mode == "regress":
    for tgt in tasks.REGRESSION:
        d = L[[tgt, "site"]].dropna()
        pts, y = d.index.to_numpy(), d[tgt].to_numpy(float)
        for s in SEEDS:
            P = np.zeros(len(d))
            for tr, te in site_folds(pts, s):
                mu, sd = y[tr].mean(), y[tr].std()
                yz = torch.tensor((y[tr] - mu) / sd, dtype=torch.float32, device=dev)
                torch.manual_seed(s)
                net = ABMIL(D).to(dev)
                train(net, list(pts[tr]), lambda n, i, x: F.mse_loss(n(x)[0].float(), yz[i]), s)
                P[te] = np.array(infer(net, list(pts[te]), s, fn=lambda o: o[0].float().item())) * sd + mu
            r = np.corrcoef(y, P)[0, 1]
            log(f"{tgt} seed {s}: n {len(y)}, r = {r:.3f}")
            ROWS.append(frame("abmil_mse", tgt, s, pts, y, P).assign(exp="reg"))
            save()

# ---------------- multitask / single (X4b) ----------------
BIN = ["ER_IHC", "PR_IHC", "HER2_STRICT"]
PAM = ["LumA", "LumB", "Her2", "Basal", "Normal"]
if a.mode in ("multitask", "single"):
    d = L[BIN + tasks.REGRESSION + ["PAM50", "site"]]
    d = d[d[BIN].notna().any(axis=1)]
    pts = d.index.to_numpy()
    YB = d[BIN].to_numpy(float)
    YP = d.PAM50.map({k: i for i, k in enumerate(PAM)}).to_numpy(float)
    YR = d[tasks.REGRESSION].to_numpy(float)
    for s in SEEDS:
        if a.mode == "multitask":
            PB, PPAM, PR = np.zeros((len(d), 3)), np.zeros((len(d), 5)), np.zeros((len(d), 3))
            for tr, te in site_folds(pts, s):
                mu, sd = np.nanmean(YR[tr], 0), np.nanstd(YR[tr], 0)
                yb = torch.tensor(YB[tr], dtype=torch.float32, device=dev)
                yp = torch.tensor(np.nan_to_num(YP[tr], nan=-1), dtype=torch.long, device=dev)
                yr = torch.tensor((YR[tr] - mu) / sd, dtype=torch.float32, device=dev)
                pw = torch.stack([bce_weight(YB[tr][~np.isnan(YB[tr][:, j]), j]) for j in range(3)])

                def step_loss(n, i, x):
                    logit, aux = n(x)
                    loss = 0.0
                    for j in range(3):
                        if not torch.isnan(yb[i, j]):
                            loss = loss + F.binary_cross_entropy_with_logits(logit[j].float(), yb[i, j], pos_weight=pw[j])
                    if yp[i] >= 0:
                        loss = loss + F.cross_entropy(aux["pam"].float()[None], yp[i][None])
                    for j in range(3):
                        if not torch.isnan(yr[i, j]):
                            loss = loss + F.mse_loss(aux["reg"][j].float(), yr[i, j])
                    return loss

                torch.manual_seed(s)
                net = train(MultiTaskABMIL(D, 3, 5, 3).to(dev), list(pts[tr]), step_loss, s)
                o = infer(net, list(pts[te]), s, fn=lambda o: (torch.sigmoid(o[0].float()).cpu().numpy(),
                                                              torch.softmax(o[1]["pam"].float(), 0).cpu().numpy(),
                                                              o[1]["reg"].float().cpu().numpy()))
                PB[te], PPAM[te] = np.stack([q[0] for q in o]), np.stack([q[1] for q in o])
                PR[te] = np.stack([q[2] for q in o]) * sd + mu
            for j, c in enumerate(BIN):
                m = ~np.isnan(YB[:, j])
                ROWS.append(frame("multitask", c, s, pts[m], YB[m, j], PB[m, j]))
            for j, c in enumerate(tasks.REGRESSION):
                m = ~np.isnan(YR[:, j])
                ROWS.append(frame("multitask", c, s, pts[m], YR[m, j], PR[m, j]).assign(exp="reg"))
            m = ~np.isnan(YP)
            ROWS.append(frame("multitask", "PAM50", s, pts[m], YP[m], PPAM[m].argmax(1),
                              **{f"p_{k}": PPAM[m, i] for i, k in enumerate(PAM)}).assign(exp="pam50"))
        else:
            for j, c in enumerate(BIN):
                m = ~np.isnan(YB[:, j])
                pj, yj = pts[m], YB[m, j]
                P = np.zeros(len(pj))
                for tr, te in site_folds(pj, s):
                    lossf = nn.BCEWithLogitsLoss(pos_weight=bce_weight(yj[tr]))
                    yt = torch.tensor(yj[tr], dtype=torch.float32, device=dev)
                    torch.manual_seed(s)
                    net = train(ABMIL(D).to(dev), list(pj[tr]), lambda n, i, x: lossf(n(x)[0].float(), yt[i]), s)
                    P[te] = infer(net, list(pj[te]), s)
                ROWS.append(frame("single", c, s, pj, yj, P))
        log(f"{a.mode} seed {s} done")
        save()

# ---------------- dann (X4c ii) ----------------
if a.mode == "dann":
    EMB = {}
    for group, cols in (("ER", ["ER_IHC", "ESR1_RNApos"]), ("HER2", ["HER2_STRICT"])):
        req = ["ER_IHC", "ESR1_RNApos", "ER_RPPApos"] if group == "ER" else ["HER2_STRICT", "ERBB2_AMP"]
        d = L[req + ["site"]].dropna()
        pts, site = d.index.to_numpy(), d.site.to_numpy(dtype=str)
        for variant in ("plain", "dann"):
            for s in SEEDS:
                for c in cols:
                    y, P = d[c].to_numpy(float), np.zeros(len(d))
                    Z = np.zeros((len(d), 256), np.float32)
                    for tr, te in site_folds(pts, s):
                        smap = {k: i for i, k in enumerate(sorted(set(site[tr])))}
                        sy = torch.tensor([smap[k] for k in site[tr]], dtype=torch.long, device=dev)
                        lossf = nn.BCEWithLogitsLoss(pos_weight=bce_weight(y[tr]))
                        yt = torch.tensor(y[tr], dtype=torch.float32, device=dev)
                        torch.manual_seed(s)
                        net = ABMILDANN(D, len(smap), use_site=(variant == "dann")).to(dev)

                        def on_step(n, i, prog):
                            n.lam = 0.1 * (2 / (1 + np.exp(-10 * prog)) - 1)
                            n.site_y = sy[i]

                        train(net, list(pts[tr]), lambda n, i, x: n.objective(*n(x), yt[i], lossf), s, on_step)
                        o = infer(net, list(pts[te]), s, fn=lambda o: (torch.sigmoid(o[0].float()).item(),
                                                                      o[1]["z"].float().cpu().numpy()))
                        P[te], Z[te] = [q[0] for q in o], np.stack([q[1] for q in o])
                    ROWS.append(frame(variant, c, s, pts, y, P))
                    if s == 0:
                        EMB[f"{variant}|{c}"] = Z
                        EMB[f"patients|{c}"] = pts
                    log(f"{variant} {c} seed {s}: AUC {roc_auc_score(y, P):.3f}")
                    save()
    np.savez(os.path.join(a.out, f"dann_{a.name}_embeddings.npz"), **EMB)

# ---------------- noisy (X4d) ----------------
if a.mode == "noisy":
    d = L[["ER_IHC", "ESR1_RNApos", "ER_RPPApos", "site"]].dropna()
    pts = d.index.to_numpy()
    yi, yr, yp = (d[c].to_numpy(float) for c in ("ER_IHC", "ESR1_RNApos", "ER_RPPApos"))

    def gce(logit, y, w, q=0.7):  # generalised cross-entropy (Zhang & Sabuncu 2018), class-weighted
        p1 = torch.sigmoid(logit)
        py = y * p1 + (1 - y) * (1 - p1)
        return w * (1 - py.clamp_min(1e-6) ** q) / q

    for variant in ("bce", "consensus", "soft", "gce", "smooth"):
        for s in SEEDS:
            P = np.zeros(len(d))
            for tr, te in site_folds(pts, s):
                tr_v, target = tr, yi[tr]
                if variant == "consensus":
                    keep = yi[tr] == yr[tr]
                    tr_v, target = tr[keep], yi[tr][keep]
                elif variant == "soft":
                    target = dawid_skene(np.c_[yi[tr], yr[tr], yp[tr]])["posterior"]
                elif variant == "smooth":
                    target = yi[tr] * 0.9 + 0.05
                pwt = bce_weight((target >= 0.5).astype(float))
                yt = torch.tensor(target, dtype=torch.float32, device=dev)
                if variant == "gce":
                    def sl(n, i, x):
                        return gce(n(x)[0].float(), yt[i], torch.where(yt[i] > 0.5, pwt, torch.ones_like(pwt)))
                else:
                    lossf = nn.BCEWithLogitsLoss(pos_weight=pwt)

                    def sl(n, i, x):
                        return lossf(n(x)[0].float(), yt[i])
                torch.manual_seed(s)
                net = train(ABMIL(D).to(dev), list(pts[tr_v]), sl, s)
                P[te] = infer(net, list(pts[te]), s)
            ROWS.append(frame(variant, "ER_trained", s, pts, yi, P, y_rna=yr, y_rppa=yp))
            log(f"{variant} seed {s}: AUC vs IHC {roc_auc_score(yi, P):.3f}, vs RNA {roc_auc_score(yr, P):.3f}")
            save()

# ---------------- reliab (X6) ----------------
if a.mode == "reliab":
    for c in BIN:
        d = L[[c, "site"]].dropna()
        pts, y = d.index.to_numpy(), d[c].to_numpy(float)
        for s in SEEDS:
            cols = {k: np.zeros(len(d)) for k in ("p", "mc_mean", "mc_sd", "hb_mean", "hb_sd")}
            for tr, te in site_folds(pts, s):
                lossf = nn.BCEWithLogitsLoss(pos_weight=bce_weight(y[tr]))
                yt = torch.tensor(y[tr], dtype=torch.float32, device=dev)
                torch.manual_seed(s)
                net = train(ABMIL(D).to(dev), list(pts[tr]), lambda n, i, x: lossf(n(x)[0].float(), yt[i]), s)
                cols["p"][te] = infer(net, list(pts[te]), s)
                mc = np.array([infer(net, list(pts[te]), s + 100 * k, mode="mc") for k in range(a.passes)])
                hb = np.array([infer(net, list(pts[te]), s + 100 * k, mode="half") for k in range(a.passes)])
                cols["mc_mean"][te], cols["mc_sd"][te] = mc.mean(0), mc.std(0)
                cols["hb_mean"][te], cols["hb_sd"][te] = hb.mean(0), hb.std(0)
            p = cols.pop("p")
            ROWS.append(frame("abmil", c, s, pts, y, p, **cols))
            log(f"{c} seed {s} done; mean MC SD {cols['mc_sd'].mean():.3f}, half-bag SD {cols['hb_sd'].mean():.3f}")
            save()

save()
log("BENCH_MIL_EXT_OK")
