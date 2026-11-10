"""X2: aggregation methods beyond ABMIL on patch features (docs/DESIGN_SPEC.md §5.2).

Same experiment definitions as train_mil.py on the fixed site-grouped folds (common.site_folds), 3 seeds:
  1a  paired label groups: ER (IHC / RNA / RPPA), PR (IHC / prevalence-matched RNA), HER2 (strict / amplification)
  1d  HER2 1+ vs 0 (patients with ERBB2 mRNA, for the mRNA comparator)
  1b  train on all concordant patients, score the discordant ones (probabilities averaged over seeds)
Methods
  GPU  abmil, clam, transmil, dsmil, acmil (scripts/mil_models.py): Adam (lr 2e-4, weight decay 1e-4), 15 epochs (20 for
       TransMIL), class-weighted BCE, 2,048 random patches per training step, all (<= 4,096) at test, fp16 autocast.
       Every method uses its paper defaults with no tuning, as the ABMIL reference was run (same budget for all).
  CPU  --method pool runs the pooling baselines on the same patients and folds: mean (patient mean of all patches,
       the X1 / train_mil linear probe), meanstd and meanmax (mean of all patches + per-feature SD / max over the
       <= 4,096 bag patches), knn25 (cosine kNN on z-scored means, class-balanced vote), proto (nearest class centroid).
Outputs (--out): <method>_<name>_oof.csv.gz (method, exp, task, group, seed, patient, y, p; 1b rows have seed -1) and
<method>_<name>_runtime.csv (seconds per fit, parameters). Comparison tables: scripts/analyse_x2.py.
Usage: python scripts/bench_mil.py --method clam --name uni_v1_20x [--smoke]
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
from common import linear_probe, site_folds  # noqa: E402

warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser()
ap.add_argument("--method", required=True, choices=["abmil", "clam", "transmil", "dsmil", "acmil", "pool"])
ap.add_argument("--name", default="uni_v1_20x")
ap.add_argument("--feats", default="data/features")
ap.add_argument("--cap", type=int, default=4096)
ap.add_argument("--seeds", type=int, default=3)
ap.add_argument("--epochs", type=int, default=0, help="0 = method default (15; TransMIL 20)")
ap.add_argument("--lr", type=float, default=2e-4)
ap.add_argument("--train-patches", type=int, default=2048)
ap.add_argument("--smoke", action="store_true", help="ER group only, 1 seed, 1 epoch, 150 patients, ER 1b")
ap.add_argument("--out", default="results/bench/x2")
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
t0 = time.time()
GROUPS = [["ER_IHC", "ESR1_RNApos", "ER_RPPApos"], ["PR_IHC", "PGR_prevmatch"], ["HER2_STRICT", "ERBB2_AMP"],
          ["HER2low_1plus"]]
PAIRS = [("ER_IHC", "ESR1_RNApos"), ("PR_IHC", "PGR_prevmatch"), ("HER2_STRICT", "ERBB2_AMP")]
SEEDS = list(range(a.seeds))
EPOCHS = a.epochs or (20 if a.method == "transmil" else 15)
if a.smoke:
    GROUPS, PAIRS, SEEDS, EPOCHS = GROUPS[:1], PAIRS[:1], [0], 1


def log(*x):
    print(f"[{time.time() - t0:6.0f}s]", *x, flush=True)


# ---------------- labels and features (patient order = label file order, as train_mil.py) ----------------
idx = pd.read_csv(os.path.join(a.feats, f"{a.name}_index.csv")).set_index("patient")
L = tasks.labels()
L = L[L.index.isin(idx.index)]
if a.smoke:
    L = L[L.index.isin(np.random.default_rng(0).choice(L.index, 150, replace=False))]
z = np.load(os.path.join(a.feats, f"{a.name}_mean.npz"), allow_pickle=True)
MEAN = pd.DataFrame(z["X"], index=z["patients"].astype(str))
t = time.time()
BAGS = np.load(os.path.join(a.feats, f"{a.name}_cap{a.cap}.npy"), mmap_mode="r" if a.smoke or a.method == "pool" else None)
log(f"method {a.method}; {len(L)} patients with features and labels; bags {BAGS.shape} in {time.time() - t:.0f}s")

if a.method == "pool":
    from sklearn.preprocessing import StandardScaler

    rows = []
    for i, p in enumerate(L.index):
        x = np.asarray(BAGS[idx.at[p, "row"], :idx.at[p, "n_kept"]], dtype=np.float32)
        rows.append((x.std(0), x.max(0)))
    SD = pd.DataFrame(np.stack([r[0] for r in rows]), index=L.index)
    MX = pd.DataFrame(np.stack([r[1] for r in rows]), index=L.index)
    M = MEAN.loc[L.index]
    FEAT = {"mean": M, "meanstd": pd.concat([M, SD.add_prefix("s")], axis=1),
            "meanmax": pd.concat([M, MX.add_prefix("m")], axis=1)}
    for k in FEAT:
        FEAT[k].columns = [str(c) for c in FEAT[k].columns]
    log("pooled features built")

    def fit_pool(sub, tr, te, y_tr):
        if sub in FEAT:
            return linear_probe().fit(FEAT[sub].loc[tr].values, y_tr).predict_proba(FEAT[sub].loc[te].values)[:, 1]
        sc = StandardScaler().fit(M.loc[tr].values)
        Xtr, Xte = sc.transform(M.loc[tr].values), sc.transform(M.loc[te].values)
        Xtr /= np.linalg.norm(Xtr, axis=1, keepdims=True)
        Xte /= np.linalg.norm(Xte, axis=1, keepdims=True)
        S = Xte @ Xtr.T
        if sub == "proto":  # cosine similarity to the class centroids, mapped to (0, 1)
            c1, c0 = Xtr[y_tr == 1].mean(0), Xtr[y_tr == 0].mean(0)
            d = Xte @ (c1 / np.linalg.norm(c1)) - Xte @ (c0 / np.linalg.norm(c0))
            return 1 / (1 + np.exp(-20 * d))
        nb = np.argsort(-S, axis=1)[:, :25]  # knn25, class-balanced vote (0.5 = no preference)
        f1, f0 = (y_tr[nb] == 1).mean(1) / y_tr.mean(), (y_tr[nb] == 0).mean(1) / (1 - y_tr.mean())
        return f1 / (f1 + f0 + 1e-12)

    SUBS = ["mean", "meanstd", "meanmax", "knn25", "proto"]
    out, rt = [], []
    for sub in SUBS:
        for group in GROUPS:
            req = ["ERBB2_log2"] if group == ["HER2low_1plus"] else []
            d = L[group + req + ["site"]].dropna()
            pts = d.index.to_numpy()
            for s in SEEDS:
                for col in group:
                    y, p = d[col].astype(int).to_numpy(), np.zeros(len(d))
                    for tr, te in site_folds(pts, s):
                        p[te] = fit_pool(sub, pts[tr], pts[te], y[tr])
                    out.append(pd.DataFrame({"method": sub, "exp": "1a", "task": col, "group": "|".join(group),
                                             "seed": s, "patient": pts, "y": y, "p": p}))
        for ihc, mol in PAIRS:
            d = L[[ihc, mol]].dropna()
            conc, disc = d[d[ihc] == d[mol]], d[d[ihc] != d[mol]]
            p = fit_pool(sub, conc.index.to_numpy(), disc.index.to_numpy(), conc[ihc].astype(int).to_numpy())
            out.append(pd.DataFrame({"method": sub, "exp": "1b", "task": f"{ihc} vs {mol}", "group": "", "seed": -1,
                                     "patient": disc.index, "y": disc[mol].values, "p": p}))
        log(f"{sub} done")
    pd.concat(out).to_csv(os.path.join(a.out, f"pool_{a.name}_oof.csv.gz"), index=False)
    log("BENCH_MIL_OK")
    sys.exit(0)

# ---------------- MIL heads (GPU) ----------------
import torch  # noqa: E402
import torch.nn as nn  # noqa: E402

from mil_models import MODELS  # noqa: E402

dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
log("device", dev, torch.cuda.get_device_name(0) if dev.type == "cuda" else "")
if dev.type != "cuda" and os.environ.get("SLURM_JOB_ID"):
    sys.exit(f"no usable GPU on {os.environ.get('SLURMD_NODENAME')}: refusing to train on CPU")
AMP = dev.type == "cuda"


def bag(p, train, rng):
    r, k = idx.at[p, "row"], idx.at[p, "n_kept"]
    x = BAGS[r, :k]
    if train and k > a.train_patches:
        x = x[np.sort(rng.choice(k, a.train_patches, replace=False))] if a.method == "transmil" else \
            x[rng.choice(k, a.train_patches, replace=False)]
    return torch.from_numpy(np.ascontiguousarray(x)).to(dev, non_blocking=True)


RT = []


def fit_predict(tr, te, y_tr, seed):
    """tr / te: lists of patient ids; returns probabilities for te (same loop as train_mil.py)"""
    t = time.time()
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    net = MODELS[a.method](BAGS.shape[2]).to(dev)
    opt = torch.optim.Adam(net.parameters(), lr=a.lr, weight_decay=1e-4)
    pos = max(y_tr.sum(), 1)
    lossf = nn.BCEWithLogitsLoss(pos_weight=torch.tensor((len(y_tr) - pos) / pos, device=dev))
    scaler = torch.cuda.amp.GradScaler(enabled=AMP)
    yt = torch.tensor(y_tr, dtype=torch.float32, device=dev)
    for _ in range(EPOCHS):
        net.train()
        for i in rng.permutation(len(tr)):
            with torch.autocast(device_type=dev.type, dtype=torch.float16, enabled=AMP):
                logit, aux = net(bag(tr[i], True, rng).float())
                loss = net.objective(logit.float(), aux, yt[i], lossf)
            opt.zero_grad(set_to_none=True)
            scaler.scale(loss).backward()
            scaler.step(opt)
            scaler.update()
    net.eval()
    out = []
    with torch.no_grad(), torch.autocast(device_type=dev.type, dtype=torch.float16, enabled=AMP):
        for p in te:
            out.append(torch.sigmoid(net(bag(p, False, rng).float())[0].float()).item())
    RT.append({"seconds": time.time() - t, "n_train": len(tr), "params": sum(q.numel() for q in net.parameters())})
    return np.array(out)


out = []
fo = os.path.join(a.out, f"{a.method}_{a.name}_oof.csv.gz")
for group in GROUPS:
    req = ["ERBB2_log2"] if group == ["HER2low_1plus"] else []
    d = L[group + req + ["site"]].dropna()
    pts = d.index.to_numpy()
    for s in SEEDS:
        P = {col: np.zeros(len(d)) for col in group}
        for tr, te in site_folds(pts, s):
            for col in group:
                P[col][te] = fit_predict(list(pts[tr]), list(pts[te]), d[col].astype(int).values[tr], s)
        for col in group:
            out.append(pd.DataFrame({"method": a.method, "exp": "1a", "task": col, "group": "|".join(group), "seed": s,
                                     "patient": pts, "y": d[col].values, "p": P[col]}))
        from sklearn.metrics import roc_auc_score
        log(f"{group} seed {s}: " + ", ".join(f"{c}={roc_auc_score(d[c].astype(int), P[c]):.3f}" for c in group))
        pd.concat(out).to_csv(fo, index=False)  # checkpoint after every seed
for ihc, mol in PAIRS:
    d = L[[ihc, mol]].dropna()
    conc, disc = d[d[ihc] == d[mol]], d[d[ihc] != d[mol]]
    p = np.mean([fit_predict(list(conc.index), list(disc.index), conc[ihc].astype(int).values, s) for s in SEEDS], 0)
    k = int(((p >= 0.5).astype(int) == disc[mol].astype(int).values).sum())
    log(f"1b {ihc} vs {mol}: n = {len(disc)}, molecular {k} / IHC {len(disc) - k}")
    out.append(pd.DataFrame({"method": a.method, "exp": "1b", "task": f"{ihc} vs {mol}", "group": "", "seed": -1,
                             "patient": disc.index, "y": disc[mol].values, "p": p}))
pd.concat(out).to_csv(fo, index=False)
pd.DataFrame(RT).assign(method=a.method).to_csv(os.path.join(a.out, f"{a.method}_{a.name}_runtime.csv"), index=False)
log(f"mean seconds per fit {np.mean([r['seconds'] for r in RT]):.0f}; parameters {RT[0]['params']:,}")
log("BENCH_MIL_OK")
