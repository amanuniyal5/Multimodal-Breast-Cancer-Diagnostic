"""X5b patch-level maps and map metrics (docs/DESIGN_SPEC.md §5.5) for the map patients (data/metadata/maps_patients.csv).

For each patient, models are trained WITHOUT that patient's hospital: the fixed seed-0 site-grouped fold that holds
the patient is the test fold; training = all other ER_IHC patients (ER = the label of interest).
  contribution  exact linear map of the mean-pool linear probe (StandardScaler + LR): logit = b + (1/N) sum_i c_i with
                c_i = sum_j (w_j / s_j)(x_ij - m_j); red = towards ER+, blue = towards ER-. Identity check (S19): the
                mean of c_i plus b equals the probe's logit on the mean of the same patches (|diff| < 1e-4).
  attention     ABMIL, CLAM-SB, DSMIL, ACMIL (scripts/mil_models.py, X2 settings, seed 0) on UNI bags, applied to ALL
                patches of the patient.
  PCA-RGB       first 3 principal components (fitted on 200k random cohort patches per encoder) -> RGB:
                UNI, CONCH (all patches) and ResNet-50 (its <= 4,096 sampled patches, subset patients only).
  concepts      CONCH zero-shot concept shares per patch (14 concepts of conch_zeroshot.py, projected CONCH space);
                maps for tubules, pleomorphism, necrosis, apocrine; each UNI patch takes the concept of the CONCH patch
                (512 px) that contains its centre.
Map metrics (results/bench/x5/map_metrics.csv): Spearman correlation of each attention map with the contribution map;
contribution share from tumour-concept patches vs stroma / adipose patches (and their area shares); deletion curve:
logit after removing the top k % patches pushing towards the predicted class vs k % random patches (exact for the
linear probe); smallest k % that flips the call.
Outputs: results/bench/x5/maps/<patient>.png (overview), results/bench/x5/maps/panels/<patient>_<map>.png (single maps
for the app), results/bench/x5/maps_npz/<patient>.npz (per-patch values), map_metrics.csv.
Usage (Ada, GPU): python scripts/maps.py --work /scratch/$USER/maps
"""
import argparse
import os
import sys
import time
import warnings
from concurrent.futures import ThreadPoolExecutor

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import requests  # noqa: E402
import torch  # noqa: E402
import torch.nn as nn  # noqa: E402
import torch.nn.functional as F  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bench import tasks  # noqa: E402
from common import linear_probe, site_folds  # noqa: E402
from mil_models import MODELS  # noqa: E402

warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser()
ap.add_argument("--patients", default="data/metadata/maps_patients.csv")
ap.add_argument("--feats", default="data/features")
ap.add_argument("--work", required=True, help="scratch dir for slide downloads")
ap.add_argument("--idc-objects", default="data/metadata/idc_slides.csv")
ap.add_argument("--out", default="results/bench/x5")
ap.add_argument("--methods", default="abmil,clam,dsmil,acmil")
ap.add_argument("--epochs", type=int, default=15)
ap.add_argument("--thumb", type=int, default=1400)
ap.add_argument("--no-slides", action="store_true", help="skip IDC thumbnails (maps on a blank canvas)")
a = ap.parse_args()
for d in ("maps", "maps/panels", "maps_npz"):
    os.makedirs(os.path.join(a.out, d), exist_ok=True)
t0 = time.time()
dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
AMP = dev.type == "cuda"


def log(*x):
    print(f"[{time.time() - t0:6.0f}s]", *x, flush=True)


TARGETS = pd.read_csv(a.patients).patient.tolist()
L = tasks.labels()
idx = pd.read_csv(os.path.join(a.feats, "uni_v1_20x_index.csv")).set_index("patient")
z = np.load(os.path.join(a.feats, "uni_v1_20x_mean.npz"), allow_pickle=True)
MEAN = pd.DataFrame(z["X"], index=z["patients"].astype(str))
er = L[["ER_IHC"]].dropna()
pts = er.index.intersection(idx.index).intersection(MEAN.index).to_numpy()
y = er.loc[pts, "ER_IHC"].astype(int).to_numpy()
fold = np.zeros(len(pts), int)
for k, (tr, te) in enumerate(site_folds(pts, 0)):
    fold[te] = k
FOLD = dict(zip(pts, fold))
TARGETS = [p for p in TARGETS if p in FOLD and os.path.exists(os.path.join(a.feats, "maps_uni_v1_20x", f"{p}.npz"))]
log(f"{len(TARGETS)} target patients; training pool {len(pts)} ER_IHC patients; device {dev}")

# ---------------- full patch sets ----------------
U = {p: np.load(os.path.join(a.feats, "maps_uni_v1_20x", f"{p}.npz"), allow_pickle=True) for p in TARGETS}
Cn = {p: (np.load(os.path.join(a.feats, "maps_conch_v1_20x", f"{p}.npz"), allow_pickle=True)
          if os.path.exists(os.path.join(a.feats, "maps_conch_v1_20x", f"{p}.npz")) else None) for p in TARGETS}

# ---------------- 1. contribution maps (linear probe per fold) ----------------
RES = {p: {} for p in TARGETS}
S19 = []
for k in sorted({FOLD[p] for p in TARGETS}):
    tr = pts[fold != k]
    m = linear_probe().fit(MEAN.loc[tr].values, y[fold != k])
    sc, lr = m[0], m[-1]
    w = lr.coef_.ravel() / sc.scale_
    b = lr.intercept_[0]
    for p in [p for p in TARGETS if FOLD[p] == k]:
        X = U[p]["X"].astype(np.float32)
        c = (X - sc.mean_) @ w
        logit_probe = m.decision_function(X.mean(0, keepdims=True))[0]
        S19.append(abs(b + c.mean() - logit_probe))
        RES[p]["contribution"] = c
        RES[p]["probe_logit"] = b + c.mean()
        RES[p]["probe_p"] = 1 / (1 + np.exp(-(b + c.mean())))
        RES[p]["probe_b"] = b
log(f"S19 contribution identity: max |b + mean(c) - logit| = {max(S19):.2e} ({'PASS' if max(S19) < 1e-4 else 'FAIL'})")

# ---------------- 2. attention maps (MIL heads per fold) ----------------
t = time.time()
BAGS = np.load(os.path.join(a.feats, "uni_v1_20x_cap4096.npy"))
log(f"UNI bags {BAGS.shape} loaded in {time.time() - t:.0f}s")


def bag(p, rng, n=2048):
    r, kk = idx.at[p, "row"], idx.at[p, "n_kept"]
    x = BAGS[r, :kk]
    if kk > n:
        x = x[rng.choice(kk, n, replace=False)]
    return torch.from_numpy(np.ascontiguousarray(x)).to(dev).float()


def attention_of(method, net, x):
    with torch.no_grad(), torch.autocast(device_type=dev.type, dtype=torch.float16, enabled=AMP):
        logit, aux = net(x)
    return torch.sigmoid(logit.float()).item(), aux["attention"].float().cpu().numpy()


for method in a.methods.split(","):
    for k in sorted({FOLD[p] for p in TARGETS}):
        tr = list(pts[fold != k])
        ytr = y[fold != k]
        torch.manual_seed(0)
        rng = np.random.default_rng(0)
        net = MODELS[method](BAGS.shape[2]).to(dev)
        opt = torch.optim.Adam(net.parameters(), lr=2e-4, weight_decay=1e-4)
        pos = max(ytr.sum(), 1)
        lossf = nn.BCEWithLogitsLoss(pos_weight=torch.tensor((len(ytr) - pos) / pos, device=dev))
        scaler = torch.cuda.amp.GradScaler(enabled=AMP)
        yt = torch.tensor(ytr, dtype=torch.float32, device=dev)
        for _ in range(a.epochs):
            net.train()
            for i in rng.permutation(len(tr)):
                with torch.autocast(device_type=dev.type, dtype=torch.float16, enabled=AMP):
                    logit, aux = net(bag(tr[i], rng))
                    loss = net.objective(logit.float(), aux, yt[i], lossf)
                opt.zero_grad(set_to_none=True)
                scaler.scale(loss).backward()
                scaler.step(opt)
                scaler.update()
        net.eval()
        for p in [p for p in TARGETS if FOLD[p] == k]:
            prob, att = attention_of(method, net, torch.from_numpy(U[p]["X"]).to(dev).float())
            RES[p][f"att_{method}"], RES[p][f"p_{method}"] = att, prob
        log(f"{method} fold {k} trained and applied")
del BAGS

# ---------------- 3. PCA-RGB ----------------
def pca_fit(path, n=200_000, seed=0):
    B = np.load(path, mmap_mode="r")
    r = np.random.default_rng(seed)
    rows = r.integers(0, B.shape[0], n)
    cols = r.integers(0, B.shape[1], n)
    S = np.asarray(B[rows, cols], dtype=np.float32)
    S = S[np.abs(S).sum(1) > 0]  # padded rows are zero
    mu = S.mean(0)
    _, _, Vt = np.linalg.svd(S[:50_000] - mu, full_matrices=False)
    proj = (S - mu) @ Vt[:3].T
    lo, hi = np.percentile(proj, 1, 0), np.percentile(proj, 99, 0)
    return mu, Vt[:3], lo, hi


def pca_rgb(X, P):
    mu, V, lo, hi = P
    return np.clip(((X.astype(np.float32) - mu) @ V.T - lo) / (hi - lo), 0, 1)


P_UNI = pca_fit(os.path.join(a.feats, "uni_v1_20x_cap4096.npy"))
for p in TARGETS:
    RES[p]["pca_uni"] = pca_rgb(U[p]["X"], P_UNI)
if os.path.exists(os.path.join(a.feats, "conch_v1_20x_cap4096.npy")):
    P_CONCH = pca_fit(os.path.join(a.feats, "conch_v1_20x_cap4096.npy"))
    for p in TARGETS:
        if Cn[p] is not None:
            RES[p]["pca_conch"] = pca_rgb(Cn[p]["X"], P_CONCH)
rn_dir = os.path.join(a.feats, "dx_resnet50")
if os.path.isdir(rn_dir):
    import glob
    files = sorted(glob.glob(os.path.join(rn_dir, "*.npz")))
    r = np.random.default_rng(0)
    samp = []
    for f in r.choice(files, min(200, len(files)), replace=False):
        bz = np.load(f)["bag"]
        samp.append(bz[r.choice(len(bz), min(1000, len(bz)), replace=False)].astype(np.float32))
    S = np.concatenate(samp)
    mu = S.mean(0)
    _, _, Vt = np.linalg.svd(S[:50_000] - mu, full_matrices=False)
    pr = (S - mu) @ Vt[:3].T
    P_RN = (mu, Vt[:3], np.percentile(pr, 1, 0), np.percentile(pr, 99, 0))
    for p in TARGETS:
        fs = [f for f in files if os.path.basename(f).startswith(p)]
        if fs:
            zz = np.load(fs[0])
            RES[p]["resnet"] = (zz["coords"], pca_rgb(zz["bag"], P_RN), os.path.basename(fs[0]).split(".")[0])
log("PCA-RGB done")

# ---------------- 4. CONCH concept maps ----------------
CONCEPTS = ["tubules", "low_grade_nuclei", "pleomorphism", "mitoses", "necrosis", "lymphocytes", "apocrine",
            "solid_sheets", "lobular", "mucin", "stroma", "adipose", "dcis", "normal_breast"]
TUMOUR = {"tubules", "low_grade_nuclei", "pleomorphism", "mitoses", "solid_sheets", "lobular", "apocrine", "mucin"}
STROMA = {"stroma", "adipose"}
try:
    from conch.open_clip_custom import create_model_from_pretrained, get_tokenizer, tokenize

    cm, _ = create_model_from_pretrained("conch_ViT-B-16", "hf_hub:MahmoodLab/conch")
    cm = cm.to(dev).eval()
    import importlib.util
    src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "conch_zeroshot.py")).read()
    ns = {}
    exec(src[src.index("TEMPLATES = ["):src.index("tok = get_tokenizer()")], ns)  # same prompts as conch_zeroshot.py
    TEMPLATES, CDEF = ns["TEMPLATES"], ns["CONCEPTS"]
    tok = get_tokenizer()
    with torch.no_grad():
        g = torch.Generator().manual_seed(0)
        xx = torch.randn(4, 3, 448, 448, generator=g).to(dev)
        pre = cm.encode_image(xx, proj_contrast=False, normalize=False)
        post = cm.encode_image(xx, proj_contrast=True, normalize=False)
        proj = None
        for name, mod in cm.visual.named_modules():
            if name and isinstance(mod, (nn.Linear, nn.Sequential)):
                try:
                    if torch.allclose(mod(pre), post, atol=1e-3, rtol=1e-3):
                        proj = mod
                        break
                except Exception:  # noqa: BLE001
                    pass
        if proj is None:
            for name, prm in cm.visual.named_parameters():
                if prm.ndim == 2 and prm.shape[0] == pre.shape[1] and torch.allclose(pre @ prm, post, atol=1e-3, rtol=1e-3):
                    proj = (lambda q, prm=prm: q @ prm)
                    break

        def temb(names):
            pr_ = [t.replace("CLASSNAME", n) for n in names for t in TEMPLATES]
            e = F.normalize(cm.encode_text(tokenize(texts=pr_, tokenizer=tok).to(dev)).float(), dim=-1).mean(0)
            return F.normalize(e, dim=-1)

        T = torch.stack([temb(CDEF[c]) for c in CONCEPTS])
        scale = cm.logit_scale.exp().item()
        for p in TARGETS:
            if Cn[p] is None:
                continue
            xc = torch.from_numpy(Cn[p]["X"]).to(dev).float()
            sh = torch.cat([torch.softmax(scale * F.normalize(proj(xc[i:i + 8192]), dim=-1) @ T.T, 1)
                            for i in range(0, len(xc), 8192)]).cpu().numpy()
            RES[p]["concepts_conch"] = sh
            # map every UNI patch to the CONCH patch containing its centre (same slide)
            uc, us = U[p]["coords"].astype(float), U[p]["slide"]
            cc, cs = Cn[p]["coords"].astype(float), Cn[p]["slide"]
            un, cn_ = [str(s) for s in U[p]["slides"]], [str(s) for s in Cn[p]["slides"]]
            out = np.full((len(uc), len(CONCEPTS)), np.nan, np.float32)
            for j, name in enumerate(un):
                if name not in cn_:
                    continue
                mu_, mc_ = us == j, cs == cn_.index(name)
                if not mc_.any():
                    continue
                stepc = np.median(np.diff(np.unique(cc[mc_, 0]))) if mc_.sum() > 1 else 1024
                stepu = np.median(np.diff(np.unique(uc[mu_, 0]))) if mu_.sum() > 1 else 512
                tree = cKDTree(cc[mc_] + stepc / 2)
                dist, nn_ = tree.query(uc[mu_] + stepu / 2)
                ok = dist <= stepc * 0.75
                tmp = np.full((mu_.sum(), len(CONCEPTS)), np.nan, np.float32)
                tmp[ok] = sh[np.where(mc_)[0][nn_[ok]]]
                out[mu_] = tmp
            RES[p]["concepts_uni"] = out
    log("CONCH concept maps done")
except Exception as e:  # noqa: BLE001
    log(f"CONCH concept maps skipped: {e}")

# ---------------- 5. map metrics ----------------
rows = []
for p in TARGETS:
    R = RES[p]
    c = R["contribution"]
    lg = R["probe_logit"]
    sgn = 1 if lg >= 0 else -1
    r = {"patient": p, "n_patches": len(c), "ER_IHC": L.at[p, "ER_IHC"], "ESR1_RNApos": L.at[p, "ESR1_RNApos"],
         "PAM50": L.at[p, "PAM50"], "probe_p": R["probe_p"]}
    for m in a.methods.split(","):
        if f"att_{m}" in R:
            r[f"p_{m}"] = R[f"p_{m}"]
            r[f"spearman_att_{m}_vs_contrib"] = spearmanr(R[f"att_{m}"], c)[0]
            r[f"spearman_att_{m}_vs_abs_contrib"] = spearmanr(R[f"att_{m}"], np.abs(c))[0]
    if "concepts_uni" in R:
        cu = R["concepts_uni"]
        okc = ~np.isnan(cu).any(1)
        top = np.array([CONCEPTS[i] for i in np.nanargmax(np.where(okc[:, None], cu, -1), 1)])
        tum, stro = okc & np.isin(top, list(TUMOUR)), okc & np.isin(top, list(STROMA))
        tot = np.abs(c[okc]).sum()
        r.update({"area_tumour": tum.sum() / okc.sum(), "area_stroma_fat": stro.sum() / okc.sum(),
                  "abs_contrib_share_tumour": np.abs(c[tum]).sum() / tot,
                  "abs_contrib_share_stroma_fat": np.abs(c[stro]).sum() / tot,
                  "mean_contrib_tumour": c[tum].mean() if tum.any() else np.nan,
                  "mean_contrib_stroma_fat": c[stro].mean() if stro.any() else np.nan,
                  "patches_with_concept": okc.mean()})
        for m in ("abmil",):
            if f"att_{m}" in R:
                at = R[f"att_{m}"]
                r[f"att_{m}_share_tumour"] = at[tum].sum() / at[okc].sum()
                r[f"att_{m}_share_stroma_fat"] = at[stro].sum() / at[okc].sum()
    order = np.argsort(-sgn * c)  # patches pushing towards the predicted class first
    rng = np.random.default_rng(0)
    for kpct in (1, 5, 10, 20):
        kk = int(round(len(c) * kpct / 100))
        keep = np.ones(len(c), bool)
        keep[order[:kk]] = False
        r[f"logit_del_top{kpct}"] = R["probe_b"] + c[keep].mean()
        rnd = []
        for _ in range(20):
            kr = np.ones(len(c), bool)
            kr[rng.choice(len(c), kk, replace=False)] = False
            rnd.append(R["probe_b"] + c[kr].mean())
        r[f"logit_del_rand{kpct}"] = np.mean(rnd)
    flip = np.nan
    cs = np.cumsum(c[order])  # sum of the removed contributions (removed in order of pushing towards the call)
    for kk in range(1, len(c)):
        if sgn * (R["probe_b"] + (c.sum() - cs[kk - 1]) / (len(c) - kk)) < 0:
            flip = 100 * kk / len(c)
            break
    r["logit"] = lg
    r["pct_removed_to_flip"] = flip
    rows.append(r)
    np.savez_compressed(os.path.join(a.out, "maps_npz", f"{p}.npz"), coords=U[p]["coords"], slide=U[p]["slide"],
                        slides=U[p]["slides"], contribution=c.astype(np.float32),
                        **{k: v.astype(np.float32) for k, v in R.items() if k.startswith("att_")},
                        pca_uni=R["pca_uni"].astype(np.float16),
                        **({"concepts": R["concepts_uni"].astype(np.float16)} if "concepts_uni" in R else {}))
M = pd.DataFrame(rows)
M.to_csv(os.path.join(a.out, "map_metrics.csv"), index=False)
pd.set_option("display.width", 250)
log("map metrics:\n" + M.round(3).to_string(index=False))

# ---------------- 6. rendering ----------------
IDC = pd.read_csv(a.idc_objects)
GCS = "https://storage.googleapis.com/idc-open-data/"


def fetch(url, path, size, workers=16, step=16 << 20):
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC)
    os.ftruncate(fd, size)

    def rng_(b0):
        e = min(b0 + step, size) - 1
        for attempt in range(6):
            try:
                rr = requests.get(url, headers={"Range": f"bytes={b0}-{e}"}, timeout=(30, 120))
                if rr.status_code == 206 and len(rr.content) == e - b0 + 1:
                    os.pwrite(fd, rr.content, b0)
                    return
            except requests.RequestException:
                pass
            time.sleep(5 * 2 ** attempt)
        raise IOError(f"range {b0}-{e} failed")

    try:
        with ThreadPoolExecutor(workers) as ex:
            list(ex.map(rng_, range(0, size, step)))
    finally:
        os.close(fd)


def thumbnail(barcode):
    import openslide
    objs = IDC[IDC.slide_barcode == barcode]
    if objs.empty:
        return None, None
    wdir = os.path.join(a.work, barcode)
    os.makedirs(wdir, exist_ok=True)
    for o in objs.itertuples():
        fp = os.path.join(wdir, os.path.basename(o.object))
        if not os.path.exists(fp) or os.path.getsize(fp) != int(o.size):
            fetch(GCS + o.object, fp, int(o.size))
    sl = openslide.OpenSlide(os.path.join(wdir, os.path.basename(objs.sort_values("size").object.iloc[-1])))
    th = sl.get_thumbnail((a.thumb, a.thumb)).convert("RGB")
    return np.asarray(th), th.width / sl.dimensions[0]


def paint(shape, xy, ext, vals, scale_):
    H, W = shape[:2]
    out = np.full((H, W) + (vals.shape[1:] if vals.ndim > 1 else ()), np.nan, np.float32)
    for (x, y_), v in zip(xy, vals):
        x0, y0 = int(x * scale_), int(y_ * scale_)
        out[y0:max(y0 + 1, int((y_ + ext) * scale_)), x0:max(x0 + 1, int((x + ext) * scale_))] = v
    return out


for p in TARGETS:
    R, u = RES[p], U[p]
    names = [str(s) for s in u["slides"]]
    j = int(np.bincount(u["slide"]).argmax())  # the slide with most patches
    barcode = names[j].split(".")[0]
    m = u["slide"] == j
    xy = u["coords"][m]
    step = np.median(np.diff(np.unique(xy[:, 0]))) if m.sum() > 1 else 512
    th, scl = (None, None) if a.no_slides else thumbnail(barcode)
    if th is None:
        W, H = xy.max(0) + step
        scl = a.thumb / max(W, H)
        th = np.full((int(H * scl) + 1, int(W * scl) + 1, 3), 255, np.uint8)
    pct = lambda v: v.argsort().argsort() / max(len(v) - 1, 1)  # noqa: E731
    c = R["contribution"][m]
    lim = np.percentile(np.abs(c), 98)
    panels = [("H&E", None, None),
              ("Contribution (probe; red = ER+)", paint(th.shape, xy, step, np.clip(c / lim, -1, 1), scl), "bwr")]
    for meth in a.methods.split(","):
        if f"att_{meth}" in R:
            panels.append((f"{meth.upper()} attention (percentile)", paint(th.shape, xy, step, pct(R[f"att_{meth}"][m]),
                                                                            scl), "inferno"))
    panels.append(("UNI PCA-RGB", paint(th.shape, xy, step, R["pca_uni"][m], scl), None))
    if "concepts_uni" in R:
        for cn in ("tubules", "pleomorphism", "necrosis", "apocrine"):
            panels.append((f"CONCH concept: {cn}", paint(th.shape, xy, step, R["concepts_uni"][m, CONCEPTS.index(cn)],
                                                          scl), "viridis"))
    if "resnet" in R and R["resnet"][2] == barcode:
        rc, rv, _ = R["resnet"]
        panels.append(("ResNet-50 PCA-RGB (sampled patches)", paint(th.shape, rc, step, rv, scl), None))
    n = len(panels)
    cols = 4
    fig, axs = plt.subplots(int(np.ceil(n / cols)), cols, figsize=(4.2 * cols, 3.6 * int(np.ceil(n / cols))))
    for ax in axs.ravel():
        ax.axis("off")
    for ax, (title, img, cmap) in zip(axs.ravel(), panels):
        ax.imshow(th)
        if img is not None:
            if img.ndim == 3:
                rgba = np.dstack([np.nan_to_num(img), (~np.isnan(img[..., 0])).astype(float) * 0.9])
                ax.imshow(rgba)
            else:
                vmin, vmax = (-1, 1) if cmap == "bwr" else (0, 1) if "attention" in title else (np.nanmin(img), np.nanmax(img))
                ax.imshow(np.ma.masked_invalid(img), cmap=cmap, alpha=0.75, vmin=vmin, vmax=vmax)
        ax.set_title(title, fontsize=9)
        sp = plt.figure(figsize=(5, 5 * th.shape[0] / th.shape[1]))
        sa = sp.add_axes([0, 0, 1, 1])
        sa.axis("off")
        sa.imshow(th)
        if img is not None:
            if img.ndim == 3:
                sa.imshow(np.dstack([np.nan_to_num(img), (~np.isnan(img[..., 0])).astype(float) * 0.9]))
            else:
                vmin, vmax = (-1, 1) if cmap == "bwr" else (0, 1) if "attention" in title else (np.nanmin(img), np.nanmax(img))
                sa.imshow(np.ma.masked_invalid(img), cmap=cmap, alpha=0.75, vmin=vmin, vmax=vmax)
        key = title.split(" (")[0].lower().replace(" ", "_").replace(":", "").replace("&", "")
        sp.savefig(os.path.join(a.out, "maps", "panels", f"{p}_{key}.png"), dpi=110)
        plt.close(sp)
    row = M.set_index("patient").loc[p]
    fig.suptitle(f"{p} ({barcode}) · ER IHC {int(row.ER_IHC)} / RNA {row.ESR1_RNApos:.0f} · PAM50 {row.PAM50} · "
                 f"probe P(ER+) {row.probe_p:.2f} · ABMIL {row.get('p_abmil', np.nan):.2f}", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(a.out, "maps", f"{p}.png"), dpi=110)
    plt.close(fig)
    log(f"{p}: rendered {n} panels")
log("MAPS_OK")
