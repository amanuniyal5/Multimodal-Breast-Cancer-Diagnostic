"""ABMIL attention heatmaps for ER-discordant patients (incl. ER-low) and concordant reference patients.

1. Train ER ABMIL (same architecture/hyper-parameters as train_mil.py) on ER-concordant patients (IHC = RNA), 3 seeds,
   holding out 10 concordant reference patients (5 ER+, 5 ER-) so they are scored out-of-sample like the discordant ones.
2. Read the FULL UNI patch features + patch coordinates of the target patients from the TRACE parquet files.
3. For each target: per-patch attention (averaged over seeds), P(ER+) (averaged), which side it takes.
4. Write: results/attention/summary.csv, one PNG heatmap per patient (one panel per DX slide), an overview grid,
   and per-patient .npz (coords, slide index, attention) for later inspection.

Heatmap = patches drawn on their grid position, coloured by attention percentile within the patient (bright = high).

Usage (Ada, GPU):  python scripts/attention_maps.py --parquet-dir /scratch/$USER/trace/data/TCGA-BRCA/uni_v1/20x
"""
import argparse
import glob
import os
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402
import torch  # noqa: E402
import torch.nn as nn  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--feats", default="data/features")
ap.add_argument("--name", default="uni_v1_20x")
ap.add_argument("--cap", type=int, default=4096)
ap.add_argument("--labels", default="data/metadata/labels_patient.csv")
ap.add_argument("--parquet-dir", required=True)
ap.add_argument("--out", default="results/attention")
ap.add_argument("--seeds", type=int, default=3)
ap.add_argument("--epochs", type=int, default=15)
ap.add_argument("--train-patches", type=int, default=2048)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
t0 = time.time()


def log(*x):
    print(f"[{time.time() - t0:6.0f}s]", *x, flush=True)


dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
L = pd.read_csv(a.labels, index_col=0)
idx = pd.read_csv(os.path.join(a.feats, f"{a.name}_index.csv")).set_index("patient")
L = L[L.index.isin(idx.index)]
er = L[["ER_IHC", "ESR1_RNApos"]].dropna()
conc, disc = er[er.ER_IHC == er.ESR1_RNApos], er[er.ER_IHC != er.ESR1_RNApos]
rng0 = np.random.default_rng(0)
refs = list(rng0.choice(conc[conc.ER_IHC == 1].index, 5, replace=False)) + \
       list(rng0.choice(conc[conc.ER_IHC == 0].index, 5, replace=False))
train = conc.drop(index=refs)
targets = list(disc.index) + refs
log(f"train {len(train)} concordant, targets {len(disc)} discordant + {len(refs)} references")


class ABMIL(nn.Module):  # identical to train_mil.py, but forward can also return the attention weights
    def __init__(self, d, h=256, att=128, p=0.25):
        super().__init__()
        self.fc = nn.Sequential(nn.Linear(d, h), nn.ReLU(), nn.Dropout(p))
        self.a = nn.Sequential(nn.Linear(h, att), nn.Tanh())
        self.b = nn.Sequential(nn.Linear(h, att), nn.Sigmoid())
        self.c = nn.Linear(att, 1)
        self.cls = nn.Linear(h, 1)

    def forward(self, x, return_attention=False):
        hdn = self.fc(x)
        w = torch.softmax(self.c(self.a(hdn) * self.b(hdn)).squeeze(-1), dim=0)
        out = self.cls((w[:, None] * hdn).sum(0)).squeeze()
        return (out, w) if return_attention else out


# ---- 1. train on capped bags (same data as train_mil.py) ----
BAGS = np.load(os.path.join(a.feats, f"{a.name}_cap{a.cap}.npy"), mmap_mode="r")
tr_pts, y_tr = list(train.index), train.ER_IHC.astype(int).values
models = []
for s in range(a.seeds):
    torch.manual_seed(s)
    rng = np.random.default_rng(s)
    net = ABMIL(BAGS.shape[2]).to(dev)
    opt = torch.optim.Adam(net.parameters(), lr=2e-4, weight_decay=1e-4)
    pos = y_tr.sum()
    lossf = nn.BCEWithLogitsLoss(pos_weight=torch.tensor((len(y_tr) - pos) / pos, device=dev))
    scaler = torch.cuda.amp.GradScaler(enabled=dev.type == "cuda")
    yt = torch.tensor(y_tr, dtype=torch.float32, device=dev)
    for ep in range(a.epochs):
        net.train()
        for i in rng.permutation(len(tr_pts)):
            r, k = idx.at[tr_pts[i], "row"], idx.at[tr_pts[i], "n_kept"]
            x = np.asarray(BAGS[r, :k])
            if k > a.train_patches:
                x = x[rng.choice(k, a.train_patches, replace=False)]
            with torch.autocast(device_type=dev.type, dtype=torch.float16, enabled=dev.type == "cuda"):
                loss = lossf(net(torch.from_numpy(np.ascontiguousarray(x)).to(dev).float()), yt[i])
            opt.zero_grad(set_to_none=True)
            scaler.scale(loss).backward()
            scaler.step(opt)
            scaler.update()
    net.eval()
    models.append(net)
    log(f"trained seed {s}")
del BAGS

# ---- 2. full patch features + coords for the targets ----
want = set(targets)
feats, coords, slide_of = {p: [] for p in targets}, {p: [] for p in targets}, {p: [] for p in targets}
slide_names = {p: [] for p in targets}
files = sorted(glob.glob(os.path.join(a.parquet_dir, "*.parquet")))
for fi, f in enumerate(files):
    pf = pq.ParquetFile(f)
    for batch in pf.iter_batches(batch_size=16, columns=["case_submitter_id", "svs_filename", "coords", "embeddings"]):
        cases = batch.column("case_submitter_id").to_pylist()
        if not want.intersection(cases):
            continue
        emb, crd = batch.column("embeddings"), batch.column("coords")
        lens = emb.value_lengths().to_numpy()
        flat = emb.flatten().flatten().to_numpy(zero_copy_only=False)
        X = flat.reshape(-1, int(flat.size // lens.sum()))
        C = crd.flatten().flatten().to_numpy(zero_copy_only=False).reshape(-1, 2)
        start = 0
        for case, name, n in zip(cases, batch.column("svs_filename").to_pylist(), lens):
            if case in want and "-DX" in name:
                if name not in slide_names[case]:
                    slide_names[case].append(name)
                feats[case].append(X[start:start + n].astype(np.float16))
                coords[case].append(C[start:start + n])
                slide_of[case].append(np.full(n, slide_names[case].index(name)))
            start += n
    log(f"[{fi + 1}/{len(files)}] read {os.path.basename(f)}")

# ---- 3. attention per patient ----
rows = []
for p in targets:
    if not feats[p]:
        log(f"no DX patches found for {p}, skipped")
        continue
    X = torch.from_numpy(np.concatenate(feats[p]).astype(np.float32)).to(dev)
    C, S = np.concatenate(coords[p]), np.concatenate(slide_of[p])
    probs, atts = [], []
    with torch.no_grad():
        for net in models:
            logit, w = net(X, return_attention=True)
            probs.append(torch.sigmoid(logit).item())
            atts.append(w.float().cpu().numpy())
    w = np.mean(atts, axis=0)
    prob = float(np.mean(probs))
    n = len(w)
    ent = float(-(w * np.log(w + 1e-12)).sum() / np.log(n))
    srt = np.sort(w)[::-1]
    frac50 = float((np.searchsorted(np.cumsum(srt), 0.5) + 1) / n)
    lab = L.loc[p]
    if p in refs:
        group = "reference ER+" if lab.ER_IHC == 1 else "reference ER-"
    else:
        group = "IHC+/RNA-" if lab.ER_IHC == 1 else "IHC-/RNA+"
    pct = str(lab.get("ER_pct_IHC", ""))
    side = ("RNA" if (prob >= 0.5) == bool(lab.ESR1_RNApos) else "IHC") if p not in refs else "-"
    rows.append({"patient": p, "group": group, "ER_IHC": int(lab.ER_IHC), "ESR1_RNApos": int(lab.ESR1_RNApos),
                 "ESR1_log2": round(float(lab.ESR1_log2), 2), "ER_pct_IHC": pct, "ER_low": (pct == "<10%") and lab.ER_IHC == 1,
                 "PAM50": lab.PAM50, "P_ER_pos": round(prob, 3), "sides_with": side, "n_patches": n,
                 "n_slides": len(slide_names[p]), "attn_entropy_norm": round(ent, 3),
                 "frac_patches_for_50pct_attention": round(frac50, 4)})
    np.savez_compressed(os.path.join(a.out, f"{p}.npz"), coords=C, slide=S, attention=w, prob=prob,
                        slides=np.array(slide_names[p], dtype=str))
    # heatmap: one panel per slide, patches on their grid, colour = attention percentile within patient
    pr = w.argsort().argsort() / max(n - 1, 1)
    k = len(slide_names[p])
    fig, axes = plt.subplots(1, k, figsize=(5 * k, 5), squeeze=False)
    for j in range(k):
        m = S == j
        xy = C[m]
        step = np.diff(np.unique(xy[:, 0]))
        step = int(step[step > 0].min()) if (step > 0).any() else 1
        gx, gy = (xy[:, 0] - xy[:, 0].min()) // step, (xy[:, 1] - xy[:, 1].min()) // step
        img = np.full((gy.max() + 1, gx.max() + 1), np.nan)
        img[gy, gx] = pr[m]
        ax = axes[0, j]
        ax.imshow(img, cmap="inferno", vmin=0, vmax=1, interpolation="nearest")
        ax.set_title(slide_names[p][j].split(".")[0][-12:], fontsize=9)
        ax.axis("off")
    fig.suptitle(f"{p} | {group} | IHC {'+' if lab.ER_IHC == 1 else '-'} ({pct}) | ESR1 log2 {lab.ESR1_log2:.1f} | "
                 f"PAM50 {lab.PAM50} | P(ER+) {prob:.2f} -> sides with {side}", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(a.out, f"{p}.png"), dpi=110)
    plt.close(fig)

D = pd.DataFrame(rows)
D.to_csv(os.path.join(a.out, "summary.csv"), index=False)
log("\n" + D.drop(columns=["n_slides"]).to_string(index=False))
log("\nby group:\n" + D.groupby("group")[["P_ER_pos", "attn_entropy_norm", "frac_patches_for_50pct_attention"]]
    .mean().round(3).to_string())
dd = D[D.group.str.contains("RNA")]
log(f"discordant: sides with RNA {int((dd.sides_with == 'RNA').sum())} / {len(dd)}; "
    f"ER-low sides with RNA {int((dd[dd.ER_low].sides_with == 'RNA').sum())} / {int(dd.ER_low.sum())}")

# overview grid: ER-low discordant first, then other discordant, then references (first panel of each patient)
order = pd.concat([D[D.ER_low], D[~D.ER_low & D.group.str.contains("RNA")], D[D.group.str.startswith("reference")]])
show = order.head(24)
cols = 6
fig, axes = plt.subplots(int(np.ceil(len(show) / cols)), cols, figsize=(3.2 * cols, 3.4 * np.ceil(len(show) / cols)),
                         squeeze=False)
for ax in axes.ravel():
    ax.axis("off")
for ax, (_, r) in zip(axes.ravel(), show.iterrows()):
    z = np.load(os.path.join(a.out, f"{r.patient}.npz"))
    m = z["slide"] == 0
    xy, w = z["coords"][m], z["attention"][m]
    step = np.diff(np.unique(xy[:, 0]))
    step = int(step[step > 0].min()) if (step > 0).any() else 1
    gx, gy = (xy[:, 0] - xy[:, 0].min()) // step, (xy[:, 1] - xy[:, 1].min()) // step
    img = np.full((gy.max() + 1, gx.max() + 1), np.nan)
    img[gy, gx] = w.argsort().argsort() / max(len(w) - 1, 1)
    ax.imshow(img, cmap="inferno", vmin=0, vmax=1, interpolation="nearest")
    ax.set_title(f"{r.patient[5:]} {r.group}\nER {r.ER_pct_IHC} P={r.P_ER_pos:.2f} -> {r.sides_with}", fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(a.out, "overview.png"), dpi=110)
log("saved", a.out)
