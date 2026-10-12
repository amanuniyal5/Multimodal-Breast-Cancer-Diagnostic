"""Option-1 experiments on patch features: ABMIL (GPU) and mean-pool logistic regression (CPU).

Same experiment definitions as scripts/extended_titan.py, but on UNI (or CONCH) patch features from build_bags.py:
  1a  paired label comparison per group (ER: IHC/RNA/RPPA, PR: IHC/RNA, HER2: strict/amp/RPPA), same patients+folds
  1b  train on concordant cases, score discordant cases, binomial test (prob averaged over seeds)
  1d  HER2 0 vs 1+ and 0 vs low, H&E model vs ERBB2 mRNA alone
  ctl shuffled ER label (negative control)
All CV is patient-level and grouped by tissue source site.

Usage: python scripts/train_mil.py --feats data/features --name uni_v1_20x --cap 4096 \
           --labels data/metadata/labels_patient.csv --out results --models mean,abmil --seeds 3
"""
import argparse
import os
import time
import warnings

import numpy as np
import pandas as pd
from scipy.stats import binomtest
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.mixture import GaussianMixture
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser()
ap.add_argument("--feats", default="data/features")
ap.add_argument("--name", default="uni_v1_20x")
ap.add_argument("--cap", type=int, default=4096)
ap.add_argument("--labels", default="data/metadata/labels_patient.csv")
ap.add_argument("--clinical", default="data/metadata/clinical_patient_brca.txt")
ap.add_argument("--out", default="results")
ap.add_argument("--models", default="mean,abmil")
ap.add_argument("--seeds", type=int, default=3)
ap.add_argument("--epochs", type=int, default=15)
ap.add_argument("--train-patches", type=int, default=2048)
ap.add_argument("--lr", type=float, default=2e-4)
ap.add_argument("--pr-cut", default="gmm", choices=["gmm", "prevmatch"],
                help="PGR RNA label: GMM cut-off or threshold matching PR IHC prevalence")
ap.add_argument("--only", default="all", choices=["all", "pr"], help="pr: run only the PR 1a group and PR 1b")
ap.add_argument("--tag", default="", help="suffix for output file names")
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
SEEDS = list(range(a.seeds))
MODELS = a.models.split(",")
t_start = time.time()


def log(*x):
    print(f"[{time.time() - t_start:6.0f}s]", *x, flush=True)


# ---------------- labels (same definitions as extended_titan.py) ----------------
L = pd.read_csv(a.labels, index_col=0)
c = pd.read_csv(a.clinical, sep="\t", skiprows=[1, 2], dtype=str).set_index("bcr_patient_barcode")


def gmm_pos(x):
    x = x.dropna()
    g = GaussianMixture(2, random_state=0).fit(x.values.reshape(-1, 1))
    return pd.Series(g.predict(x.values.reshape(-1, 1)) == np.argmax(g.means_.ravel()), index=x.index).astype(float)


if a.pr_cut == "prevmatch":
    q = L.PGR_log2.quantile(1 - L.PR_IHC.mean())
    L["PGR_RNApos"] = (L.PGR_log2 >= q).astype(float).where(L.PGR_log2.notna())
L["ER_RPPApos"] = gmm_pos(L.ER_RPPA).reindex(L.index)
L["HER2_RPPApos"] = gmm_pos(L.HER2_RPPA).reindex(L.index)
L["HER2low_1plus"] = L.HER2_IHCscore.map({"0": 0, "1+": 1})
fish = c.her2_fish_status.reindex(L.index)
low = (L.HER2_IHCscore == "1+") | ((L.HER2_IHCscore == "2+") & (fish == "Negative"))
L["HER2_0_vs_low"] = np.where(L.HER2_IHCscore == "0", 0.0, np.where(low, 1.0, np.nan))
L["ER_IHC_shuffled"] = np.random.default_rng(0).permutation(L.ER_IHC.values)

# ---------------- features ----------------
idx = pd.read_csv(os.path.join(a.feats, f"{a.name}_index.csv")).set_index("patient")
mean_npz = np.load(os.path.join(a.feats, f"{a.name}_mean.npz"), allow_pickle=True)
MEAN = pd.DataFrame(mean_npz["X"], index=mean_npz["patients"])
L = L[L.index.isin(idx.index)]
log(f"patients with features and labels: {len(L)}")

BAGS = None
if "abmil" in MODELS:
    import torch
    import torch.nn as nn

    dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    log("device", dev, torch.cuda.get_device_name(0) if dev.type == "cuda" else "")
    t = time.time()
    BAGS = np.load(os.path.join(a.feats, f"{a.name}_cap{a.cap}.npy"))  # fully into RAM (~9 GB for 1,060 x 4096 x 1024)
    log(f"bags loaded {BAGS.shape} in {time.time() - t:.0f}s")

    class ABMIL(nn.Module):  # gated attention MIL (Ilse et al. 2018), CLAM-style sizes
        def __init__(self, d, h=256, att=128, p=0.25):
            super().__init__()
            self.fc = nn.Sequential(nn.Linear(d, h), nn.ReLU(), nn.Dropout(p))
            self.a = nn.Sequential(nn.Linear(h, att), nn.Tanh())
            self.b = nn.Sequential(nn.Linear(h, att), nn.Sigmoid())
            self.c = nn.Linear(att, 1)
            self.cls = nn.Linear(h, 1)

        def forward(self, x):
            hdn = self.fc(x)
            w = torch.softmax(self.c(self.a(hdn) * self.b(hdn)).squeeze(-1), dim=0)
            return self.cls((w[:, None] * hdn).sum(0)).squeeze()


def bag(p, train, rng):
    r, k = idx.at[p, "row"], idx.at[p, "n_kept"]
    x = BAGS[r, :k]
    if train and k > a.train_patches:
        x = x[rng.choice(k, a.train_patches, replace=False)]
    return torch.from_numpy(np.ascontiguousarray(x)).to(dev, non_blocking=True)


def fit_predict(model, tr, te, y_tr, seed):
    """tr/te: lists of patient ids; returns probabilities for te"""
    if model == "mean":
        clf = make_pipeline(StandardScaler(), LogisticRegression(C=0.01, max_iter=3000, class_weight="balanced"))
        clf.fit(MEAN.loc[tr].values, y_tr)
        return clf.predict_proba(MEAN.loc[te].values)[:, 1]
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    net = ABMIL(BAGS.shape[2]).to(dev)
    opt = torch.optim.Adam(net.parameters(), lr=a.lr, weight_decay=1e-4)
    pos = max(y_tr.sum(), 1)
    lossf = nn.BCEWithLogitsLoss(pos_weight=torch.tensor((len(y_tr) - pos) / pos, device=dev))
    scaler = torch.cuda.amp.GradScaler(enabled=dev.type == "cuda")
    yt = torch.tensor(y_tr, dtype=torch.float32, device=dev)
    for _ in range(a.epochs):
        net.train()
        for i in rng.permutation(len(tr)):
            with torch.autocast(device_type=dev.type, dtype=torch.float16, enabled=dev.type == "cuda"):
                loss = lossf(net(bag(tr[i], True, rng).float()), yt[i])
            opt.zero_grad(set_to_none=True)
            scaler.scale(loss).backward()
            scaler.step(opt)
            scaler.update()
    net.eval()
    out = []
    with torch.no_grad(), torch.autocast(device_type=dev.type, dtype=torch.float16, enabled=dev.type == "cuda"):
        for p in te:
            out.append(torch.sigmoid(net(bag(p, False, rng).float())).item())
    return np.array(out)


def oof(cols, model, seed, require=()):
    d = L[cols + list(require) + ["site"]].dropna()
    pts, g = d.index.to_numpy(), d.site.values
    folds = list(StratifiedGroupKFold(5, shuffle=True, random_state=seed).split(pts, d[cols[0]].astype(int), g))
    out = {col: np.zeros(len(d)) for col in cols}
    for tr, te in folds:
        for col in cols:
            out[col][te] = fit_predict(model, list(pts[tr]), list(pts[te]), d[col].astype(int).values[tr], seed)
    return d, out


rows, preds = [], []
for model in MODELS:
    log(f"===== model: {model} =====")
    # 1a paired groups (also gives the per-label AUC table)
    groups = [["ER_IHC", "ESR1_RNApos", "ER_RPPApos"], ["PR_IHC", "PGR_RNApos"],
              ["HER2_STRICT", "ERBB2_AMP", "HER2_RPPApos"]]
    for group in (groups if a.only == "all" else [["PR_IHC", "PGR_RNApos"]]):
        per = {col: [] for col in group}
        for s in SEEDS:
            d, o = oof(group, model, s)
            for col in group:
                per[col].append(roc_auc_score(d[col].astype(int), o[col]))
                preds.append(pd.DataFrame({"model": model, "exp": "1a", "label": col, "seed": s, "patient": d.index,
                                           "y": d[col].values, "p": o[col]}))
            log(f"1a {group} seed {s}: " + ", ".join(f"{k}={v[-1]:.3f}" for k, v in per.items()))
        for col in group:
            rows.append({"model": model, "exp": "1a", "label": col, "n": len(d), "pos": int(d[col].sum()),
                         "AUC_mean": np.mean(per[col]), "AUC_sd": np.std(per[col]),
                         "delta_vs_IHC_mean": np.mean(np.subtract(per[col], per[group[0]])),
                         "delta_vs_IHC_sd": np.std(np.subtract(per[col], per[group[0]]))})
    # negative control and 1d
    for col in (["ER_IHC_shuffled", "HER2low_1plus", "HER2_0_vs_low"] if a.only == "all" else []):
        aucs = []
        for s in SEEDS:
            # 1d: same patients for the H&E model and the ERBB2-mRNA comparator
            d, o = oof([col], model, s, require=["ERBB2_log2"] if col.startswith("HER2") else [])
            aucs.append(roc_auc_score(d[col].astype(int), o[col]))
        extra = {"AUC_ERBB2_mRNA": roc_auc_score(d[col].astype(int), L.loc[d.index, "ERBB2_log2"])} \
            if col.startswith("HER2") else {}
        rows.append({"model": model, "exp": "ctl" if "shuffled" in col else "1d", "label": col, "n": len(d),
                     "pos": int(d[col].sum()), "AUC_mean": np.mean(aucs), "AUC_sd": np.std(aucs), **extra})
        log(f"{col}: AUC {np.mean(aucs):.3f} ± {np.std(aucs):.3f} {extra}")
    # 1b which side
    pairs = [("ER_IHC", "ESR1_RNApos"), ("ER_IHC", "ER_RPPApos"), ("PR_IHC", "PGR_RNApos"),
             ("HER2_STRICT", "ERBB2_AMP")]
    for ihc, mol in (pairs if a.only == "all" else [("PR_IHC", "PGR_RNApos")]):
        d = L[[ihc, mol]].dropna()
        conc, disc = d[d[ihc] == d[mol]], d[d[ihc] != d[mol]]
        prob = np.mean([fit_predict(model, list(conc.index), list(disc.index), conc[ihc].astype(int).values, s)
                        for s in SEEDS], axis=0)
        k = int(((prob >= 0.5).astype(int) == disc[mol].astype(int).values).sum())
        pval = binomtest(k, len(disc)).pvalue
        rows.append({"model": model, "exp": "1b", "label": f"{ihc} vs {mol}", "n": len(disc), "side_molecular": k,
                     "side_IHC": len(disc) - k, "p": pval})
        preds.append(pd.DataFrame({"model": model, "exp": "1b", "label": f"{ihc} vs {mol}", "seed": -1,
                                   "patient": disc.index, "y": disc[mol].values, "p": prob}))
        log(f"1b {ihc} vs {mol}: n={len(disc)} molecular {k} / IHC {len(disc) - k}, p={pval:.3g}")
    pd.DataFrame(rows).to_csv(os.path.join(a.out, f"mil_{a.name}{a.tag}_summary.csv"), index=False)
    pd.concat(preds).to_csv(os.path.join(a.out, f"mil_{a.name}{a.tag}_oof.csv.gz"), index=False)

log("saved", os.path.join(a.out, f"mil_{a.name}{a.tag}_summary.csv"))
