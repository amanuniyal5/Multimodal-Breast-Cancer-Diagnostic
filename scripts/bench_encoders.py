"""X1: encoder benchmark with the mean-pool linear probe (docs/DESIGN_SPEC.md §5.1).

Same experiment definitions as train_mil.py (1a label groups on the same patients, 1b on discordant patients, 1d
HER2-low, shuffled control), run for every encoder on ONE common cohort (patients with features for all encoders), so
encoders are compared on identical patients and folds:
  aucs.csv        per encoder x label: AUC (seed-averaged predictions) with DeLong CI, SD over seeds, AUPRC,
                  difference vs the reference encoder (DeLong, same label and patients)
  preference.csv  per encoder x label pair: ΔAUC molecular − IHC (generalised DeLong) and 1b counts (binomial)
  site.csv        hospital identification accuracy (10 largest sites, as baseline_titan.py)
  curves.csv      label efficiency: AUC vs number of training patients (seed-0 folds, 5 stratified repeats)
  oof.csv.gz      seed-averaged out-of-fold probabilities (for fusion / calibration / reliability)
--cohort subset restricts to the pre-set random ResNet-50 subset (data/metadata/dx_manifest.csv, subset == True).
Encoders marked tcga_pretrained in configs/encoders.json are leakage contrasts (†), never "best".

Usage: python scripts/bench_encoders.py --out results/bench/x1 [--encoders uni_v1,conch_v1,...] [--cohort subset]
"""
import argparse
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.stats import binomtest
from sklearn.model_selection import StratifiedKFold, cross_val_predict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bench import folds, registry, stats, tasks  # noqa: E402
from bench.metrics import classification  # noqa: E402
from common import linear_probe  # noqa: E402

warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser()
ap.add_argument("--encoders", default="", help="comma list; default: every encoder with features")
ap.add_argument("--ref", default="uni_v1")
ap.add_argument("--cohort", default="all", choices=["all", "subset"])
ap.add_argument("--seeds", type=int, default=5)
ap.add_argument("--sizes", default="25,50,100,200,400,0", help="training-set sizes for the curves (0 = all)")
ap.add_argument("--reps", type=int, default=5)
ap.add_argument("--curve-tasks", default="ER_IHC,ESR1_RNApos,PR_IHC,HER2_STRICT")
ap.add_argument("--no-curves", action="store_true")
ap.add_argument("--std", action="store_true", help="also the mean+SD variant for encoders that store SD")
ap.add_argument("--n-jobs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)))
ap.add_argument("--out", default="results/bench/x1")
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
t0 = time.time()


def log(*x):
    print(f"[{time.time() - t0:6.0f}s]", *x, flush=True)


L = tasks.labels()
ENCS = a.encoders.split(",") if a.encoders else registry.available()
FE = {e: registry.load_means(e) for e in ENCS}
if a.std:
    for e in list(ENCS):
        s = registry.load_means(e, std=True)
        if s is not None:
            FE[e + "+sd"] = s
ENCS = list(FE)
common = set(L.index)
for e in ENCS:
    common &= set(FE[e].index)
if a.cohort == "subset":
    m = pd.read_csv("data/metadata/dx_manifest.csv")
    common &= set(m[m.subset].patient)
common = sorted(common & set(pd.read_csv("data/splits/folds.csv").patient))
log(f"encoders {ENCS}; common cohort {len(common)} patients")
REF = a.ref if a.ref in ENCS else ENCS[0]

GROUPS = [["ER_IHC", "ESR1_RNApos", "ER_RPPApos"], ["PR_IHC", "PGR_prevmatch"], ["HER2_STRICT", "ERBB2_AMP"],
          ["HER2low_1plus"], ["HER2_0_vs_low"], ["ER_IHC_shuffled"]]
L["ER_IHC_shuffled"] = np.random.default_rng(0).permutation(L.ER_IHC.values)


def oof(enc, pts, y, seed, train_idx_fn=None):
    """seed's out-of-fold probabilities for patients pts (positional)"""
    X = FE[enc].loc[pts].values
    p = np.zeros(len(pts))
    for tr, te in folds.outer(pts, seed):
        if train_idx_fn is not None:
            tr = tr[train_idx_fn(y[tr])]
        p[te] = linear_probe().fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]
    return p


def run_group(enc, group):
    req = group + (["ERBB2_log2"] if group[0].startswith("HER2") and len(group) == 1 else [])
    d = L.loc[common, req].dropna()
    pts = d.index.to_numpy()
    out = []
    for col in group:
        y = d[col].astype(int).to_numpy()
        per = [oof(enc, pts, y, s) for s in range(a.seeds)]
        out.append((enc, tuple(group), col, pts, y, np.mean(per, 0), [stats.safe_auc(y, q) for q in per]))
    return out


def run_1b(enc, ihc, mol):
    d = L.loc[common, [ihc, mol]].dropna()
    conc, disc = d[d[ihc] == d[mol]], d[d[ihc] != d[mol]]
    clf = linear_probe().fit(FE[enc].loc[conc.index].values, conc[ihc].astype(int).values)
    prob = clf.predict_proba(FE[enc].loc[disc.index].values)[:, 1]
    k = int(((prob >= 0.5).astype(int) == disc[mol].astype(int).values).sum())
    return {"encoder": enc, "pair": f"{ihc} vs {mol}", "n_disc": len(disc), "side_molecular": k,
            "side_IHC": len(disc) - k, "p_binom": binomtest(k, len(disc)).pvalue}


def run_curve(enc, task, size, rep):
    d = L.loc[common, [task]].dropna()
    pts, y = d.index.to_numpy(), d[task].astype(int).to_numpy()
    p = oof(enc, pts, y, 0, train_idx_fn=lambda yt: folds.subsample(yt, size if size else len(yt), rep, 0))
    return {"encoder": enc, "task": task, "n_train": size if size else "all", "rep": rep, "auc": stats.safe_auc(y, p)}


# ---- 1a groups, 1d, control ----
res = Parallel(n_jobs=a.n_jobs)(delayed(run_group)(e, g) for e in ENCS for g in GROUPS)
res = [r for g in res for r in g]
by = {(r[0], r[2], r[1]): r for r in res}
rows, oofs = [], []
for enc, group, col, pts, y, p, seed_aucs in res:
    auc, lo, hi = stats.auc_ci(y, p)
    row = {"encoder": enc, "label": registry.label(enc.replace("+sd", "")) + (" +SD" if enc.endswith("+sd") else ""),
           "task": col, "group": "|".join(group), "n": len(y), "n_pos": int(y.sum()), "auc": auc, "lo": lo, "hi": hi,
           "auc_sd_seeds": np.std(seed_aucs), **{k: v for k, v in classification(y, p).items() if k != "auc"}}
    if enc != REF:
        r = by[(REF, col, group)]
        dl = stats.delong(y, r[5], y, p)
        row.update(diff_vs_ref=dl["diff"], diff_lo=dl["lo"], diff_hi=dl["hi"], diff_p=dl["p"])
    if col.startswith("HER2") and len(group) == 1:
        row["auc_ERBB2_mRNA"] = stats.safe_auc(y, L.loc[pts, "ERBB2_log2"].values)
    rows.append(row)
    oofs.append(pd.DataFrame({"encoder": enc, "task": col, "group": "|".join(group), "patient": pts, "y": y, "p": p}))
A = pd.DataFrame(rows)
A.to_csv(os.path.join(a.out, "aucs.csv"), index=False)
pd.concat(oofs).to_csv(os.path.join(a.out, "oof.csv.gz"), index=False)
log("AUC table:\n" + A.pivot(index="encoder", columns="task", values="auc").round(3).to_string())

# ---- label preference: 1a ΔAUC and 1b ----
pref = []
for e in ENCS:
    for g in GROUPS[:3]:
        for mol in g[1:]:
            ri, rm = by[(e, g[0], tuple(g))], by[(e, mol, tuple(g))]
            dl = stats.delong(ri[4], ri[5], rm[4], rm[5])
            pref.append({"encoder": e, "pair": f"{g[0]} vs {mol}", "n": len(ri[4]), "auc_ihc": dl["auc_ref"],
                         "auc_mol": dl["auc_alt"], "delta": dl["diff"], "lo": dl["lo"], "hi": dl["hi"], "p": dl["p"]})
B = pd.DataFrame(Parallel(n_jobs=a.n_jobs)(delayed(run_1b)(e, i, m) for e in ENCS for i, m in tasks.PAIRS))
P = pd.DataFrame(pref).merge(B, on=["encoder", "pair"], how="outer")
P.to_csv(os.path.join(a.out, "preference.csv"), index=False)
log("label preference:\n" + P[["encoder", "pair", "delta", "lo", "hi", "side_molecular", "side_IHC", "p_binom"]]
    .round(3).to_string(index=False))

# ---- hospital identification ----
top = L.loc[common, "site"].value_counts().index[:10]
d = L.loc[common][L.loc[common, "site"].isin(top)]
site_rows = []
for e in ENCS:
    pred = cross_val_predict(linear_probe(), FE[e].loc[d.index].values, d.site.to_numpy(dtype=str),
                             cv=StratifiedKFold(5, shuffle=True, random_state=0), n_jobs=a.n_jobs)
    acc = (pred == d.site.to_numpy(dtype=str)).mean()
    lo, hi = stats.wilson(int(round(acc * len(d))), len(d))
    site_rows.append({"encoder": e, "n": len(d), "accuracy": acc, "lo": lo, "hi": hi,
                      "chance": d.site.value_counts(normalize=True).max()})
S = pd.DataFrame(site_rows)
S.to_csv(os.path.join(a.out, "site.csv"), index=False)
log("hospital identification:\n" + S.round(3).to_string(index=False))

# ---- label-efficiency curves ----
if not a.no_curves:
    sizes = [int(s) for s in a.sizes.split(",")]
    jobs = [(e, t, s, r) for e in ENCS if not e.endswith("+sd") for t in a.curve_tasks.split(",")
            for s in sizes for r in (range(a.reps) if s else [0])]
    C = pd.DataFrame(Parallel(n_jobs=a.n_jobs)(delayed(run_curve)(*j) for j in jobs))
    C.to_csv(os.path.join(a.out, "curves.csv"), index=False)
    log("curves (mean AUC over repeats):\n" +
        C.groupby(["task", "encoder", "n_train"], sort=False).auc.mean().unstack().round(3).to_string())
log("BENCH_ENCODERS_OK")
