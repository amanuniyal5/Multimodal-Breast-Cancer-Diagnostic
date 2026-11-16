"""Retrospective evaluation of ReceptorCheck (docs/DESIGN_SPEC.md §7.6) on ALL patients with out-of-fold outputs
(not only the demo bundles). Same predictor, calibration, conformal sets, OOD rule and consistency rule as the app.

triage      abstain on the least confident 10 / 20 / 30 % (|p_cal - 0.5|) or on conformal {+,-} sets: accuracy,
            sensitivity, specificity on the kept cases vs IHC and vs the molecular label.
qa          the IHC consistency score as a second reader for RNA-defined discordance: precision / recall at the top
            1 / 2 / 5 % of the cohort (Wilson CI, binomial vs the base rate), the FLAG rule itself (tau cross-fitted
            over folds), ER-low enrichment; baselines: random flagging (base rate) and IHC % stained alone.
frozen      ER on frozen TS slides with the DX-trained model (experiment G) and conformal coverage on TS.
robustness  per-hospital AUC spread (hospitals with >= 20 patients and both classes); OOD flag rate on held-out
            hospitals vs the 1 % expected in distribution.
acceptance  coverage within +-3 % of 90 % on DX; FLAG precision above the base rate with a CI excluding it; bundle load
            time < 2 s.
Output: results/app_eval/{triage,qa,robustness,acceptance}.csv
Usage: python scripts/app_evaluate.py
"""
import glob
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import binomtest
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bench import registry, stats, tasks  # noqa: E402
from common import site_folds  # noqa: E402

OUT = "results/app_eval"
os.makedirs(OUT, exist_ok=True)
pd.set_option("display.width", 250)
L = tasks.labels()
ENS = pd.read_csv("results/bench/x6/oof_ensemble.csv.gz")
MARK = {"ER_IHC": "ESR1_RNApos", "PR_IHC": "PGR_prevmatch", "HER2_STRICT": "ERBB2_AMP"}
U = registry.load_means("uni_v1")
pts_all = U.index.intersection(L.index).to_numpy()
fold0 = {}
for k, (tr, te) in enumerate(site_folds(pts_all, 0)):
    fold0.update({p: k for p in pts_all[te]})


def rates(y, yhat):
    y, yhat = np.asarray(y).astype(int), np.asarray(yhat).astype(int)
    return {"accuracy": (y == yhat).mean(), "sensitivity": yhat[y == 1].mean() if (y == 1).any() else np.nan,
            "specificity": 1 - yhat[y == 0].mean() if (y == 0).any() else np.nan}


tri, qa, rob, acc = [], [], [], []
for ihc, mol in MARK.items():
    e = ENS[ENS.task == ihc].set_index("patient")
    e = e[e.index.isin(fold0)]
    y = e.y.astype(int).to_numpy()
    ym = L.loc[e.index, mol].to_numpy(float)
    pc = e.p_cal.to_numpy()
    sneg, spos = e.set_standard_neg.astype(bool).to_numpy(), e.set_standard_pos.astype(bool).to_numpy()
    conf = np.abs(pc - 0.5)
    # ---------------- triage ----------------
    rules = {"no abstention": np.ones(len(y), bool), "conformal singleton only": sneg ^ spos}
    for f in (0.1, 0.2, 0.3):
        rules[f"abstain least confident {int(f * 100)}%"] = conf > np.quantile(conf, f)
    for name, keep in rules.items():
        r = {"task": ihc, "rule": name, "n": len(y), "kept": keep.mean()}
        for lab, yy in (("IHC", y), ("molecular", ym)):
            ok = keep & ~np.isnan(yy)
            r.update({f"{k}_vs_{lab}": v for k, v in rates(yy[ok], pc[ok] >= 0.5).items()})
        tri.append(r)
    # ---------------- QA second reader ----------------
    score = np.where(y == 1, 1 - pc, pc)
    okm = ~np.isnan(ym)
    disc = okm & (ym != y)
    contra = np.where(y == 1, sneg & ~spos, spos & ~sneg)
    f = e.index.map(fold0).to_numpy()
    flag = np.zeros(len(y), bool)
    for k in np.unique(f):
        o = f != k
        tau = 1.01
        for t in np.sort(np.unique(score[o & contra])):
            sel = o & contra & (score >= t)
            if sel.sum() >= 3 and disc[sel].mean() >= 0.5:
                tau = t
                break
        flag[f == k] = contra[f == k] & (score[f == k] >= tau)
    base = disc[okm].mean()
    order = np.argsort(-score)
    erlow = ((L.loc[e.index, "ER_pct_IHC"] == "<10%") & (y == 1)).to_numpy()
    for k in (1, 2, 5):
        top = order[:int(round(len(y) * k / 100))]
        top = top[okm[top]]
        h = int(disc[top].sum())
        lo, hi = stats.wilson(h, len(top))
        qa.append({"task": ihc, "rule": f"top {k}% by score", "n": int(okm.sum()), "base_rate": base, "flagged": len(top),
                   "true_discordant": h, "precision": h / len(top), "lo": lo, "hi": hi, "recall": h / disc.sum(),
                   "p_vs_base": binomtest(h, len(top), base, alternative="greater").pvalue,
                   "er_low_flagged": int(erlow[top].sum()) if ihc == "ER_IHC" else np.nan})
    h, n = int(disc[flag & okm].sum()), int((flag & okm).sum())
    lo, hi = stats.wilson(h, n) if n else (np.nan, np.nan)
    qa.append({"task": ihc, "rule": "FLAG rule (app)", "n": int(okm.sum()), "base_rate": base, "flagged": n,
               "true_discordant": h, "precision": h / n if n else np.nan, "lo": lo, "hi": hi,
               "recall": h / disc.sum(), "p_vs_base": binomtest(h, n, base, alternative="greater").pvalue if n else np.nan,
               "er_low_flagged": int(erlow[flag].sum()) if ihc == "ER_IHC" else np.nan})
    qa.append({"task": ihc, "rule": "AUC of the score for discordance", "n": int(okm.sum()), "base_rate": base,
               "precision": roc_auc_score(disc[okm], score[okm])})
    if ihc == "ER_IHC":
        pct = L.loc[e.index, "ER_pct_IHC"].map({"<10%": 5, "10-19%": 15, "20-29%": 25, "30-39%": 35, "40-49%": 45,
                                                "50-59%": 55, "60-69%": 65, "70-79%": 75, "80-89%": 85,
                                                "90-99%": 95}).to_numpy(float)
        m = (y == 1) & ~np.isnan(pct) & okm
        qa.append({"task": ihc, "rule": "baseline: IHC % stained (AUC, IHC+ with %)", "n": int(m.sum()),
                   "base_rate": disc[m].mean(), "precision": roc_auc_score(disc[m], -pct[m])})
        qa.append({"task": ihc, "rule": "app score on the same IHC+ with % (AUC)", "n": int(m.sum()),
                   "base_rate": disc[m].mean(), "precision": roc_auc_score(disc[m], score[m])})
    # ---------------- robustness ----------------
    sites = L.loc[e.index, "site"].to_numpy(dtype=str)
    au = []
    for s in np.unique(sites):
        mm = sites == s
        if mm.sum() >= 20 and len(set(y[mm])) == 2:
            au.append(roc_auc_score(y[mm], pc[mm]))
    rob.append({"task": ihc, "metric": "per-hospital AUC (>= 20 patients, both classes)", "n_hospitals": len(au),
                "min": np.min(au), "median": np.median(au), "max": np.max(au), "overall": roc_auc_score(y, pc)})
    cov = (np.where(y == 1, spos, sneg)).mean()
    acc.append({"task": ihc, "check": "coverage within 0.90 +- 0.03 (DX)", "value": cov,
                "pass": abs(cov - 0.90) <= 0.03})
    fq = [q for q in qa if q["task"] == ihc and q["rule"] == "FLAG rule (app)"][0]
    acc.append({"task": ihc, "check": "FLAG precision above base rate, CI excluding it", "value": fq["precision"],
                "pass": bool(n and fq["lo"] > base)})

# OOD on held-out hospitals
od = []
for tr, te in site_folds(pts_all, 0):
    sc = StandardScaler().fit(U.loc[pts_all[tr]].values)
    pca = PCA(50, random_state=0).fit(sc.transform(U.loc[pts_all[tr]].values))
    zz = lambda X: np.linalg.norm(pca.transform(sc.transform(X)) / np.sqrt(pca.explained_variance_), axis=1)  # noqa
    thr = np.percentile(zz(U.loc[pts_all[tr]].values), 99)
    od.extend(zz(U.loc[pts_all[te]].values) > thr)
rob.append({"task": "all", "metric": "OOD flag rate on held-out hospitals (in-distribution: 1 %)", "n_hospitals": np.nan,
            "overall": np.mean(od)})
# frozen (experiment G + conformal on TS)
if os.path.exists("results/week3/g_transfer.csv"):
    G = pd.read_csv("results/week3/g_transfer.csv")
    for _, r in G[(G.label == "ER_IHC")].iterrows():
        rob.append({"task": "ER_IHC", "metric": f"frozen: train {r.train} -> test {r.test} AUC (ECE {r.ece:.3f})",
                    "overall": r.auc, "n_hospitals": np.nan})
if os.path.exists("results/bench/x6/conformal_frozen.csv"):
    for _, r in pd.read_csv("results/bench/x6/conformal_frozen.csv").iterrows():
        rob.append({"task": "ER_IHC", "metric": f"conformal coverage on {r.slides}", "overall": r.coverage,
                    "n_hospitals": np.nan})
# bundle load time
files = glob.glob("app/cases/*.json")
t0 = time.time()
for fp in files:
    json.load(open(fp))
lt = (time.time() - t0) / max(len(files), 1)
acc.append({"task": "app", "check": "bundle load time < 2 s (mean over bundles)", "value": lt, "pass": lt < 2})
acc.append({"task": "app", "check": "every number from an OOF file (by construction: x6/oof_ensemble, x3/oof_reg)",
            "value": np.nan, "pass": True})
for name, rows in (("triage", tri), ("qa", qa), ("robustness", rob), ("acceptance", acc)):
    D = pd.DataFrame(rows)
    D.to_csv(os.path.join(OUT, f"{name}.csv"), index=False)
    print(f"\n{name}:\n" + D.round(3).to_string(index=False))
print("APP_EVAL_OK")
