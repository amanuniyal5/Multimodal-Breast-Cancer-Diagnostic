"""Baselines B2-B4 + preliminary 1b on TITAN slide embeddings (laptop CPU, ~5 min).

  - AUC per label under random vs site-grouped 5-fold CV (+ shuffled-label negative control)
  - site-only (no image) baseline
  - site identification from the embedding (10 largest sites)
  - 1b preliminary: train on ER-concordant cases, score the IHC != RNA cases

Usage: python scripts/baseline_titan.py [--emb data/features/titan_brca_20x.parquet]
                                        [--labels data/metadata/labels_patient.csv] [--out results]
"""
import argparse
import os
import warnings

import numpy as np
import pandas as pd
from scipy.stats import binomtest
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import OneHotEncoder

from common import bootstrap_auc_ci, linear_probe, load_titan

warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser()
ap.add_argument("--emb", default="data/features/titan_brca_20x.parquet")
ap.add_argument("--labels", default="data/metadata/labels_patient.csv")
ap.add_argument("--out", default="results")
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)

E = load_titan(a.emb)
L = pd.read_csv(a.labels, index_col=0)
L["HER2low_1plus"] = L.HER2_IHCscore.map({"0": 0, "1+": 1})
L["ER_IHC_shuffled"] = np.random.default_rng(0).permutation(L.ER_IHC.values)
tasks = ["ER_IHC", "ESR1_RNApos", "PR_IHC", "PGR_RNApos", "HER2_STRICT", "ERBB2_AMP",
         "HER2low_1plus", "ER_IHC_shuffled"]

res = []
for task in tasks:
    d = L[[task, "site"]].dropna().join(E, how="inner")
    y, X, g = d[task].astype(int).values, d[E.columns].values, d.site.values
    row = {"task": task, "n": len(y), "pos": int(y.sum())}
    for name, cv in [("random", StratifiedKFold(5, shuffle=True, random_state=0)),
                     ("site", StratifiedGroupKFold(5, shuffle=True, random_state=0))]:
        p = cross_val_predict(linear_probe(), X, y, cv=cv, groups=g if name == "site" else None,
                              method="predict_proba")[:, 1]
        lo, hi = bootstrap_auc_ci(y, p)
        row[f"AUC_{name}"] = round(roc_auc_score(y, p), 3)
        row[f"CI_{name}"] = f"{lo:.2f}-{hi:.2f}"
    res.append(row)
R = pd.DataFrame(res)
print(R.to_string(index=False))
R.to_csv(os.path.join(a.out, "baseline_results.csv"), index=False)

for task in ["ER_IHC", "ESR1_RNApos", "HER2_STRICT"]:
    d = L[[task, "site"]].dropna().join(E[[]], how="inner")
    y = d[task].astype(int).values
    p = cross_val_predict(make_pipeline(OneHotEncoder(handle_unknown="ignore"), LogisticRegression(max_iter=2000)),
                          d[["site"]], y, cv=StratifiedKFold(5, shuffle=True, random_state=0),
                          method="predict_proba")[:, 1]
    print(f"site-only (no image) AUC {task}: {roc_auc_score(y, p):.3f}")

top = L.site.value_counts().index[:10]
d = L[L.site.isin(top)][["site"]].join(E, how="inner")
pred = cross_val_predict(linear_probe(), d[E.columns].values, d.site.values,
                         cv=StratifiedKFold(5, shuffle=True, random_state=0))
print(f"site identification (10 largest sites, n={len(d)}): accuracy {(pred == d.site.values).mean():.3f} "
      f"vs chance {d.site.value_counts(normalize=True).max():.3f}")

d = L[["ER_IHC", "ESR1_RNApos"]].dropna().join(E, how="inner")
conc, disc = d[d.ER_IHC == d.ESR1_RNApos], d[d.ER_IHC != d.ESR1_RNApos]
m = linear_probe().fit(conc[E.columns].values, conc.ER_IHC.astype(int).values)
pred = (m.predict_proba(disc[E.columns].values)[:, 1] >= 0.5).astype(int)
rna = int((pred == disc.ESR1_RNApos.astype(int).values).sum())
print(f"1b prelim: discordant n={len(disc)}; sides with RNA {rna}, with IHC {len(disc) - rna}; "
      f"binomial p={binomtest(rna, len(disc)).pvalue:.3g}")
