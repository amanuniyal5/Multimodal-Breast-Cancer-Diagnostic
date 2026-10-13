"""Preliminary 1a (paired IHC-vs-RNA label comparison) + 1b robustness, on TITAN slide embeddings.

1a: same patients, same site-grouped folds; train one model on the IHC label and one on the RNA label;
    paired bootstrap of AUC(RNA) - AUC(IHC).
1b robustness: drop discordant cases whose ESR1 sits within 0.5 / 1.0 log2 of the GMM cut-off.

Usage: python scripts/paired_ihc_vs_rna.py [--emb ...] [--labels ...] [--esr1-cut 10.17]
"""
import argparse
import warnings

import numpy as np
import pandas as pd
from scipy.stats import binomtest
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import cross_val_predict

from common import linear_probe, load_titan, site_folds

warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser()
ap.add_argument("--emb", default="data/features/titan_brca_20x.parquet")
ap.add_argument("--labels", default="data/metadata/labels_patient.csv")
ap.add_argument("--esr1-cut", type=float, default=10.17, help="GMM cut-off printed by build_labels.py")
a = ap.parse_args()

E = load_titan(a.emb)
L = pd.read_csv(a.labels, index_col=0)


def paired(lab_a, lab_b, name):
    d = L[[lab_a, lab_b, "site"]].dropna().join(E, how="inner")
    X, g = d[E.columns].values, d.site.values
    cv = site_folds(d.index, 0)  # fixed site-grouped folds (data/splits)
    ya, yb = d[lab_a].astype(int).values, d[lab_b].astype(int).values
    pa = cross_val_predict(linear_probe(), X, ya, cv=cv, method="predict_proba")[:, 1]
    pb = cross_val_predict(linear_probe(), X, yb, cv=cv, method="predict_proba")[:, 1]
    r = np.random.default_rng(0)
    diffs = []
    for _ in range(2000):
        i = r.integers(0, len(d), len(d))
        if len(set(ya[i])) < 2 or len(set(yb[i])) < 2:
            continue
        diffs.append(roc_auc_score(yb[i], pb[i]) - roc_auc_score(ya[i], pa[i]))
    diffs = np.array(diffs)
    print(f"{name}: n={len(d)} | AUC {lab_a}={roc_auc_score(ya, pa):.3f} {lab_b}={roc_auc_score(yb, pb):.3f} "
          f"delta={diffs.mean():+.3f} 95%CI [{np.percentile(diffs, 2.5):+.3f}, {np.percentile(diffs, 97.5):+.3f}]")


paired("ER_IHC", "ESR1_RNApos", "1a ER")
paired("PR_IHC", "PGR_RNApos", "1a PR")
paired("HER2_STRICT", "ERBB2_AMP", "1a HER2")

d = L[["ER_IHC", "ESR1_RNApos", "ESR1_log2"]].dropna().join(E, how="inner")
conc, disc = d[d.ER_IHC == d.ESR1_RNApos], d[d.ER_IHC != d.ESR1_RNApos]
m = linear_probe().fit(conc[E.columns].values, conc.ER_IHC.astype(int))
disc = disc.assign(p=m.predict_proba(disc[E.columns].values)[:, 1])
for lab, sub in [("all", disc),
                 ("|ESR1-cut|>0.5", disc[(disc.ESR1_log2 - a.esr1_cut).abs() > 0.5]),
                 ("|ESR1-cut|>1.0", disc[(disc.ESR1_log2 - a.esr1_cut).abs() > 1.0])]:
    rna = int(((sub.p >= 0.5).astype(int) == sub.ESR1_RNApos.astype(int)).sum())
    print(f"1b {lab}: n={len(sub)} sides with RNA {rna}, IHC {len(sub) - rna}, "
          f"p={binomtest(rna, len(sub)).pvalue:.3g}")
