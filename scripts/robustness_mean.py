"""Robustness baselines across encoders (TITAN slide embedding, UNI mean-pool, CONCH mean-pool). Laptop CPU.

R1  PR cut-off sensitivity: PGR GMM / prevalence-matched / confident-only labels -> 1a (paired AUC vs PR IHC)
    and 1b (which side on discordant cases). Tests whether the PR "sides with IHC" reversal is a cut-off artefact.
R2  ER 1b broken down by discordance direction (IHC+/RNA- vs IHC-/RNA+) and by IHC % stained (ER-low < 10%,
    10-39%, >= 40%, unknown): which discordant cases drive the "sides with RNA" result.

Usage: python scripts/robustness_mean.py [--seeds 3] [--out results]
"""
import argparse
import os
import warnings

import numpy as np
import pandas as pd
from scipy.stats import binomtest
from sklearn.metrics import roc_auc_score
from sklearn.mixture import GaussianMixture
from sklearn.model_selection import cross_val_predict

from common import linear_probe, load_titan, site_folds

warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser()
ap.add_argument("--labels", default="data/metadata/labels_patient.csv")
ap.add_argument("--feats", default="data/features")
ap.add_argument("--seeds", type=int, default=3)
ap.add_argument("--out", default="results")
a = ap.parse_args()


def load_mean(name):
    z = np.load(os.path.join(a.feats, f"{name}_mean.npz"), allow_pickle=True)
    return pd.DataFrame(z["X"], index=[str(p) for p in z["patients"]])


ENC = {"TITAN": load_titan(os.path.join(a.feats, "titan_brca_20x.parquet")),
       "UNI_mean": load_mean("uni_v1_20x"), "CONCH_mean": load_mean("conch_v1_20x")}
L = pd.read_csv(a.labels, index_col=0)

# PR labels under three cut-off schemes
x = L.PGR_log2.dropna()
g = GaussianMixture(2, random_state=0).fit(x.values.reshape(-1, 1))
post = pd.Series(g.predict_proba(x.values.reshape(-1, 1))[:, np.argmax(g.means_.ravel())], index=x.index)
q = x.quantile(1 - L.PR_IHC.mean())
L["PGR_prevmatch"] = (L.PGR_log2 >= q).astype(float).where(L.PGR_log2.notna())
L["PGR_confident"] = pd.Series(np.where(post > 0.9, 1.0, np.where(post < 0.1, 0.0, np.nan)), index=x.index).reindex(L.index)
print(f"PGR cut-offs: GMM pos {int(L.PGR_RNApos.sum())}, prevalence-matched (log2 >= {q:.2f}) pos "
      f"{int(L.PGR_prevmatch.sum())}, confident-only n={int(L.PGR_confident.notna().sum())}")

rows = []
for enc, E in ENC.items():
    # R1a: paired PR comparison per cut-off
    for mol in ["PGR_RNApos", "PGR_prevmatch", "PGR_confident"]:
        d = L[["PR_IHC", mol, "site"]].dropna().join(E, how="inner")
        X = d[E.columns].values
        deltas, base = [], []
        for s in range(a.seeds):
            folds = site_folds(d.index, s)  # fixed site-grouped folds (data/splits)
            pi = cross_val_predict(linear_probe(), X, d.PR_IHC.astype(int), cv=folds, method="predict_proba")[:, 1]
            pm = cross_val_predict(linear_probe(), X, d[mol].astype(int), cv=folds, method="predict_proba")[:, 1]
            ai, am = roc_auc_score(d.PR_IHC.astype(int), pi), roc_auc_score(d[mol].astype(int), pm)
            base.append(ai)
            deltas.append(am - ai)
        # R1b: which side
        dd = L[["PR_IHC", mol]].dropna().join(E, how="inner")
        conc, disc = dd[dd.PR_IHC == dd[mol]], dd[dd.PR_IHC != dd[mol]]
        m = linear_probe().fit(conc[E.columns].values, conc.PR_IHC.astype(int))
        k = int(((m.predict_proba(disc[E.columns].values)[:, 1] >= 0.5).astype(int) == disc[mol].astype(int).values).sum())
        p = binomtest(k, len(disc)).pvalue
        rows.append({"exp": "R1 PR cut-off", "encoder": enc, "label": mol, "n_1a": len(d), "AUC_IHC": np.mean(base),
                     "delta_mean": np.mean(deltas), "delta_sd": np.std(deltas), "n_disc": len(disc),
                     "side_molecular": k, "side_IHC": len(disc) - k, "p": p})
        print(f"R1 {enc:10s} {mol:14s} 1a: IHC {np.mean(base):.3f}, delta {np.mean(deltas):+.3f} ± {np.std(deltas):.3f} | "
              f"1b: n={len(disc)} RNA {k} / IHC {len(disc) - k} p={p:.3g}")

    # R2: ER discordant breakdown (model trained on all concordant ER cases)
    dd = L[["ER_IHC", "ESR1_RNApos", "ER_pct_IHC"]].dropna(subset=["ER_IHC", "ESR1_RNApos"]).join(E, how="inner")
    conc, disc = dd[dd.ER_IHC == dd.ESR1_RNApos], dd[dd.ER_IHC != dd.ESR1_RNApos].copy()
    m = linear_probe().fit(conc[E.columns].values, conc.ER_IHC.astype(int))
    disc["rna_side"] = (m.predict_proba(disc[E.columns].values)[:, 1] >= 0.5).astype(int) == disc.ESR1_RNApos.astype(int)
    disc["direction"] = np.where(disc.ER_IHC == 1, "IHC+/RNA-", "IHC-/RNA+")
    pct = disc.ER_pct_IHC.fillna("unknown")
    low = pct.isin(["<10%"])
    mid = pct.isin(["10-19%", "20-29%", "30-39%"])
    disc["stain"] = np.where(pct.str.contains("Not|unknown"), "unknown", np.where(low, "<10% (ER-low)",
                             np.where(mid, "10-39%", ">=40%")))
    for col in ["direction", "stain"]:
        for grp, sub in disc.groupby(col):
            k = int(sub.rna_side.sum())
            p = binomtest(k, len(sub)).pvalue
            rows.append({"exp": f"R2 ER by {col}", "encoder": enc, "label": grp, "n_disc": len(sub),
                         "side_molecular": k, "side_IHC": len(sub) - k, "p": p})
            print(f"R2 {enc:10s} ER {col}={grp:14s} n={len(sub):3d} RNA {k} / IHC {len(sub) - k} p={p:.3g}")

pd.DataFrame(rows).to_csv(os.path.join(a.out, "robustness_mean.csv"), index=False)
print("saved", os.path.join(a.out, "robustness_mean.csv"))
