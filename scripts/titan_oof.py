"""TITAN linear-probe out-of-fold predictions on the fixed splits, saved in the same format as train_mil.py.

extended_titan.py only writes summary rows; the paired tests in stats_delong_1b.py need patient-level predictions.
Label definitions and experiment groups are those of train_mil.py, so the five model variants can be compared directly.

Usage: python scripts/titan_oof.py [--emb data/features/titan_brca_20x.parquet] [--seeds 5] [--pr-cut gmm|prevmatch]
Output: results/titan_linear{tag}_oof.csv.gz (columns model, exp, label, seed, patient, y, p)
"""
import argparse
import os
import warnings

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.mixture import GaussianMixture
from sklearn.model_selection import cross_val_predict

from common import linear_probe, load_titan, site_folds

warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser()
ap.add_argument("--emb", default="data/features/titan_brca_20x.parquet")
ap.add_argument("--labels", default="data/metadata/labels_patient.csv")
ap.add_argument("--out", default="results")
ap.add_argument("--seeds", type=int, default=5)
ap.add_argument("--pr-cut", default="gmm", choices=["gmm", "prevmatch"])
ap.add_argument("--tag", default="_fixedsplits")
a = ap.parse_args()

E = load_titan(a.emb)
L = pd.read_csv(a.labels, index_col=0)


def gmm_pos(x):  # as in train_mil.py
    x = x.dropna()
    g = GaussianMixture(2, random_state=0).fit(x.values.reshape(-1, 1))
    return pd.Series(g.predict(x.values.reshape(-1, 1)) == np.argmax(g.means_.ravel()), index=x.index).astype(float)


if a.pr_cut == "prevmatch":
    q = L.PGR_log2.quantile(1 - L.PR_IHC.mean())
    L["PGR_RNApos"] = (L.PGR_log2 >= q).astype(float).where(L.PGR_log2.notna())
L["ER_RPPApos"] = gmm_pos(L.ER_RPPA).reindex(L.index)
L["HER2_RPPApos"] = gmm_pos(L.HER2_RPPA).reindex(L.index)
L = L[L.index.isin(E.index)]
print(f"patients with TITAN embeddings and labels: {len(L)}", flush=True)

preds = []
groups = [["ER_IHC", "ESR1_RNApos", "ER_RPPApos"], ["PR_IHC", "PGR_RNApos"],
          ["HER2_STRICT", "ERBB2_AMP", "HER2_RPPApos"]]
for group in groups:
    d = L[group + ["site"]].dropna()
    X = E.loc[d.index].values
    for s in range(a.seeds):
        folds = site_folds(d.index, s)
        msg = []
        for col in group:
            y = d[col].astype(int).values
            p = cross_val_predict(linear_probe(), X, y, cv=folds, method="predict_proba")[:, 1]
            preds.append(pd.DataFrame({"model": "titan_linear", "exp": "1a", "label": col, "seed": s,
                                       "patient": d.index, "y": d[col].values, "p": p}))
            msg.append(f"{col}={roc_auc_score(y, p):.3f}")
        print(f"1a n={len(d)} seed {s}: " + ", ".join(msg), flush=True)

for ihc, mol in [("ER_IHC", "ESR1_RNApos"), ("ER_IHC", "ER_RPPApos"), ("PR_IHC", "PGR_RNApos"),
                 ("HER2_STRICT", "ERBB2_AMP")]:
    d = L[[ihc, mol]].dropna()
    conc, disc = d[d[ihc] == d[mol]], d[d[ihc] != d[mol]]
    m = linear_probe().fit(E.loc[conc.index].values, conc[ihc].astype(int))
    p = m.predict_proba(E.loc[disc.index].values)[:, 1]
    k = int(((p >= 0.5).astype(int) == disc[mol].astype(int).values).sum())
    preds.append(pd.DataFrame({"model": "titan_linear", "exp": "1b", "label": f"{ihc} vs {mol}", "seed": -1,
                               "patient": disc.index, "y": disc[mol].values, "p": p}))
    print(f"1b {ihc} vs {mol}: n={len(disc)} molecular {k} / IHC {len(disc) - k}", flush=True)

out = os.path.join(a.out, f"titan_linear{a.tag}_oof.csv.gz")
pd.concat(preds).to_csv(out, index=False)
print("saved", out)
