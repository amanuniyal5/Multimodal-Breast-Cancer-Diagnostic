"""Experiment 1c: does the answer to "which truth does the slide see?" change between the FFPE diagnostic slide (DX) and
the frozen slide cut next to the sequenced tissue (TS)?

Same patients, same fixed site-grouped folds (data/splits), same linear probe (mean-pooled UNI v1, C = 0.01, balanced):
  AUC[DX, IHC], AUC[DX, RNA], AUC[TS, IHC], AUC[TS, RNA]   (probabilities averaged over 5 seeds per patient)
  interaction = (AUC[TS, mol] - AUC[TS, IHC]) - (AUC[DX, mol] - AUC[DX, IHC]); paired bootstrap over patients (2,000)
  1b on TS: train on concordant patients (TS features), score discordant ones, binomial test; compared with DX on the
  same discordant patients.
  % tumour cells (biospecimen table, TS slide): interaction within the lower / upper half.
DX features = TRACE UNI v1 (patient mean over DX slides); TS features = our extraction (primary TS slide per patient).
The extraction-equivalence check (results/equivalence_check.csv) showed both come from the same pipeline.

Usage: python scripts/frozen_1c.py [--frozen data/features/frozen_uni_v1] [--seeds 5] [--boot 2000]
Output: results/frozen_1c.csv (summary), results/frozen_1c_oof.csv.gz (patient-level predictions)
"""
import argparse
import os
import warnings

import numpy as np
import pandas as pd
from scipy.stats import binomtest
from sklearn.metrics import roc_auc_score
from sklearn.mixture import GaussianMixture

from common import linear_probe, load_frozen_means, site_folds

warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser()
ap.add_argument("--feats", default="data/features")
ap.add_argument("--name", default="uni_v1_20x")
ap.add_argument("--frozen", default="data/features/frozen_uni_v1")
ap.add_argument("--manifest", default="data/metadata/frozen_manifest.csv")
ap.add_argument("--labels", default="data/metadata/labels_patient.csv")
ap.add_argument("--seeds", type=int, default=5)
ap.add_argument("--boot", type=int, default=2000)
ap.add_argument("--out", default="results")
a = ap.parse_args()

z = np.load(os.path.join(a.feats, f"{a.name}_mean.npz"), allow_pickle=True)
DX = pd.DataFrame(z["X"], index=z["patients"])
man = pd.read_csv(a.manifest)
TSm = load_frozen_means(a.frozen, man[man.primary_TS])  # one primary TS slide per patient
TS = TSm.set_index("patient")[[c for c in TSm.columns if isinstance(c, int)]]
tumour = TSm.set_index("patient").percent_tumor_cells
L = pd.read_csv(a.labels, index_col=0)


def gmm_pos(x):  # as in train_mil.py
    x = x.dropna()
    g = GaussianMixture(2, random_state=0).fit(x.values.reshape(-1, 1))
    return pd.Series(g.predict(x.values.reshape(-1, 1)) == np.argmax(g.means_.ravel()), index=x.index).astype(float)


L["ER_RPPApos"] = gmm_pos(L.ER_RPPA).reindex(L.index)
q = L.PGR_log2.quantile(1 - L.PR_IHC.mean())
L["PGR_prevmatch"] = (L.PGR_log2 >= q).astype(float).where(L.PGR_log2.notna())
both = DX.index.intersection(TS.index)
print(f"patients with DX features {len(DX)}, primary TS features {len(TS)}, both {len(both)}", flush=True)

COMPARISONS = [("ER", "ER_IHC", "ESR1_RNApos"), ("ER", "ER_IHC", "ER_RPPApos"), ("PR", "PR_IHC", "PGR_prevmatch"),
               ("HER2", "HER2_STRICT", "ERBB2_AMP")]
rows, preds = [], []
rng = np.random.default_rng(0)
for marker, ihc, mol in COMPARISONS:
    d = L.loc[L.index.intersection(both), [ihc, mol, "site"]].dropna()
    pts = d.index
    P = {}
    for src, F in (("DX", DX), ("TS", TS)):
        X = F.loc[pts].values
        for lab in (ihc, mol):
            y = d[lab].astype(int).values
            pr = np.zeros(len(pts))
            for s in range(a.seeds):
                for tr, te in site_folds(pts, s):
                    pr[te] += linear_probe().fit(X[tr], y[tr]).predict_proba(X[te])[:, 1] / a.seeds
            P[src, lab] = pr
            preds.append(pd.DataFrame({"marker": marker, "slide": src, "label": lab, "patient": pts, "y": y, "p": pr}))
    Y = {lab: d[lab].astype(int).values for lab in (ihc, mol)}

    def aucs(i):
        return {k: roc_auc_score(Y[k[1]][i], P[k][i]) for k in P}

    def inter(A):
        return (A["TS", mol] - A["TS", ihc]) - (A["DX", mol] - A["DX", ihc])

    A = aucs(np.arange(len(pts)))
    bs = []
    for _ in range(a.boot):
        i = rng.integers(0, len(pts), len(pts))
        if min(len(set(Y[ihc][i])), len(set(Y[mol][i]))) < 2:
            continue
        B = aucs(i)
        bs.append([inter(B), B["DX", mol] - B["DX", ihc], B["TS", mol] - B["TS", ihc]])
    bs = np.array(bs)
    lo, hi = np.percentile(bs, [2.5, 97.5], axis=0)
    p_int = 2 * min((bs[:, 0] <= 0).mean(), (bs[:, 0] >= 0).mean())
    row = {"marker": marker, "comparison": f"{mol} vs {ihc}", "n": len(pts), "pos_ihc": int(Y[ihc].sum()),
           "pos_mol": int(Y[mol].sum()), "AUC_DX_ihc": A["DX", ihc], "AUC_DX_mol": A["DX", mol],
           "AUC_TS_ihc": A["TS", ihc], "AUC_TS_mol": A["TS", mol],
           "delta_DX": A["DX", mol] - A["DX", ihc], "delta_DX_lo": lo[1], "delta_DX_hi": hi[1],
           "delta_TS": A["TS", mol] - A["TS", ihc], "delta_TS_lo": lo[2], "delta_TS_hi": hi[2],
           "interaction": inter(A), "interaction_lo": lo[0], "interaction_hi": hi[0], "interaction_p_boot": p_int}
    # % tumour cells: interaction within lower / upper half of the TS slide's tumour content
    t = tumour.reindex(pts)
    med = t.median()
    for half, msk in (("low", (t < med).values), ("high", (t >= med).values)):
        if msk.sum() > 30 and min(Y[ihc][msk].sum(), Y[mol][msk].sum(), (1 - Y[ihc][msk]).sum()) >= 5:
            Ah = aucs(np.where(msk)[0])
            row[f"interaction_tumour_{half}"] = inter(Ah)
            row[f"n_tumour_{half}"] = int(msk.sum())
    row["tumour_median_pct"] = med
    rows.append(row)
    print(f"1c {marker} {mol} vs {ihc}: n={len(pts)} | DX {A['DX', ihc]:.3f}->{A['DX', mol]:.3f} "
          f"(d {row['delta_DX']:+.3f} [{lo[1]:+.3f}, {hi[1]:+.3f}]) | TS {A['TS', ihc]:.3f}->{A['TS', mol]:.3f} "
          f"(d {row['delta_TS']:+.3f} [{lo[2]:+.3f}, {hi[2]:+.3f}]) | interaction {row['interaction']:+.3f} "
          f"[{lo[0]:+.3f}, {hi[0]:+.3f}] p={p_int:.3f}", flush=True)

# 1b on frozen slides: which side on discordant patients (same patients for DX and TS)
for ihc, mol in [("ER_IHC", "ESR1_RNApos"), ("HER2_STRICT", "ERBB2_AMP"), ("PR_IHC", "PGR_prevmatch")]:
    d = L.loc[L.index.intersection(both), [ihc, mol]].dropna()
    conc, disc = d[d[ihc] == d[mol]], d[d[ihc] != d[mol]]
    out = {"marker": ihc.split("_")[0], "comparison": f"1b {ihc} vs {mol}", "n": len(disc)}
    for src, F in (("DX", DX), ("TS", TS)):
        p = linear_probe().fit(F.loc[conc.index].values, conc[ihc].astype(int)).predict_proba(F.loc[disc.index].values)[:, 1]
        k = int(((p >= 0.5).astype(int) == disc[mol].astype(int).values).sum())
        out.update({f"side_mol_{src}": k, f"side_ihc_{src}": len(disc) - k, f"p_{src}": binomtest(k, len(disc)).pvalue})
        preds.append(pd.DataFrame({"marker": out["marker"], "slide": src, "label": f"1b {ihc} vs {mol}",
                                   "patient": disc.index, "y": disc[mol].values, "p": p}))
    rows.append(out)
    print(f"1b {ihc} vs {mol}: n={len(disc)} | DX molecular {out['side_mol_DX']}:{out['side_ihc_DX']} "
          f"(p={out['p_DX']:.3g}) | TS molecular {out['side_mol_TS']}:{out['side_ihc_TS']} (p={out['p_TS']:.3g})", flush=True)

os.makedirs(a.out, exist_ok=True)
pd.DataFrame(rows).to_csv(os.path.join(a.out, "frozen_1c.csv"), index=False)
pd.concat(preds).to_csv(os.path.join(a.out, "frozen_1c_oof.csv.gz"), index=False)
print("saved results/frozen_1c.csv, results/frozen_1c_oof.csv.gz")
