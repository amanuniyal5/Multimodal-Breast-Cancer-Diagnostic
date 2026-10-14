"""Experiment 1e: do the two frozen slides of one tumour (top TS and bottom BS, cut from the same vial + portion,
bracketing the sequenced tissue) give the same prediction?

Model: ER linear probe on mean-pooled UNI v1 features (C = 0.01, balanced), trained on the primary TS slides of the
patients outside the test fold (fixed site-grouped folds, 5 seeds), applied to BOTH slides of each test pair; P(ER+)
averaged over seeds. Done for the IHC label and the RNA label.
Agreement over pairs: ICC(2,1) of logit P (two-way random, absolute agreement), Pearson r of P, Cohen's kappa and
% flipped calls at 0.5, AUC on TS vs BS; bootstrap 95% CIs over pairs (2,000). Also: feature cosine TS vs BS (same
tumour vs different patients), and whether |dP| grows with the difference in % tumour cells (Spearman).

Usage: python scripts/frozen_1e.py [--frozen data/features/frozen_uni_v1] [--seeds 5] [--boot 2000]
Output: results/frozen_1e.csv (summary), results/frozen_1e_pairs.csv (per pair)
"""
import argparse
import os
import warnings

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import cohen_kappa_score, roc_auc_score

from common import linear_probe, load_frozen_means, site_folds

warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser()
ap.add_argument("--frozen", default="data/features/frozen_uni_v1")
ap.add_argument("--manifest", default="data/metadata/frozen_manifest.csv")
ap.add_argument("--labels", default="data/metadata/labels_patient.csv")
ap.add_argument("--seeds", type=int, default=5)
ap.add_argument("--boot", type=int, default=2000)
ap.add_argument("--out", default="results")
a = ap.parse_args()

man = pd.read_csv(a.manifest)
man["pair_id"] = man.pair_id.fillna("")
F = load_frozen_means(a.frozen, man[man.primary_TS | (man.pair_id != "")])
feat = [c for c in F.columns if isinstance(c, int)]
L = pd.read_csv(a.labels, index_col=0)

prim = F[F.primary_TS].set_index("patient")
pairs = F[F.pair_id != ""]
ts = pairs[pairs.section_type == "TS"].drop_duplicates("pair_id").set_index("pair_id")
bs = pairs[pairs.section_type == "BS"].drop_duplicates("pair_id").set_index("pair_id")
pid = ts.index.intersection(bs.index)
ts, bs = ts.loc[pid], bs.loc[pid]
pair_patient = ts.patient
print(f"primary TS slides with features: {len(prim)}; complete TS+BS pairs: {len(pid)}", flush=True)


def icc21(x, y):
    """ICC(2,1), two-way random effects, absolute agreement, single rater (Shrout & Fleiss)."""
    M = np.column_stack([x, y])
    n, k = M.shape
    gm = M.mean()
    msr = k * ((M.mean(1) - gm) ** 2).sum() / (n - 1)
    msc = n * ((M.mean(0) - gm) ** 2).sum() / (k - 1)
    sse = ((M - M.mean(1, keepdims=True) - M.mean(0, keepdims=True) + gm) ** 2).sum()
    mse = sse / ((n - 1) * (k - 1))
    return (msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n)


def logit(p):
    p = np.clip(p, 1e-4, 1 - 1e-4)
    return np.log(p / (1 - p))


def stats(pt, pb, y):
    ct, cb = (pt >= 0.5).astype(int), (pb >= 0.5).astype(int)
    out = {"ICC_logit": icc21(logit(pt), logit(pb)), "pearson_P": np.corrcoef(pt, pb)[0, 1],
           "kappa": cohen_kappa_score(ct, cb), "pct_flipped": 100 * (ct != cb).mean(),
           "median_abs_dP": np.median(np.abs(pt - pb))}
    if len(set(y)) == 2:
        out.update(AUC_TS=roc_auc_score(y, pt), AUC_BS=roc_auc_score(y, pb))
    return out


rows, per_pair = [], []
rng = np.random.default_rng(0)
for lab in ["ER_IHC", "ESR1_RNApos"]:
    tr_pts = prim.index[prim.index.isin(L.index[L[lab].notna() & L.site.notna()])]
    test = pid[pair_patient.isin(L.index[L[lab].notna()]).values]
    PT, PB = pd.Series(0.0, index=test), pd.Series(0.0, index=test)
    y_all = L.loc[tr_pts, lab].astype(int).values
    for s in range(a.seeds):
        for tr, te in site_folds(tr_pts, s):
            m = linear_probe().fit(prim.loc[tr_pts[tr], feat].values, y_all[tr])
            te_pairs = test[pair_patient.loc[test].isin(tr_pts[te]).values]
            if len(te_pairs):
                PT[te_pairs] += m.predict_proba(ts.loc[te_pairs, feat].values)[:, 1] / a.seeds
                PB[te_pairs] += m.predict_proba(bs.loc[te_pairs, feat].values)[:, 1] / a.seeds
    keep = test[pair_patient.loc[test].isin(tr_pts).values]  # pairs whose patient got an out-of-fold prediction
    pt, pb = PT[keep].values, PB[keep].values
    y = L.loc[pair_patient.loc[keep], lab].astype(int).values
    est = stats(pt, pb, y)
    bs_ = [stats(pt[i], pb[i], y[i]) for i in (rng.integers(0, len(keep), len(keep)) for _ in range(a.boot))]
    bs_ = pd.DataFrame(bs_)
    row = {"label": lab, "n_pairs": len(keep), "pos": int(y.sum())}
    for k, v in est.items():
        row[k] = v
        row[k + "_lo"], row[k + "_hi"] = np.nanpercentile(bs_[k], [2.5, 97.5])
    dt = (ts.loc[keep, "percent_tumor_cells"].values - bs.loc[keep, "percent_tumor_cells"].values)
    ok = ~np.isnan(dt)
    rho, p_rho = spearmanr(np.abs(dt[ok]), np.abs(pt - pb)[ok])
    row.update(spearman_absdP_vs_absdTumour=rho, spearman_p=p_rho)
    rows.append(row)
    per_pair.append(pd.DataFrame({"label": lab, "pair_id": keep, "patient": pair_patient.loc[keep].values, "y": y,
                                  "P_TS": pt, "P_BS": pb, "tumour_TS": ts.loc[keep, "percent_tumor_cells"].values,
                                  "tumour_BS": bs.loc[keep, "percent_tumor_cells"].values}))
    print(f"1e {lab}: {len(keep)} pairs | ICC(logit P) {row['ICC_logit']:.3f} [{row['ICC_logit_lo']:.3f}, "
          f"{row['ICC_logit_hi']:.3f}] | r(P) {row['pearson_P']:.3f} | kappa {row['kappa']:.3f} "
          f"[{row['kappa_lo']:.3f}, {row['kappa_hi']:.3f}] | flipped {row['pct_flipped']:.1f}% "
          f"[{row['pct_flipped_lo']:.1f}, {row['pct_flipped_hi']:.1f}] | AUC TS {row.get('AUC_TS', np.nan):.3f} "
          f"BS {row.get('AUC_BS', np.nan):.3f} | |dP| vs |d tumour %| rho {rho:.2f} (p={p_rho:.2g})", flush=True)

# feature-level similarity: same tumour vs different patients
X_t, X_b = ts[feat].values, bs[feat].values
cos = (X_t * X_b).sum(1) / np.linalg.norm(X_t, axis=1) / np.linalg.norm(X_b, axis=1)
perm = rng.permutation(len(pid))
cos_other = (X_t * X_b[perm]).sum(1) / np.linalg.norm(X_t, axis=1) / np.linalg.norm(X_b[perm], axis=1)
rows.append({"label": "feature cosine TS vs BS", "n_pairs": len(pid), "median_same_tumour": np.median(cos),
             "median_other_patient": np.median(cos_other)})
print(f"feature cosine TS vs BS: same tumour median {np.median(cos):.3f}, other patient {np.median(cos_other):.3f}")

os.makedirs(a.out, exist_ok=True)
pd.DataFrame(rows).to_csv(os.path.join(a.out, "frozen_1e.csv"), index=False)
pd.concat(per_pair).to_csv(os.path.join(a.out, "frozen_1e_pairs.csv"), index=False)
print("saved", os.path.join(a.out, "frozen_1e.csv"), os.path.join(a.out, "frozen_1e_pairs.csv"))
