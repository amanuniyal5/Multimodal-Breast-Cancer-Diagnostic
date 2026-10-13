"""Extended baselines on TITAN slide embeddings (laptop CPU, ~10 min).

E1  repeated site-grouped CV (5 seeds) for every label incl. the shuffled control
E2  1a paired IHC vs RNA vs RPPA-protein label (same patients, same folds, per seed)
E3  1b "which side" for ER, PR and HER2 (train on concordant, score discordant)
E4  RNA cut-off sensitivity for ER: GMM / prevalence-matched / confident-only
E5  1d HER2-low: 0 vs 1+ and 0 vs low (1+ or 2+ FISH-), H&E model vs ERBB2 mRNA alone

Usage: python scripts/extended_titan.py [--emb ...] [--labels ...] [--out results] [--seeds 5]
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
ap.add_argument("--emb", default="data/features/titan_brca_20x.parquet")
ap.add_argument("--labels", default="data/metadata/labels_patient.csv")
ap.add_argument("--out", default="results")
ap.add_argument("--seeds", type=int, default=5)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
SEEDS = range(a.seeds)

E = load_titan(a.emb)
L = pd.read_csv(a.labels, index_col=0)
c = pd.read_csv("data/metadata/clinical_patient_brca.txt", sep="\t", skiprows=[1, 2], dtype=str).set_index(
    "bcr_patient_barcode")


def gmm(x, log=True):
    x = x.dropna()
    v = (np.log2(x + 1) if log else x).values.reshape(-1, 1)
    g = GaussianMixture(2, random_state=0).fit(v)
    hi = np.argmax(g.means_.ravel())
    post = g.predict_proba(v)[:, hi]
    return pd.Series(post, index=x.index)


# extra labels
L["ER_RPPApos"] = (gmm(L.ER_RPPA, log=False) > 0.5).astype(float).reindex(L.index)
L["HER2_RPPApos"] = (gmm(L.HER2_RPPA, log=False) > 0.5).astype(float).reindex(L.index)
post = gmm(np.power(2, L.ESR1_log2) - 1).reindex(L.index)
ihc_prev = L.ER_IHC.mean()
q = L.ESR1_log2.quantile(1 - ihc_prev)
L["ESR1_prevmatch"] = (L.ESR1_log2 >= q).astype(float).where(L.ESR1_log2.notna())
L["ESR1_confident"] = np.where(post > 0.9, 1.0, np.where(post < 0.1, 0.0, np.nan))
L["HER2low_1plus"] = L.HER2_IHCscore.map({"0": 0, "1+": 1})
fish = c.her2_fish_status.reindex(L.index)
low = (L.HER2_IHCscore == "1+") | ((L.HER2_IHCscore == "2+") & (fish == "Negative"))
L["HER2_0_vs_low"] = np.where(L.HER2_IHCscore == "0", 0.0, np.where(low, 1.0, np.nan))
L["ER_IHC_shuffled"] = np.random.default_rng(0).permutation(L.ER_IHC.values)


def oof(cols, seed):
    """site-grouped out-of-fold probabilities for each label in cols, on the patients that have all of them"""
    d = L[cols + ["site"]].dropna().join(E, how="inner")
    X, g = d[E.columns].values, d.site.values
    ref = d[cols[0]].astype(int).values
    folds = site_folds(d.index, seed)  # fixed site-grouped folds (data/splits)
    out = {}
    for col in cols:
        y = d[col].astype(int).values
        out[col] = cross_val_predict(linear_probe(), X, y, cv=folds, method="predict_proba")[:, 1]
    return d, out


rows = []
# E1: repeated CV per label
for col in ["ER_IHC", "ESR1_RNApos", "ER_RPPApos", "PR_IHC", "PGR_RNApos", "HER2_STRICT", "ERBB2_AMP",
            "HER2_RPPApos", "HER2low_1plus", "HER2_0_vs_low", "ER_IHC_shuffled"]:
    aucs = []
    for s in SEEDS:
        d, o = oof([col], s)
        aucs.append(roc_auc_score(d[col].astype(int), o[col]))
    rows.append({"exp": "E1 repeated CV", "label": col, "n": len(d), "pos": int(d[col].sum()),
                 "AUC_mean": np.mean(aucs), "AUC_sd": np.std(aucs)})
    print(f"E1 {col:16s} n={len(d):4d} AUC {np.mean(aucs):.3f} ± {np.std(aucs):.3f}")


# E2/E4: paired label comparisons, per seed, with paired bootstrap on seed 0
def paired(cols, tag):
    deltas = {c2: [] for c2 in cols[1:]}
    base = []
    for s in SEEDS:
        d, o = oof(cols, s)
        a0 = roc_auc_score(d[cols[0]].astype(int), o[cols[0]])
        base.append(a0)
        for c2 in cols[1:]:
            deltas[c2].append(roc_auc_score(d[c2].astype(int), o[c2]) - a0)
    d, o = oof(cols, 0)
    r = np.random.default_rng(0)
    for c2 in cols[1:]:
        y0, y2 = d[cols[0]].astype(int).values, d[c2].astype(int).values
        bs = []
        for _ in range(2000):
            i = r.integers(0, len(d), len(d))
            if len(set(y0[i])) < 2 or len(set(y2[i])) < 2:
                continue
            bs.append(roc_auc_score(y2[i], o[c2][i]) - roc_auc_score(y0[i], o[cols[0]][i]))
        lo, hi = np.percentile(bs, [2.5, 97.5])
        rows.append({"exp": tag, "label": f"{c2} - {cols[0]}", "n": len(d), "AUC_ref": np.mean(base),
                     "delta_mean": np.mean(deltas[c2]), "delta_sd": np.std(deltas[c2]), "CI_lo": lo, "CI_hi": hi})
        print(f"{tag} {c2} vs {cols[0]}: n={len(d)} ref AUC {np.mean(base):.3f}; "
              f"delta {np.mean(deltas[c2]):+.3f} ± {np.std(deltas[c2]):.3f} (seed-0 95% CI {lo:+.3f}, {hi:+.3f})")


paired(["ER_IHC", "ESR1_RNApos", "ER_RPPApos"], "E2 1a ER")
paired(["PR_IHC", "PGR_RNApos"], "E2 1a PR")
paired(["HER2_STRICT", "ERBB2_AMP", "HER2_RPPApos"], "E2 1a HER2")
paired(["ER_IHC", "ESR1_prevmatch"], "E4 ER cut-off prevalence-matched")
paired(["ER_IHC", "ESR1_confident"], "E4 ER cut-off confident-only")


# E3: which side on discordant cases
def side(ihc, mol, tag):
    d = L[[ihc, mol]].dropna().join(E, how="inner")
    conc, disc = d[d[ihc] == d[mol]], d[d[ihc] != d[mol]]
    m = linear_probe().fit(conc[E.columns].values, conc[ihc].astype(int))
    pred = (m.predict_proba(disc[E.columns].values)[:, 1] >= 0.5).astype(int)
    k = int((pred == disc[mol].astype(int).values).sum())
    p = binomtest(k, len(disc)).pvalue
    rows.append({"exp": tag, "label": f"{ihc} vs {mol}", "n": len(disc), "side_molecular": k,
                 "side_IHC": len(disc) - k, "p": p})
    print(f"{tag} {ihc} vs {mol}: discordant n={len(disc)}, sides with molecular {k}, IHC {len(disc) - k}, p={p:.3g}")


side("ER_IHC", "ESR1_RNApos", "E3 1b")
side("ER_IHC", "ER_RPPApos", "E3 1b")
side("PR_IHC", "PGR_RNApos", "E3 1b")
side("HER2_STRICT", "ERBB2_AMP", "E3 1b")
side("ER_IHC", "ESR1_prevmatch", "E4 1b cut-off")
side("ER_IHC", "ESR1_confident", "E4 1b cut-off")

# E5: HER2-low, H&E model vs ERBB2 mRNA alone (same patients)
for col in ["HER2low_1plus", "HER2_0_vs_low"]:
    d = L[[col, "ERBB2_log2", "site"]].dropna().join(E, how="inner")
    y = d[col].astype(int).values
    he = [roc_auc_score(y, cross_val_predict(linear_probe(), d[E.columns].values, y,
                                             cv=site_folds(d.index, s), method="predict_proba")[:, 1]) for s in SEEDS]
    rna = roc_auc_score(y, d.ERBB2_log2)
    rows.append({"exp": "E5 1d", "label": col, "n": len(d), "pos": int(y.sum()), "AUC_HE_mean": np.mean(he),
                 "AUC_HE_sd": np.std(he), "AUC_ERBB2_mRNA": rna})
    print(f"E5 {col}: n={len(d)} (pos {y.sum()}), H&E AUC {np.mean(he):.3f} ± {np.std(he):.3f}, ERBB2 mRNA AUC {rna:.3f}")

pd.DataFrame(rows).to_csv(os.path.join(a.out, "extended_titan.csv"), index=False)
print("saved", os.path.join(a.out, "extended_titan.csv"))
