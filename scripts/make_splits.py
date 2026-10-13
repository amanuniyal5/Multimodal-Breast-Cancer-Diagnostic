"""Create fixed, library-version-independent, site-grouped 5-fold splits for all experiments.

Why: scikit-learn 1.5 and 1.9 shuffle StratifiedGroupKFold differently, so the same seed gave different folds on the
laptop and on Ada. These splits are written once to data/splits/folds.csv and every script reads them.

Algorithm (per seed): shuffle the 40 tissue-source sites with numpy's legacy RandomState(seed) (stable across numpy
versions), then assign sites largest-first to the fold that currently has the fewest patients. All patients of a site
share one fold, so every experiment (any subset of patients) is automatically site-grouped. A seed is accepted only if
every fold contains at least 2 patients of each class for every task; otherwise a salt is added and it retries.

Output: data/splits/folds.csv with columns patient, site, fold_s0 ... fold_s4
Usage:  python scripts/make_splits.py [--labels data/metadata/labels_patient.csv] [--seeds 5]
"""
import argparse
import os

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--labels", default="data/metadata/labels_patient.csv")
ap.add_argument("--clinical", default="data/metadata/clinical_patient_brca.txt")
ap.add_argument("--out", default="data/splits/folds.csv")
ap.add_argument("--seeds", type=int, default=5)
ap.add_argument("--k", type=int, default=5)
a = ap.parse_args()

L = pd.read_csv(a.labels, index_col=0)
c = pd.read_csv(a.clinical, sep="\t", skiprows=[1, 2], dtype=str).set_index("bcr_patient_barcode")
fish = c.her2_fish_status.reindex(L.index)
tasks = {col: L[col] for col in ["ER_IHC", "ESR1_RNApos", "PR_IHC", "PGR_RNApos", "HER2_STRICT", "ERBB2_AMP"]}
tasks["HER2low_1plus"] = L.HER2_IHCscore.map({"0": 0, "1+": 1})
low = (L.HER2_IHCscore == "1+") | ((L.HER2_IHCscore == "2+") & (fish == "Negative"))
tasks["HER2_0_vs_low"] = pd.Series(np.where(L.HER2_IHCscore == "0", 0.0, np.where(low, 1.0, np.nan)), index=L.index)

sizes = L.site.value_counts()
out = pd.DataFrame({"patient": L.index, "site": L.site.values}).set_index("patient")


her2_0 = L[L.HER2_IHCscore == "0"].site.value_counts()  # HER2-0 patients come from only ~7 sites


def assign(seed):
    rs = np.random.RandomState(seed)
    order = list(sizes.index[rs.permutation(len(sizes))])
    load, h0, fold_of = np.zeros(a.k), np.zeros(a.k), {}
    # phase 1: sites holding HER2-0 patients, most first -> fold with fewest HER2-0 (ties: fewest patients, random)
    first = sorted([s for s in order if her2_0.get(s, 0) > 0], key=lambda s: -her2_0[s])
    for s in first:
        f = int(np.lexsort((rs.random(a.k), load, h0))[0])
        fold_of[s] = f
        load[f] += sizes[s]
        h0[f] += her2_0[s]
    # phase 2: remaining sites largest first (random among equal sizes) -> fold with fewest patients
    rest = sorted([s for s in order if s not in fold_of], key=lambda s: -sizes[s])
    for s in rest:
        f = int(np.argmin(load))
        fold_of[s] = f
        load[f] += sizes[s]
    return L.site.map(fold_of)


for s in range(a.seeds):
    for salt in range(1000):
        f = assign(s * 1000 + salt)
        ok = True
        for name, y in tasks.items():
            t = pd.crosstab(f[y.notna()], y.dropna())
            if t.shape != (a.k, 2) or t.values.min() < 2:
                ok = False
                break
        if ok:
            break
    out[f"fold_s{s}"] = f.values
    print(f"seed {s}: salt {salt}, fold sizes {np.bincount(f, minlength=a.k).tolist()}")

os.makedirs(os.path.dirname(a.out), exist_ok=True)
out.reset_index().to_csv(a.out, index=False)
# report the smallest class count per fold for each task (seed 0)
for name, y in tasks.items():
    t = pd.crosstab(out.fold_s0[y.notna()], y.dropna())
    print(f"  {name:15s} min class per fold (seed 0): {int(t.values.min())}")
print("saved", a.out)
