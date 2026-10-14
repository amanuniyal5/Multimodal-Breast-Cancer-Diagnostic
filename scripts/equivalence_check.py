"""Tier 4 extraction-equivalence check: does our pipeline (IDC DICOM copy -> TRIDENT -> UNI v1) reproduce the TRACE UNI
features that every DX result so far is based on? If it does, DX (TRACE) and frozen (ours) features are comparable in 1c.

--select N : pick N patients with exactly one DX slide (seed 0) that are in the IDC listing; print their DX barcodes
--compare  : for those slides, compare our slide mean (data/features/dx_check_uni_v1/<barcode>.npz) with the TRACE
             patient mean (uni_v1_20x_mean.npz):
             - cosine similarity, vs the cosine between DIFFERENT patients (what "not the same slide" looks like);
             - ER: logistic regression (C = 0.01, balanced) trained on TRACE means of all OTHER patients, applied to both
               feature versions of the check slides: Pearson r of P(ER+), agreement of calls, |dP|.

Usage: python scripts/equivalence_check.py --select 50 > /tmp/dx_check.txt
       python scripts/equivalence_check.py --compare --ours data/features/dx_check_uni_v1
Output (compare): results/equivalence_check.csv (per slide) and a printed summary.
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

from common import linear_probe

ap = argparse.ArgumentParser()
ap.add_argument("--feats", default="data/features")
ap.add_argument("--name", default="uni_v1_20x")
ap.add_argument("--idc-objects", default="data/metadata/idc_slides.csv")
ap.add_argument("--labels", default="data/metadata/labels_patient.csv")
ap.add_argument("--select", type=int, default=0)
ap.add_argument("--compare", action="store_true")
ap.add_argument("--ours", default="data/features/dx_check_uni_v1")
ap.add_argument("--out", default="results/equivalence_check.csv")
a = ap.parse_args()

idx = pd.read_csv(os.path.join(a.feats, f"{a.name}_index.csv")).set_index("patient")
z = np.load(os.path.join(a.feats, f"{a.name}_mean.npz"), allow_pickle=True)
TRACE = pd.DataFrame(z["X"], index=z["patients"])

if a.select:
    idc = pd.read_csv(a.idc_objects)
    dx = idc[idc.slide_type == "DX"].drop_duplicates("slide_barcode")
    dx = dx.assign(patient=dx.slide_barcode.str[:12])
    one = set(idx.index[idx.n_slides == 1])  # patient mean == that one slide's mean
    n_dx = dx.groupby("patient").size()
    dx = dx[dx.patient.isin(one) & dx.patient.map(n_dx).eq(1)]
    L = pd.read_csv(a.labels, index_col=0)
    dx = dx[dx.patient.isin(L.index[L.ER_IHC.notna()])]
    pick = dx.sample(a.select, random_state=0).slide_barcode
    print(",".join(sorted(pick)))
    sys.exit()

if a.compare:
    files = sorted(glob.glob(os.path.join(a.ours, "*.npz")))
    ours = pd.DataFrame({os.path.basename(f)[:-4]: np.load(f)["mean"] for f in files}).T
    ours.index = ours.index.str[:12]
    ours = ours[ours.index.isin(TRACE.index)]
    t = TRACE.loc[ours.index]

    def cos(x, y):
        return (x * y).sum(1) / np.linalg.norm(x, axis=1) / np.linalg.norm(y, axis=1)

    same = cos(ours.values, t.values)
    rng = np.random.default_rng(0)
    other = TRACE.drop(ours.index).sample(len(ours), random_state=0).values
    diff = cos(ours.values, other)
    L = pd.read_csv(a.labels, index_col=0)
    train = TRACE.drop(ours.index)
    train = train[train.index.isin(L.index[L.ER_IHC.notna()])]
    clf = linear_probe().fit(train.values, L.loc[train.index, "ER_IHC"].astype(int))
    p_trace, p_ours = clf.predict_proba(t.values)[:, 1], clf.predict_proba(ours.values)[:, 1]
    r = pd.DataFrame({"patient": ours.index, "cos_same_slide": same, "cos_other_patient": diff,
                      "P_ER_trace": p_trace, "P_ER_ours": p_ours})
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    r.to_csv(a.out, index=False)
    print(f"n = {len(r)} DX slides (one DX slide per patient)")
    print(f"cosine(ours, TRACE) same slide: median {np.median(same):.4f}, min {same.min():.4f}")
    print(f"cosine(ours, TRACE) other patient: median {np.median(diff):.4f}, max {diff.max():.4f}")
    print(f"P(ER+) ours vs TRACE: Pearson r {np.corrcoef(p_trace, p_ours)[0, 1]:.3f}, "
          f"median |dP| {np.median(np.abs(p_trace - p_ours)):.3f}, "
          f"same call {int(((p_trace >= 0.5) == (p_ours >= 0.5)).sum())}/{len(r)}")
