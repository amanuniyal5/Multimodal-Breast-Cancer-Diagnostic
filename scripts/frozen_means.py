"""Slide-level UNI means of the frozen TS/BS slides (data/features/frozen_uni_v1/<barcode>.npz, key 'mean') in one small
file, so experiment G (DX <-> frozen transfer) and X6 can run on the laptop.

Output: data/features/frozen_uni_v1_slides.npz with slides, patients, section (TS/BS), primary_TS, pair_id, n_patches, X
Usage: python scripts/frozen_means.py [--dir data/features/frozen_uni_v1] [--manifest data/metadata/frozen_manifest.csv]
"""
import argparse
import os

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--dir", default="data/features/frozen_uni_v1")
ap.add_argument("--manifest", default="data/metadata/frozen_manifest.csv")
ap.add_argument("--out", default="data/features/frozen_uni_v1_slides.npz")
a = ap.parse_args()
M = pd.read_csv(a.manifest).drop_duplicates("slide_barcode")
rows, X = [], []
for _, r in M.iterrows():
    f = os.path.join(a.dir, f"{r.slide_barcode}.npz")
    if not os.path.exists(f):
        continue
    z = np.load(f)
    if int(z["n_patches"]) == 0:
        continue
    X.append(z["mean"].astype(np.float32))
    rows.append((r.slide_barcode, r.patient, r.section_type, bool(r.primary_TS), str(r.pair_id), int(z["n_patches"])))
R = pd.DataFrame(rows, columns=["slide", "patient", "section", "primary_TS", "pair_id", "n_patches"])
np.savez(a.out, X=np.stack(X), slides=R.slide.to_numpy(dtype=str), patients=R.patient.to_numpy(dtype=str),
         section=R.section.to_numpy(dtype=str), pair_id=R.pair_id.to_numpy(dtype=str),
         primary_TS=R.primary_TS.to_numpy(dtype=bool), n_patches=R.n_patches.to_numpy(dtype=np.int64))
print(f"{len(R)} slides ({(R.section == 'TS').sum()} TS, {(R.section == 'BS').sum()} BS), "
      f"{R.patient.nunique()} patients -> {a.out}")
print("FROZEN_MEANS_OK")
