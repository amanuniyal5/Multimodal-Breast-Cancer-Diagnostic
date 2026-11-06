"""X1: manifest of DX (FFPE) slides in IDC for new patch encoders (ResNet-50), in the format extract_frozen.py reads.

Patients (those in data/splits/folds.csv with a DX slide in IDC) are put in a fixed random order (seed 0) and cut into
batches of 25 patients, so the first batches form a pre-set random subset: the ResNet-50 comparison subset is the first
--subset patients of this order (all their DX slides), fixed before any result is seen.

Usage: python scripts/make_dx_manifest.py --out data/metadata/dx_manifest.csv --subset 450
"""
import argparse

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--idc-objects", default="data/metadata/idc_slides.csv")
ap.add_argument("--folds", default="data/splits/folds.csv")
ap.add_argument("--out", default="data/metadata/dx_manifest.csv")
ap.add_argument("--per-batch", type=int, default=25)
ap.add_argument("--subset", type=int, default=450)
a = ap.parse_args()

x = pd.read_csv(a.idc_objects)
x = x[x.slide_type == "DX"]
s = x.groupby("slide_barcode", as_index=False)["size"].sum()
s["patient"] = s.slide_barcode.str[:12]
pts = sorted(set(s.patient) & set(pd.read_csv(a.folds).patient))
order = {p: i for i, p in enumerate(np.random.default_rng(0).permutation(pts))}
s = s[s.patient.isin(order)].copy()
s["patient_order"] = s.patient.map(order)
s["batch"] = s.patient_order // a.per_batch
s["subset"] = s.patient_order < a.subset
s["section_type"] = "DX"
s = s.sort_values(["patient_order", "slide_barcode"]).reset_index(drop=True)
s.to_csv(a.out, index=False)
sub = s[s.subset]
print(f"{len(s)} DX slides, {s.patient.nunique()} patients, {s.batch.nunique()} batches, {s['size'].sum() / 1e9:.0f} GB")
print(f"subset: {sub.patient.nunique()} patients, {len(sub)} slides, batches 0-{sub.batch.max()}, "
      f"{sub['size'].sum() / 1e9:.0f} GB")
