"""Tier 4: manifest of primary-tumour frozen slides (TS = top, BS = bottom) for experiments 1c and 1e.

Keeps open-access TCGA-BRCA slides with sample code 01 (primary tumour) and a TS/BS section, for patients in
labels_patient.csv; adds % tumour cells etc. from the biospecimen slide table; marks TS+BS pairs cut from the same
vial + portion (they bracket the sequenced tissue); assigns download batches of roughly equal size for a job array.

Usage: python scripts/make_frozen_manifest.py [--batch-gb 15]
Output: data/metadata/frozen_manifest.csv (one row per slide, with column `batch`)
"""
import argparse

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--slides", default="data/metadata/gdc_slides.csv")
ap.add_argument("--biospecimen", default="data/metadata/biospecimen_slide_brca.txt")
ap.add_argument("--labels", default="data/metadata/labels_patient.csv")
ap.add_argument("--batch-gb", type=float, default=15.0, help="target raw size per download batch")
ap.add_argument("--out", default="data/metadata/frozen_manifest.csv")
a = ap.parse_args()

s = pd.read_csv(a.slides)
bc = s.file_name.str.split(".").str[0]  # e.g. TCGA-A7-A0DC-01A-01-TS1
parts = bc.str.split("-", expand=True)
s = s.assign(slide_barcode=bc, patient=s.case, sample=parts[3], portion=parts[4], section=parts[5])
s["sample_type"] = s["sample"].str[:2]
s["section_type"] = s.section.str[:2]
s["vial_portion"] = s.patient + "-" + s["sample"] + "-" + s.portion
print(f"all slides: {len(s)}; by section type: {s.section_type.value_counts().to_dict()}")

L = pd.read_csv(a.labels, index_col=0)
f = s[(s.section_type.isin(["TS", "BS"])) & (s.sample_type == "01") & (s.access == "open") & s.patient.isin(L.index)]
print(f"primary-tumour TS/BS slides of labelled patients: {len(f)} ({f['size'].sum() / 1e9:.0f} GB); "
      f"{f.section_type.value_counts().to_dict()}; patients {f.patient.nunique()}")

b = pd.read_csv(a.biospecimen, sep="\t", skiprows=[1], dtype=str)
keep = ["percent_tumor_cells", "percent_tumor_nuclei", "percent_stromal_cells", "percent_necrosis",
        "percent_lymphocyte_infiltration", "percent_normal_cells", "section_location"]
b = b.drop_duplicates("bcr_slide_barcode").set_index("bcr_slide_barcode")[keep]
for k in keep[:-1]:
    b[k] = pd.to_numeric(b[k], errors="coerce")
f = f.join(b, on="slide_barcode")
print(f"with % tumour cells: {f.percent_tumor_cells.notna().sum()} of {len(f)}; "
      f"median TS {f[f.section_type == 'TS'].percent_tumor_cells.median():.0f}%, "
      f"BS {f[f.section_type == 'BS'].percent_tumor_cells.median():.0f}%")

# 1e pairs: a TS and a BS from the same vial + portion (first slide of each type if there are several)
ts = f[f.section_type == "TS"].sort_values("section").groupby("vial_portion").slide_barcode.first()
bs = f[f.section_type == "BS"].sort_values("section").groupby("vial_portion").slide_barcode.first()
pairs = ts.index.intersection(bs.index)
f["pair_id"] = np.where(f.vial_portion.isin(pairs) & f.slide_barcode.isin(set(ts[pairs]) | set(bs[pairs])),
                        f.vial_portion, "")
print(f"TS+BS pairs from the same vial+portion: {len(pairs)} ({pd.Index(pairs).str[:12].nunique()} patients)")

# 1c: one TS slide per patient (first by section name) is the primary frozen slide
f["primary_TS"] = f.slide_barcode.isin(f[f.section_type == "TS"].sort_values(["sample", "portion", "section"])
                                       .groupby("patient").slide_barcode.first())
print(f"patients with a primary TS slide: {int(f.primary_TS.sum())}")

# batches: greedy, in a fixed order (pairs and primary TS first so the essential slides finish early)
f = f.assign(_prio=np.where(f.primary_TS | (f.pair_id != ""), 0, 1)).sort_values(["_prio", "patient", "section"])
batch, acc, cur = [], 0.0, 0
for sz in f["size"]:
    if acc + sz / 1e9 > a.batch_gb and acc > 0:
        cur, acc = cur + 1, 0.0
    batch.append(cur)
    acc += sz / 1e9
f["batch"] = batch
cols = ["file_id", "file_name", "size", "slide_barcode", "patient", "section_type", "vial_portion", "pair_id",
        "primary_TS", "batch"] + keep
f[cols].to_csv(a.out, index=False)
print(f"batches: {f.batch.nunique()} (~{a.batch_gb:.0f} GB each; essential slides in batches 0-"
      f"{f[f._prio == 0].batch.max()}); saved {a.out}")
