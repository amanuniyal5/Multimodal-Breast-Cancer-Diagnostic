"""Download TRACE TITAN slide embeddings for TCGA-BRCA DX slides (~7 MB, 1,125 slides x 768-d).

Source: https://huggingface.co/datasets/IUCompPath/TRACE (licence not stated; TITAN weights are CC-BY-NC-ND).
Academic use only - do NOT commit the parquet to git (data/features/ is git-ignored).

Usage: python scripts/fetch_titan.py [--out data/features/titan_brca_20x.parquet]
"""
import argparse
import os

import requests

URL = ("https://huggingface.co/datasets/IUCompPath/TRACE/resolve/main/"
       "data/TCGA-BRCA/titan/20x/part-00000.parquet")

ap = argparse.ArgumentParser()
ap.add_argument("--out", default="data/features/titan_brca_20x.parquet")
a = ap.parse_args()
os.makedirs(os.path.dirname(a.out), exist_ok=True)
r = requests.get(URL, timeout=300)
r.raise_for_status()
with open(a.out, "wb") as f:
    f.write(r.content)
print(f"saved {len(r.content) / 1e6:.1f} MB -> {a.out}")
