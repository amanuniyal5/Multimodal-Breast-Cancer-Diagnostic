"""Encoder registry (configs/encoders.json) and patient-level feature loading."""
import glob
import json
import os

import numpy as np
import pandas as pd

from common import load_titan

ENC = json.load(open("configs/encoders.json"))


def available(feats="data/features"):
    """encoders whose patient-level features exist (ResNet-50: built from the per-slide files on first use)"""
    out = []
    for e, c in ENC.items():
        f = os.path.join(feats, f"{c['file']}_mean.npz")
        if e == "titan":
            f = os.path.join(feats, "titan_brca_20x.parquet")
        if os.path.exists(f) or (e == "resnet50" and glob.glob(os.path.join(feats, "dx_resnet50", "*.npz"))):
            out.append(e)
    return out


def slides_to_patients(slide_dir, out_npz):
    """per-slide extract_frozen.py outputs (DX) -> patient mean weighted by patch count, like TRACE's patient means"""
    rows = {}
    for f in sorted(glob.glob(os.path.join(slide_dir, "*.npz"))):
        z = np.load(f)
        p = os.path.basename(f)[:12]
        s, n = rows.get(p, (0.0, 0))
        rows[p] = (s + z["mean"].astype(np.float64) * int(z["n_patches"]), n + int(z["n_patches"]))
    pts = sorted(p for p in rows if rows[p][1] > 0)
    X = np.stack([rows[p][0] / rows[p][1] for p in pts]).astype(np.float32)
    np.savez(out_npz, patients=np.array(pts, dtype=str), X=X, n=np.array([rows[p][1] for p in pts]))
    return pts


def load_means(enc, feats="data/features", std=False):
    """DataFrame patients x d (patient-level mean features); std=True appends per-feature patch SD where stored"""
    c = ENC[enc]
    if enc == "titan":
        return load_titan(os.path.join(feats, "titan_brca_20x.parquet"))
    f = os.path.join(feats, f"{c['file']}_mean.npz")
    if enc == "resnet50":  # rebuilt every time: the extraction may still be adding slides
        slides_to_patients(os.path.join(feats, "dx_resnet50"), f)
    z = np.load(f, allow_pickle=True)
    X = pd.DataFrame(z["X"], index=z["patients"].astype(str))
    if std:
        if "S" not in z.files:
            return None
        X = pd.concat([X, pd.DataFrame(z["S"], index=X.index).add_prefix("s")], axis=1)
    X.columns = [str(k) for k in X.columns]
    return X


def label(enc):
    return ENC[enc]["label"] + (" †" if ENC[enc]["tcga_pretrained"] else "")
