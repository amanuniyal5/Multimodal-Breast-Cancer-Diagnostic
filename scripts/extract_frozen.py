"""Tier 4: download frozen (TS/BS) or DX slides from GDC, run TRIDENT (tissue segmentation -> 256 px patches at 20x
-> patch encoder), and keep per-slide features only; the raw slide is deleted right after encoding.

Per slide, written to --out:
  <slide_barcode>.npz   mean (exact mean of all patch features, float32), bag (uniform random sample of --cap patches,
                        float16, seed 0), coords (of the bag patches, level-0 pixels), n_patches
and one row per slide (seconds for download / segmentation / patching / encoding, patch count, mpp) in
--out/timing_<job>.csv. Slides already present in --out are skipped, so a job can be resubmitted after a time-out.

--encoder uni_v1        UNI v1 (gated: needs Hugging Face access + `huggingface-cli login` on Ada)
--encoder vitl_random   same architecture as UNI v1 (ViT-L/16, 224 px) with random weights: timing only

Usage: python scripts/extract_frozen.py --manifest data/metadata/frozen_manifest.csv --batch 0 \
           --work /scratch/$USER/frozen --out data/features/frozen_uni_v1 --encoder uni_v1
       (--slides TCGA-A7-A0DC-01A-01-TS1,... selects slides by barcode instead of --batch; --limit N for tests)
"""
import argparse
import os
import shutil
import sys
import time

import h5py
import numpy as np
import pandas as pd
import requests
import torch

ap = argparse.ArgumentParser()
ap.add_argument("--manifest", default="data/metadata/frozen_manifest.csv")
ap.add_argument("--batch", default=None, help="manifest batch: 3, a range 0-4, or a list 0,5,7")
ap.add_argument("--slides", default="", help="comma-separated slide barcodes (overrides --batch)")
ap.add_argument("--limit", type=int, default=0)
ap.add_argument("--work", required=True, help="node-local scratch for raw slides and TRIDENT outputs")
ap.add_argument("--out", required=True)
ap.add_argument("--encoder", default="uni_v1", choices=["uni_v1", "vitl_random"])
ap.add_argument("--segmenter", default="hest")
ap.add_argument("--mag", type=int, default=20)
ap.add_argument("--patch-size", type=int, default=256)
ap.add_argument("--cap", type=int, default=4096)
ap.add_argument("--batch-limit", type=int, default=256, help="patches per GPU batch")
ap.add_argument("--keep-h5", action="store_true", help="also keep TRIDENT's full feature .h5 (large)")
a = ap.parse_args()

from trident import OpenSlideWSI  # noqa: E402  (after argparse so --help works without TRIDENT)
from trident.patch_encoder_models import encoder_factory  # noqa: E402
from trident.segmentation_models import segmentation_model_factory  # noqa: E402

os.makedirs(a.work, exist_ok=True)
os.makedirs(a.out, exist_ok=True)
dev = "cuda:0" if torch.cuda.is_available() else "cpu"
job = os.environ.get("SLURM_JOB_ID", "local")
NCPU = int(os.environ.get("SLURM_CPUS_PER_TASK", 4))
if a.encoder == "uni_v1" and dev == "cpu":
    sys.exit("no usable GPU (check the torch build vs the driver): refusing to run UNI on CPU")
t_start = time.time()


def log(*x):
    print(f"[{time.time() - t_start:6.0f}s]", *x, flush=True)


m = pd.read_csv(a.manifest)
if a.slides:
    m = m[m.slide_barcode.isin(a.slides.split(","))]
elif a.batch is not None:
    sel = set()
    for part in str(a.batch).split(","):
        lo, _, hi = part.partition("-")
        sel |= set(range(int(lo), int(hi or lo) + 1))
    m = m[m.batch.isin(sel)]
if a.limit:
    m = m.head(a.limit)
todo = m[~m.slide_barcode.apply(lambda s: os.path.exists(os.path.join(a.out, f"{s}.npz")))]
log(f"{len(m)} slides selected, {len(todo)} to do ({todo['size'].sum() / 1e9:.1f} GB); device {dev} "
    f"{torch.cuda.get_device_name(0) if dev != 'cpu' else ''}")

seg_model = segmentation_model_factory(a.segmenter)
if a.encoder == "uni_v1":
    encoder = encoder_factory("uni_v1")
else:  # timing surrogate: UNI v1 is timm vit_large_patch16_224 (init_values=1e-5, no class token pooling change)
    import timm
    from trident.patch_encoder_models import CustomInferenceEncoder
    from torchvision import transforms as T
    net = timm.create_model("vit_large_patch16_224", pretrained=False, init_values=1e-5, num_classes=0,
                            dynamic_img_size=True)
    tf = T.Compose([T.Resize(224), T.ToTensor(), T.Normalize((0.485, 0.456, 0.406), (0.229, 0.224, 0.225))])
    encoder = CustomInferenceEncoder(enc_name="vitl_random", model=net, transforms=tf, precision=torch.float16)
log("encoder ready:", a.encoder)


def download(file_id, path, size):
    with requests.get(f"https://api.gdc.cancer.gov/data/{file_id}", stream=True, timeout=600) as r:
        r.raise_for_status()
        with open(path, "wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 22):
                f.write(chunk)
    got = os.path.getsize(path)
    if got != size:  # GDC downloads can truncate silently (see HANDOFF gotchas)
        raise IOError(f"truncated download: {got} of {size} bytes")


rows = []
timing_csv = os.path.join(a.out, f"timing_{job}.csv")
for _, s in todo.iterrows():
    wdir = os.path.join(a.work, s.slide_barcode)
    os.makedirs(wdir, exist_ok=True)
    svs = os.path.join(wdir, s.file_name)
    row = {"slide_barcode": s.slide_barcode, "size_GB": s["size"] / 1e9, "node": os.uname().nodename}
    try:
        t = time.time()
        for attempt in range(3):
            try:
                download(s.file_id, svs, int(s["size"]))
                break
            except Exception as e:  # noqa: BLE001
                log(f"  download attempt {attempt + 1} failed: {e}")
                if attempt == 2:
                    raise
        row["t_download"] = time.time() - t

        t = time.time()
        # TRIDENT sizes its DataLoaders from os.cpu_count() (40 on a node), not the job allocation: cap them
        slide = OpenSlideWSI(slide_path=svs, lazy_init=False, max_workers=NCPU)
        row["mpp"] = getattr(slide, "mpp", np.nan)
        slide.segment_tissue(segmentation_model=seg_model, target_mag=10, job_dir=wdir, device=dev,
                             num_workers=NCPU)
        row["t_segment"] = time.time() - t

        t = time.time()
        coords_path = slide.extract_tissue_coords(target_mag=a.mag, patch_size=a.patch_size, save_coords=wdir)
        row["t_patch"] = time.time() - t

        t = time.time()
        feats_path = slide.extract_patch_features(patch_encoder=encoder, coords_path=coords_path,
                                                  save_features=os.path.join(wdir, "features"), device=dev,
                                                  batch_limit=a.batch_limit)
        if dev != "cpu":
            torch.cuda.synchronize()
        row["t_encode"] = time.time() - t
        slide.release() if hasattr(slide, "release") else None

        with h5py.File(feats_path, "r") as h:
            X = h["features"][:]
            C = h["coords"][:]
        n = len(X)
        rng = np.random.default_rng(0)
        keep = np.sort(rng.choice(n, a.cap, replace=False)) if n > a.cap else np.arange(n)
        np.savez(os.path.join(a.out, f"{s.slide_barcode}.npz"), mean=X.astype(np.float64).mean(0).astype(np.float32),
                 bag=X[keep].astype(np.float16), coords=C[keep], n_patches=n)
        if a.keep_h5:
            shutil.copy(feats_path, os.path.join(a.out, f"{s.slide_barcode}.h5"))
        row.update(n_patches=n, status="ok")
        log(f"{s.slide_barcode}: {n} patches, {row['size_GB']:.2f} GB | download {row['t_download']:.0f}s, "
            f"segment {row['t_segment']:.0f}s, patch {row['t_patch']:.0f}s, encode {row['t_encode']:.0f}s")
    except Exception as e:  # noqa: BLE001  (log and continue with the next slide)
        row["status"] = f"error: {e}"
        log(f"{s.slide_barcode}: FAILED {e}")
    finally:
        shutil.rmtree(wdir, ignore_errors=True)  # never keep raw slides
    rows.append(row)
    pd.DataFrame(rows).to_csv(timing_csv, index=False)

ok = [r for r in rows if r.get("status") == "ok"]
if ok:
    tot = pd.DataFrame(ok)[["t_download", "t_segment", "t_patch", "t_encode"]]
    log(f"done: {len(ok)}/{len(rows)} ok; mean seconds per slide " +
        ", ".join(f"{k[2:]} {v:.0f}" for k, v in tot.mean().items()) + f"; total {tot.sum(1).mean():.0f}")
if rows and not ok:
    sys.exit("no slide succeeded")  # non-zero exit, so dependent jobs (afterok) do not start
