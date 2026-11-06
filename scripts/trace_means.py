"""X1: patient-level mean (and per-patch SD) features from TRACE parquet, streamed one file at a time.

Unlike build_bags.py this keeps no patch bags and never holds the whole encoder folder on disk: each parquet file is
downloaded to --work (a few files ahead, in parallel), reduced to per-patient sums, then deleted. DX slides only, the
same rule as build_bags.py, so the means are directly comparable (check: rerun on conch_v1/20x and compare with
conch_v1_20x_mean.npz).

Outputs (in --out):
  <name>_mean.npz    patients, X (mean of all patch features, float32), S (per-feature SD over patches), n (patches)
  <name>_index.csv   patient, n_patches_total, n_slides

--slide-level: rows are one embedding per slide (TRACE chief/10x, gigapath_slide/20x): patient X = mean over the
patient's DX slides; S and n are omitted.

Usage: python scripts/trace_means.py --path data/TCGA-BRCA/uni_v2/20x --work /scratch/$USER/tm --out data/features \
           --name uni_v2_20x
"""
import argparse
import os
import shutil
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from huggingface_hub import HfApi, hf_hub_download

ap = argparse.ArgumentParser()
ap.add_argument("--repo", default="IUCompPath/TRACE")
ap.add_argument("--path", required=True, help="folder inside the dataset repo, e.g. data/TCGA-BRCA/uni_v2/20x")
ap.add_argument("--work", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--name", required=True)
ap.add_argument("--ahead", type=int, default=4, help="parquet files downloaded in parallel ahead of processing")
ap.add_argument("--slide-level", action="store_true")
ap.add_argument("--max-files", type=int, default=0, help="debug: only the first N files")
a = ap.parse_args()
os.makedirs(a.work, exist_ok=True)
os.makedirs(a.out, exist_ok=True)
t0 = time.time()


def log(*x):
    print(f"[{time.time() - t0:6.0f}s]", *x, flush=True)


files = sorted(f for f in HfApi().list_repo_files(a.repo, repo_type="dataset")
               if f.startswith(a.path.rstrip("/") + "/") and f.endswith(".parquet"))
if a.max_files:
    files = files[:a.max_files]
log(f"{len(files)} parquet files under {a.path}")


def is_list(t):
    return pa.types.is_list(t) or pa.types.is_large_list(t) or pa.types.is_fixed_size_list(t)


def fetch(f):
    for attempt in range(6):
        try:
            return hf_hub_download(a.repo, f, repo_type="dataset", local_dir=a.work)
        except Exception as e:  # noqa: BLE001  (HF/network hiccups: back off and retry)
            log(f"retry {attempt + 1} {f}: {e}")
            time.sleep(30 * (attempt + 1))
    raise IOError(f"could not download {f}")


sums, sq, cnt, slides, D = {}, {}, {}, {}, None
pool = ThreadPoolExecutor(a.ahead)
futs = [pool.submit(fetch, f) for f in files[:a.ahead]]
for i, f in enumerate(files):
    path = futs[i].result()
    if i + a.ahead < len(files):
        futs.append(pool.submit(fetch, files[i + a.ahead]))
    for batch in pq.ParquetFile(path).iter_batches(batch_size=16,
                                                    columns=["case_submitter_id", "svs_filename", "embeddings"]):
        emb = batch.column("embeddings")
        nested = is_list(emb.type.value_type)  # list<list<float>> = patches x D; list<float> = one vector per row
        lens = emb.value_lengths().to_numpy() if nested else np.ones(len(emb), dtype=np.int64)
        flat = emb
        while is_list(flat.type):
            flat = flat.flatten()
        flat = flat.to_numpy(zero_copy_only=False)
        if D is None:
            D = int(flat.size // lens.sum())
            log("embedding dim", D)
        X = flat.reshape(-1, D).astype(np.float64)
        start = 0
        for case, name, n in zip(batch.column("case_submitter_id").to_pylist(),
                                 batch.column("svs_filename").to_pylist(), lens):
            chunk = X[start:start + n]
            start += n
            if "-DX" not in name or n == 0:
                continue
            if case not in sums:
                sums[case], sq[case], cnt[case], slides[case] = np.zeros(D), np.zeros(D), 0, {}
            if a.slide_level:  # one embedding per slide (possibly repeated rows): average per slide, then per patient
                slides[case].setdefault(name, []).append(chunk.mean(0))
                continue
            sums[case] += chunk.sum(0)
            sq[case] += (chunk ** 2).sum(0)
            cnt[case] += n
            slides[case][name] = True
    os.remove(path)
    log(f"[{i + 1}/{len(files)}] {os.path.basename(f)} patients {len(sums)}")

pts = sorted(sums)
if a.slide_level:
    X = np.stack([np.mean([np.mean(v, 0) for v in slides[p].values()], 0) for p in pts]).astype(np.float32)
    np.savez(os.path.join(a.out, f"{a.name}_mean.npz"), patients=np.array(pts, dtype=str), X=X)
else:
    X = np.stack([sums[p] / cnt[p] for p in pts])
    S = np.sqrt(np.maximum(np.stack([sq[p] / cnt[p] for p in pts]) - X ** 2, 0))
    np.savez(os.path.join(a.out, f"{a.name}_mean.npz"), patients=np.array(pts, dtype=str), X=X.astype(np.float32),
             S=S.astype(np.float32), n=np.array([cnt[p] for p in pts]))
pd.DataFrame({"patient": pts, "n_patches_total": [cnt[p] for p in pts],
              "n_slides": [len(slides[p]) for p in pts]}).to_csv(os.path.join(a.out, f"{a.name}_index.csv"), index=False)
shutil.rmtree(a.work, ignore_errors=True)
log(f"done: {len(pts)} patients, D={D}")
print("TRACE_MEANS_OK", flush=True)
