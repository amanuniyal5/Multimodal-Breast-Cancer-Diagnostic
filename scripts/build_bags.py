"""Turn TRACE parquet chunks (rows of up to 256 patches) into patient-level features.

Outputs (in --out):
  <name>_mean.npz      exact mean of ALL patch features per patient (DX slides only) -> mean-pool baselines
  <name>_cap<K>.npy    float16 array [n_patients, K, D]: a uniform random sample of K patches per patient
                       (all patches if fewer; unused rows are zero) -> ABMIL
  <name>_index.csv     patient, row, n_patches_total, n_kept, n_slides

Only DX (FFPE diagnostic) slides are used. Sampling = vectorised reservoir sampling (seed 0), written straight
into a disk memmap so RAM stays small. Put --work on fast local disk (/scratch); the script copies results to --out.

Usage: python scripts/build_bags.py --parquet-dir /scratch/$USER/trace/data/TCGA-BRCA/uni_v1/20x \
           --work /scratch/$USER/bags --out ~/project/data/features --name uni_v1_20x --cap 4096
"""
import argparse
import glob
import os
import shutil
import time

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

ap = argparse.ArgumentParser()
ap.add_argument("--parquet-dir", required=True)
ap.add_argument("--work", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--name", required=True)
ap.add_argument("--cap", type=int, default=4096)
ap.add_argument("--max-patients", type=int, default=1200)
ap.add_argument("--max-files", type=int, default=0, help="debug: only read the first N parquet files")
a = ap.parse_args()
os.makedirs(a.work, exist_ok=True)
os.makedirs(a.out, exist_ok=True)
files = sorted(glob.glob(os.path.join(a.parquet_dir, "*.parquet")))
if a.max_files:
    files = files[:a.max_files]
print(f"{len(files)} parquet files", flush=True)

rng = np.random.default_rng(0)
row_of, sums, seen, filled, slides = {}, {}, {}, {}, {}
mm, D = None, None
t0 = time.time()
for fi, f in enumerate(files):
    pf = pq.ParquetFile(f)
    for batch in pf.iter_batches(batch_size=16, columns=["case_submitter_id", "svs_filename", "embeddings"]):
        emb = batch.column("embeddings")
        lens = emb.value_lengths().to_numpy()
        flat = emb.flatten().flatten().to_numpy(zero_copy_only=False)
        if D is None:
            D = int(flat.size // lens.sum())
            mm = np.lib.format.open_memmap(os.path.join(a.work, "tmp_cap.npy"), mode="w+", dtype=np.float16,
                                           shape=(a.max_patients, a.cap, D))
            print("embedding dim", D, flush=True)
        X = flat.reshape(-1, D)
        start = 0
        for case, name, n in zip(batch.column("case_submitter_id").to_pylist(),
                                 batch.column("svs_filename").to_pylist(), lens):
            chunk = X[start:start + n]
            start += n
            if "-DX" not in name or n == 0:
                continue
            if case not in row_of:
                row_of[case] = len(row_of)
                sums[case], seen[case], filled[case], slides[case] = np.zeros(D), 0, 0, set()
            r = row_of[case]
            sums[case] += chunk.sum(0, dtype=np.float64)
            slides[case].add(name)
            k = min(a.cap - filled[case], n)  # fill free slots first
            if k > 0:
                mm[r, filled[case]:filled[case] + k] = chunk[:k]
                filled[case] += k
            rest = chunk[k:]
            seen[case] += k
            if len(rest):  # reservoir: item number t (1-based) replaces slot j ~ U[0, t) if j < cap
                t = seen[case] + 1 + np.arange(len(rest))
                j = (rng.random(len(rest)) * t).astype(np.int64)
                keep = j < a.cap
                if keep.any():
                    mm[r, j[keep]] = rest[keep]  # later duplicates win, as in sequential sampling
                seen[case] += len(rest)
    print(f"[{fi + 1}/{len(files)}] {os.path.basename(f)} patients {len(row_of)} ({time.time() - t0:.0f}s)", flush=True)

patients = sorted(row_of)
mean = np.stack([sums[p] / seen[p] for p in patients]).astype(np.float32)
np.savez(os.path.join(a.out, f"{a.name}_mean.npz"), patients=np.array(patients, dtype=str), X=mean)
final = np.lib.format.open_memmap(os.path.join(a.work, f"{a.name}_cap{a.cap}.npy"), mode="w+", dtype=np.float16,
                                  shape=(len(patients), a.cap, D))
rows = []
for i, p in enumerate(patients):
    final[i] = mm[row_of[p]]
    rows.append({"patient": p, "row": i, "n_patches_total": seen[p], "n_kept": filled[p], "n_slides": len(slides[p])})
final.flush()
del final, mm
os.remove(os.path.join(a.work, "tmp_cap.npy"))
pd.DataFrame(rows).to_csv(os.path.join(a.out, f"{a.name}_index.csv"), index=False)
shutil.copy(os.path.join(a.work, f"{a.name}_cap{a.cap}.npy"), os.path.join(a.out, f"{a.name}_cap{a.cap}.npy"))
print(f"done: {len(patients)} patients, D={D}, median patches {int(np.median([seen[p] for p in patients]))}, "
      f"{time.time() - t0:.0f}s", flush=True)
