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

Opt-in extras (X1 ABMIL on a second encoder, X5b maps); both stream in the same single pass:
  --bags K       also write <name>_cap<K>.npy (fp16, patients x K x D, uniform reservoir sample, seed 0, unused rows 0),
                 <name>_cap<K>_coords.npy (int32, -1 padded) and <name>_cap<K>_slide.npy (DX slide number per patch);
                 the index CSV gains row / n_kept. Reservoirs live in RAM (about 7 GB for K = 2,048, D = 1,536).
  --keep FILE    patient list (one per line or a CSV with a 'patient' column): ALL their DX patches are written to
                 --keep-dir/<patient>.npz (X fp16, coords int32, slide int, slides) for patch-level maps.

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
ap.add_argument("--bags", type=int, default=0, help="also keep a K-patch reservoir sample per patient (0 = off)")
ap.add_argument("--keep", default="", help="patients whose full DX patch features + coords are written")
ap.add_argument("--keep-dir", default="")
a = ap.parse_args()
KEEP = set()
if a.keep:
    k = pd.read_csv(a.keep) if a.keep.endswith(".csv") else pd.DataFrame({"patient": open(a.keep).read().split()})
    KEEP = set(k.patient.astype(str))
    os.makedirs(a.keep_dir, exist_ok=True)
RES = {}  # patient -> [features (K, D) fp16, coords (K, 2), slide (K,), seen]
KEPT = {p: {"X": [], "C": [], "S": [], "slides": []} for p in KEEP}
RNG = np.random.default_rng(0)
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
    cols = ["case_submitter_id", "svs_filename", "embeddings"] + (["coords"] if (a.bags or KEEP) else [])
    for batch in pq.ParquetFile(path).iter_batches(batch_size=16, columns=cols):
        emb = batch.column("embeddings")
        CRD = None
        if "coords" in cols:
            c = batch.column("coords")
            while is_list(c.type):
                c = c.flatten()
            CRD = c.to_numpy(zero_copy_only=False).reshape(-1, 2).astype(np.int32)
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
            sid = list(slides[case]).index(name)
            crd = CRD[start - n:start] if CRD is not None else None
            if case in KEEP:
                KEPT[case]["X"].append(chunk.astype(np.float16))
                KEPT[case]["C"].append(crd)
                KEPT[case]["S"].append(np.full(n, sid, dtype=np.int16))
                if name not in KEPT[case]["slides"]:
                    KEPT[case]["slides"].append(name)
            if a.bags:  # reservoir sampling (Algorithm R, vectorised; later writes win, as in the sequential version)
                if case not in RES:
                    RES[case] = [np.zeros((a.bags, D), np.float16), np.full((a.bags, 2), -1, np.int32),
                                 np.full(a.bags, -1, np.int16), 0]
                r = RES[case]
                j = r[3] + np.arange(n)
                fill = j < a.bags
                r[0][j[fill]], r[1][j[fill]], r[2][j[fill]] = chunk[fill], crd[fill], sid
                rest = np.where(~fill)[0]
                if len(rest):
                    t = RNG.integers(0, j[rest] + 1)
                    ok = t < a.bags
                    r[0][t[ok]], r[1][t[ok]], r[2][t[ok]] = chunk[rest[ok]], crd[rest[ok]], sid
                r[3] += n
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
I = pd.DataFrame({"patient": pts, "n_patches_total": [cnt[p] for p in pts], "n_slides": [len(slides[p]) for p in pts]})
if a.bags:
    I["row"] = np.arange(len(pts))
    I["n_kept"] = [min(cnt[p], a.bags) for p in pts]
    B = np.lib.format.open_memmap(os.path.join(a.out, f"{a.name}_cap{a.bags}.npy"), mode="w+", dtype=np.float16,
                                  shape=(len(pts), a.bags, D))
    for i, p in enumerate(pts):
        B[i] = RES[p][0]
    B.flush()
    del B
    np.save(os.path.join(a.out, f"{a.name}_cap{a.bags}_coords.npy"), np.stack([RES[p][1] for p in pts]))
    np.save(os.path.join(a.out, f"{a.name}_cap{a.bags}_slide.npy"), np.stack([RES[p][2] for p in pts]))
    log(f"bags written: {len(pts)} x {a.bags} x {D}")
I.to_csv(os.path.join(a.out, f"{a.name}_index.csv"), index=False)
for p, k in KEPT.items():
    if k["X"]:
        np.savez(os.path.join(a.keep_dir, f"{p}.npz"), X=np.concatenate(k["X"]), coords=np.concatenate(k["C"]),
                 slide=np.concatenate(k["S"]), slides=np.array(k["slides"]))
        log(f"kept {p}: {sum(len(x) for x in k['X'])} patches, {len(k['slides'])} slides")
    else:
        log(f"keep: no DX patches found for {p}")
shutil.rmtree(a.work, ignore_errors=True)
log(f"done: {len(pts)} patients, D={D}")
print("TRACE_MEANS_OK", flush=True)
