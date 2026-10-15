"""Tier 3: overlay ABMIL attention (scripts/attention_maps.py output) on the real H&E slide, to check where attention sits.

For each patient: download its DX slide(s) from the IDC copy (Google Cloud), and draw
  (1) the slide thumbnail, (2) the thumbnail with every patch coloured by its attention percentile within the patient,
  (3) the 6 highest- and 6 lowest-attention patches at 20x (256 px), so one can see what the model looks at.
TRACE coordinates are level-0 pixels of the original slide; IDC holds the same scan (the extraction-equivalence check
reproduced TRACE patch counts), so they map directly. Patch extent at level 0 = 256 px at 20x = 256 * 0.5 / mpp.

Usage: python scripts/attention_overlay.py --att-dir results/attention_v2 --patients TCGA-A2-A4S1,TCGA-AR-A0TX \
           --work /scratch/$USER/overlay --out results/attention_overlay
Output: <out>/<patient>.png and <out>/overlay_index.csv
"""
import argparse
import os
import shutil
import time
from concurrent.futures import ThreadPoolExecutor

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import openslide  # noqa: E402
import pandas as pd  # noqa: E402
import requests  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--att-dir", default="results/attention_v2")
ap.add_argument("--patients", required=True)
ap.add_argument("--idc-objects", default="data/metadata/idc_slides.csv")
ap.add_argument("--work", required=True)
ap.add_argument("--out", default="results/attention_overlay")
ap.add_argument("--thumb", type=int, default=1600, help="thumbnail size (longest side, px)")
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
IDC = pd.read_csv(a.idc_objects)
summary = pd.read_csv(os.path.join(a.att_dir, "summary.csv")).set_index("patient")
GCS = "https://storage.googleapis.com/idc-open-data/"


def fetch(url, path, size, workers=16, step=16 << 20):
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC)
    os.ftruncate(fd, size)

    def rng(b):
        e = min(b + step, size) - 1
        for attempt in range(6):
            try:
                r = requests.get(url, headers={"Range": f"bytes={b}-{e}"}, timeout=(30, 120))
                if r.status_code == 206 and len(r.content) == e - b + 1:
                    os.pwrite(fd, r.content, b)
                    return
            except requests.RequestException:
                pass
            time.sleep(5 * 2 ** attempt)
        raise IOError(f"range {b}-{e} failed")

    try:
        with ThreadPoolExecutor(workers) as ex:
            list(ex.map(rng, range(0, size, step)))
    finally:
        os.close(fd)


def open_slide(barcode, wdir):
    objs = IDC[IDC.slide_barcode == barcode]
    if objs.empty:
        raise IOError("not in IDC listing")
    os.makedirs(wdir, exist_ok=True)
    for o in objs.itertuples():
        fetch(GCS + o.object, os.path.join(wdir, os.path.basename(o.object)), int(o.size))
    return openslide.OpenSlide(os.path.join(wdir, os.path.basename(objs.sort_values("size").object.iloc[-1])))


rows = []
for p in a.patients.split(","):
    t0 = time.time()
    f = os.path.join(a.att_dir, f"{p}.npz")
    if not os.path.exists(f):
        print(f"{p}: no attention file, skipped", flush=True)
        continue
    z = np.load(f, allow_pickle=True)
    C, S, w, slides = z["coords"], z["slide"], z["attention"], [str(x) for x in z["slides"]]
    pct = w.argsort().argsort() / max(len(w) - 1, 1)  # attention percentile within the patient
    info = summary.loc[p]
    for j, name in enumerate(slides):
        barcode = name.split(".")[0]
        wdir = os.path.join(a.work, barcode)
        try:
            sl = open_slide(barcode, wdir)
            mpp = float(sl.properties[openslide.PROPERTY_NAME_MPP_X])
            ext = int(round(256 * 0.5 / mpp))  # patch side at level 0
            W, H = sl.dimensions
            m = S == j
            xy, pc, ww = C[m], pct[m], w[m]
            assert xy[:, 0].max() < W and xy[:, 1].max() < H, "coordinates outside the slide"
            thumb = sl.get_thumbnail((a.thumb, a.thumb)).convert("RGB")
            sc = thumb.width / W
            heat = np.full((thumb.height, thumb.width), np.nan)
            for (x, y), v in zip(xy, pc):
                x0, y0 = int(x * sc), int(y * sc)
                heat[y0:max(y0 + 1, int((y + ext) * sc)), x0:max(x0 + 1, int((x + ext) * sc))] = v
            order = np.argsort(ww)
            top, low = order[::-1][:6], order[:6]

            def crop(k):
                x, y = map(int, xy[k])
                return np.asarray(sl.read_region((x, y), 0, (ext, ext)).convert("RGB").resize((256, 256)))

            fig = plt.figure(figsize=(18, 10))
            gs = fig.add_gridspec(3, 6, height_ratios=[3, 1, 1])
            ax = fig.add_subplot(gs[0, :3])
            ax.imshow(thumb)
            ax.set_title("H&E (DX slide)")
            ax.axis("off")
            ax = fig.add_subplot(gs[0, 3:])
            ax.imshow(thumb)
            ax.imshow(heat, cmap="inferno", vmin=0, vmax=1, alpha=0.55, interpolation="nearest")
            ax.set_title("ABMIL attention percentile (bright = high)")
            ax.axis("off")
            for r, (ks, lab) in enumerate(((top, "highest"), (low, "lowest"))):
                for c, k in enumerate(ks):
                    ax = fig.add_subplot(gs[1 + r, c])
                    ax.imshow(crop(k))
                    ax.set_title(f"{lab} #{c + 1}", fontsize=8)
                    ax.axis("off")
            fig.suptitle(f"{p} | {info.group} | ER IHC {'+' if info.ER_IHC == 1 else '-'} ({info.ER_pct_IHC}) | "
                         f"ESR1 log2 {info.ESR1_log2} | PAM50 {info.PAM50} | P(ER+) {info.P_ER_pos:.2f} -> "
                         f"sides with {info.sides_with} | {barcode}", fontsize=11)
            fig.tight_layout()
            png = os.path.join(a.out, f"{p}{'' if len(slides) == 1 else '_' + str(j + 1)}.png")
            fig.savefig(png, dpi=80)
            plt.close(fig)
            rows.append({"patient": p, "slide": barcode, "mpp": mpp, "patch_ext_level0": ext, "n_patches": int(m.sum()),
                         "png": os.path.basename(png), "seconds": round(time.time() - t0), "status": "ok"})
            print(f"{p} {barcode}: ok ({int(m.sum())} patches, mpp {mpp:.4f}, {time.time() - t0:.0f}s)", flush=True)
        except Exception as e:  # noqa: BLE001
            rows.append({"patient": p, "slide": barcode, "status": f"error: {e}"})
            print(f"{p} {barcode}: FAILED {e}", flush=True)
        finally:
            shutil.rmtree(wdir, ignore_errors=True)
pd.DataFrame(rows).to_csv(os.path.join(a.out, "overlay_index.csv"), index=False)
print("saved", a.out)
