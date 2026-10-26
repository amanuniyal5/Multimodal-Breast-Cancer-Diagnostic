"""Experiments 1c and 1e with the ABMIL head (frozen_1c.py / frozen_1e.py use the mean-pool linear probe).

Same patients, labels, fixed site-grouped folds and comparisons as frozen_1c.py / frozen_1e.py; the classifier is the
gated-attention ABMIL of train_mil.py with the same settings (15 epochs, 2,048 random patches per training bag, Adam
lr 2e-4, weight decay 1e-4, positive-class weight), 3 seeds, probabilities averaged over seeds per patient.
Bags: DX = TRACE UNI v1 (uni_v1_20x_cap4096.npy, ≤ 4,096 patches per patient); TS/BS = our extraction
(frozen_uni_v1/<slide>.npz "bag", ≤ 4,096 patches per slide). Same encoder and patching (results/equivalence_check.csv).

Training runs as separate tasks so each task holds only one set of bags in memory:
  --slide DX|TS --part er       1c ER: ESR1 RNA vs IHC, ER RPPA vs IHC; 1b ER (train concordant, score discordant)
  --slide DX|TS --part prher2   1c PR (prevalence-matched PGR) and HER2 (amplification vs strict); 1b PR and HER2
  --slide TS    --part 1e       1e: ER model (IHC label, RNA label) trained on primary TS slides of out-of-fold
                                patients, applied to both slides of each TS+BS pair
  each writes results/frozen_abmil_raw_<slide>_<part>.csv.gz (per seed)
  --analyse                     seed-averaged predictions in the format of frozen_1c_oof.csv.gz / frozen_1e_pairs.csv
                                (frozen_abmil_1c_oof.csv.gz, frozen_abmil_1e_pairs.csv) and summaries
                                (frozen_abmil_1c.csv: AUCs, ΔAUC, interaction with paired bootstrap CI;
                                frozen_abmil_1e.csv: ICC, r, kappa, % flipped with bootstrap CIs)

Usage: python scripts/frozen_abmil.py --slide TS --part er
       python scripts/frozen_abmil.py --analyse
"""
import argparse
import os
import time
import warnings

import numpy as np
import pandas as pd
from scipy.stats import binomtest
from sklearn.metrics import cohen_kappa_score, roc_auc_score
from sklearn.mixture import GaussianMixture

from common import site_folds

warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser()
ap.add_argument("--feats", default="data/features")
ap.add_argument("--name", default="uni_v1_20x")
ap.add_argument("--cap", type=int, default=4096)
ap.add_argument("--frozen", default="data/features/frozen_uni_v1")
ap.add_argument("--manifest", default="data/metadata/frozen_manifest.csv")
ap.add_argument("--labels", default="data/metadata/labels_patient.csv")
ap.add_argument("--slide", choices=["DX", "TS"])
ap.add_argument("--part", choices=["er", "prher2", "1e"])
ap.add_argument("--analyse", action="store_true")
ap.add_argument("--seeds", type=int, default=3)
ap.add_argument("--epochs", type=int, default=15)
ap.add_argument("--train-patches", type=int, default=2048)
ap.add_argument("--lr", type=float, default=2e-4)
ap.add_argument("--boot", type=int, default=2000)
ap.add_argument("--out", default="results")
a = ap.parse_args()
t_start = time.time()


def log(*x):
    print(f"[{time.time() - t_start:6.0f}s]", *x, flush=True)


def gmm_pos(x):  # as in train_mil.py / frozen_1c.py
    x = x.dropna()
    g = GaussianMixture(2, random_state=0).fit(x.values.reshape(-1, 1))
    return pd.Series(g.predict(x.values.reshape(-1, 1)) == np.argmax(g.means_.ravel()), index=x.index).astype(float)


L = pd.read_csv(a.labels, index_col=0)
L["ER_RPPApos"] = gmm_pos(L.ER_RPPA).reindex(L.index)
q = L.PGR_log2.quantile(1 - L.PR_IHC.mean())
L["PGR_prevmatch"] = (L.PGR_log2 >= q).astype(float).where(L.PGR_log2.notna())
COMPARISONS = {"er": [("ER", "ER_IHC", "ESR1_RNApos"), ("ER", "ER_IHC", "ER_RPPApos")],
               "prher2": [("PR", "PR_IHC", "PGR_prevmatch"), ("HER2", "HER2_STRICT", "ERBB2_AMP")]}
PAIRS_1B = {"er": [("ER_IHC", "ESR1_RNApos")], "prher2": [("HER2_STRICT", "ERBB2_AMP"), ("PR_IHC", "PGR_prevmatch")]}

man = pd.read_csv(a.manifest).drop_duplicates("slide_barcode")
man["pair_id"] = man.pair_id.fillna("")
man = man[man.slide_barcode.apply(lambda s: os.path.exists(os.path.join(a.frozen, f"{s}.npz")))]
prim = man[man.primary_TS].drop_duplicates("patient").set_index("patient")  # one primary TS slide per patient
idx = pd.read_csv(os.path.join(a.feats, f"{a.name}_index.csv")).set_index("patient")
both = idx.index.intersection(prim.index)

if not a.analyse:
    import torch
    import torch.nn as nn

    dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    log("device", dev, torch.cuda.get_device_name(0) if dev.type == "cuda" else "")

    class ABMIL(nn.Module):  # gated attention MIL (Ilse et al. 2018), as in train_mil.py
        def __init__(self, d, h=256, att=128, p=0.25):
            super().__init__()
            self.fc = nn.Sequential(nn.Linear(d, h), nn.ReLU(), nn.Dropout(p))
            self.a = nn.Sequential(nn.Linear(h, att), nn.Tanh())
            self.b = nn.Sequential(nn.Linear(h, att), nn.Sigmoid())
            self.c = nn.Linear(att, 1)
            self.cls = nn.Linear(h, 1)

        def forward(self, x):
            hdn = self.fc(x)
            w = torch.softmax(self.c(self.a(hdn) * self.b(hdn)).squeeze(-1), dim=0)
            return self.cls((w[:, None] * hdn).sum(0)).squeeze()

    # bags in RAM, keyed by patient (1c/1b: DX or primary TS) or by slide barcode (1e: TS and BS of the pairs)
    BAGS = {}
    t = time.time()
    if a.slide == "DX":
        arr = np.load(os.path.join(a.feats, f"{a.name}_cap{a.cap}.npy"), mmap_mode="r")
        for p in both:
            BAGS[p] = np.array(arr[idx.at[p, "row"], :idx.at[p, "n_kept"]])
        del arr
    else:  # one copy per slide: a patient's primary TS is often also the TS of its pair
        by_slide = {}

        def load(s):
            if s not in by_slide:
                by_slide[s] = np.load(os.path.join(a.frozen, f"{s}.npz"))["bag"]
            return by_slide[s]

        for p in (prim.index if a.part == "1e" else both):
            BAGS[p] = load(prim.at[p, "slide_barcode"])
        if a.part == "1e":
            for s in man[man.pair_id != ""].slide_barcode:
                BAGS[s] = load(s)
    D = next(iter(BAGS.values())).shape[1]
    log(f"{len(BAGS)} bags ({a.slide}, part {a.part}) loaded in {time.time() - t:.0f}s, "
        f"{sum(b.nbytes for b in {id(b): b for b in BAGS.values()}.values()) / 1e9:.1f} GB")

    def bag(k, train, rng):
        x = BAGS[k]
        if train and len(x) > a.train_patches:
            x = x[rng.choice(len(x), a.train_patches, replace=False)]
        return torch.from_numpy(np.ascontiguousarray(x)).to(dev, non_blocking=True)

    def fit_predict(tr, y_tr, te_sets, seed):
        """train on bags tr (keys) with labels y_tr; return one probability array per list of keys in te_sets"""
        torch.manual_seed(seed)
        rng = np.random.default_rng(seed)
        net = ABMIL(D).to(dev)
        opt = torch.optim.Adam(net.parameters(), lr=a.lr, weight_decay=1e-4)
        pos = max(y_tr.sum(), 1)
        lossf = nn.BCEWithLogitsLoss(pos_weight=torch.tensor((len(y_tr) - pos) / pos, device=dev))
        scaler = torch.cuda.amp.GradScaler(enabled=dev.type == "cuda")
        yt = torch.tensor(y_tr, dtype=torch.float32, device=dev)
        for _ in range(a.epochs):
            net.train()
            for i in rng.permutation(len(tr)):
                with torch.autocast(device_type=dev.type, dtype=torch.float16, enabled=dev.type == "cuda"):
                    loss = lossf(net(bag(tr[i], True, rng).float()), yt[i])
                opt.zero_grad(set_to_none=True)
                scaler.scale(loss).backward()
                scaler.step(opt)
                scaler.update()
        net.eval()
        with torch.no_grad(), torch.autocast(device_type=dev.type, dtype=torch.float16, enabled=dev.type == "cuda"):
            return [np.array([torch.sigmoid(net(bag(k, False, rng).float())).item() for k in te]) for te in te_sets]

    preds = []
    if a.part in COMPARISONS:
        for marker, ihc, mol in COMPARISONS[a.part]:
            d = L.loc[L.index.intersection(both), [ihc, mol, "site"]].dropna()
            pts = d.index.to_numpy()
            for lab in (ihc, mol):
                y = d[lab].astype(int).values
                for s in range(a.seeds):
                    pr = np.zeros(len(pts))
                    for tr, te in site_folds(pts, s):
                        pr[te] = fit_predict(list(pts[tr]), y[tr], [list(pts[te])], s)[0]
                    preds.append(pd.DataFrame({"exp": "1c", "marker": marker, "comparison": f"{mol} vs {ihc}",
                                               "slide": a.slide, "label": lab, "seed": s, "patient": pts, "y": y,
                                               "p": pr}))
                    log(f"1c {a.slide} {mol} vs {ihc} label {lab} seed {s}: AUC {roc_auc_score(y, pr):.3f} (n={len(pts)})")
        for ihc, mol in PAIRS_1B[a.part]:
            d = L.loc[L.index.intersection(both), [ihc, mol]].dropna()
            conc, disc = d[d[ihc] == d[mol]], d[d[ihc] != d[mol]]
            for s in range(a.seeds):
                p = fit_predict(list(conc.index), conc[ihc].astype(int).values, [list(disc.index)], s)[0]
                preds.append(pd.DataFrame({"exp": "1b", "marker": ihc.split("_")[0], "comparison": f"1b {ihc} vs {mol}",
                                           "slide": a.slide, "label": f"1b {ihc} vs {mol}", "seed": s,
                                           "patient": disc.index, "y": disc[mol].values, "p": p}))
                k = int(((p >= 0.5).astype(int) == disc[mol].astype(int).values).sum())
                log(f"1b {a.slide} {ihc} vs {mol} seed {s}: molecular {k} / IHC {len(disc) - k} (n={len(disc)})")
    else:  # 1e, as frozen_1e.py
        pairs = man[man.pair_id != ""]
        ts = pairs[pairs.section_type == "TS"].drop_duplicates("pair_id").set_index("pair_id")
        bs = pairs[pairs.section_type == "BS"].drop_duplicates("pair_id").set_index("pair_id")
        pid = ts.index.intersection(bs.index)
        ts, bs = ts.loc[pid], bs.loc[pid]
        log(f"1e: {len(pid)} complete TS+BS pairs, {len(prim)} primary TS slides")
        for lab in ["ER_IHC", "ESR1_RNApos"]:
            tr_pts = prim.index[prim.index.isin(L.index[L[lab].notna() & L.site.notna()])].to_numpy()
            y_all = L.loc[tr_pts, lab].astype(int).values
            for s in range(a.seeds):
                for tr, te in site_folds(tr_pts, s):
                    te_pairs = pid[ts.patient.isin(tr_pts[te]).values]
                    if not len(te_pairs):
                        continue
                    pt, pb = fit_predict(list(tr_pts[tr]), y_all[tr], [list(ts.loc[te_pairs, "slide_barcode"]),
                                                                        list(bs.loc[te_pairs, "slide_barcode"])], s)
                    for sec, pp in (("TS", pt), ("BS", pb)):
                        preds.append(pd.DataFrame({"exp": "1e", "marker": "ER", "label": lab, "seed": s,
                                                   "pair_id": te_pairs, "patient": ts.loc[te_pairs, "patient"].values,
                                                   "section": sec, "y": L.loc[ts.loc[te_pairs, "patient"], lab].values,
                                                   "p": pp}))
                log(f"1e {lab} seed {s} done")
    os.makedirs(a.out, exist_ok=True)
    f = os.path.join(a.out, f"frozen_abmil_raw_{a.slide}_{a.part}.csv.gz")
    pd.concat(preds).to_csv(f, index=False)
    log("saved", f)

else:
    R = lambda f: os.path.join(a.out, f)  # noqa: E731
    raw = pd.concat([pd.read_csv(R(f)) for f in sorted(os.listdir(a.out)) if f.startswith("frozen_abmil_raw_")])
    rng = np.random.default_rng(0)
    rows = []
    # ---- 1c and 1b: seed-averaged predictions in the frozen_1c_oof.csv.gz format ----
    c = raw[raw.exp.isin(["1c", "1b"])]
    oof = c.groupby(["marker", "comparison", "slide", "label", "patient"], sort=False).agg(y=("y", "first"),
                                                                                         p=("p", "mean")).reset_index()
    oof.to_csv(R("frozen_abmil_1c_oof.csv.gz"), index=False)
    for (marker, comp), g in oof[~oof.comparison.str.startswith("1b")].groupby(["marker", "comparison"], sort=False):
        mol, ihc = comp.split(" vs ")
        w = g.pivot_table(index="patient", columns=["slide", "label"], values="p")
        if w.shape[1] < 4:
            log(f"1c {comp}: missing slide/label cells ({list(w.columns)}), skipped")
            continue
        yy = g.pivot_table(index="patient", columns="label", values="y").loc[w.index]
        Y = {lab: yy[lab].astype(int).values for lab in (ihc, mol)}
        P = {(s, lab): w[s, lab].values for s in ("DX", "TS") for lab in (ihc, mol)}
        n = len(w)

        def aucs(i):
            return {k: roc_auc_score(Y[k[1]][i], P[k][i]) for k in P}

        def inter(A):
            return (A["TS", mol] - A["TS", ihc]) - (A["DX", mol] - A["DX", ihc])

        A = aucs(np.arange(n))
        bs = []
        for _ in range(a.boot):
            i = rng.integers(0, n, n)
            if min(len(set(Y[ihc][i])), len(set(Y[mol][i]))) == 2:
                B = aucs(i)
                bs.append([inter(B), B["DX", mol] - B["DX", ihc], B["TS", mol] - B["TS", ihc]])
        bs = np.array(bs)
        lo, hi = np.percentile(bs, [2.5, 97.5], axis=0)
        # interaction per seed (spread over training runs, no bootstrap)
        per_seed = []
        for s, gs in c[(c.comparison == comp)].groupby("seed"):
            ws = gs.pivot_table(index="patient", columns=["slide", "label"], values="p").loc[w.index]
            As = {k: roc_auc_score(Y[k[1]], ws[k].values) for k in P}
            per_seed.append(inter(As))
        row = {"marker": marker, "comparison": comp, "model": "UNI ABMIL", "n": n, "pos_ihc": int(Y[ihc].sum()),
               "pos_mol": int(Y[mol].sum()), "AUC_DX_ihc": A["DX", ihc], "AUC_DX_mol": A["DX", mol],
               "AUC_TS_ihc": A["TS", ihc], "AUC_TS_mol": A["TS", mol],
               "delta_DX": A["DX", mol] - A["DX", ihc], "delta_DX_lo": lo[1], "delta_DX_hi": hi[1],
               "delta_TS": A["TS", mol] - A["TS", ihc], "delta_TS_lo": lo[2], "delta_TS_hi": hi[2],
               "interaction": inter(A), "interaction_lo": lo[0], "interaction_hi": hi[0],
               "interaction_p_boot": 2 * min((bs[:, 0] <= 0).mean(), (bs[:, 0] >= 0).mean()),
               "interaction_seed_mean": np.mean(per_seed), "interaction_seed_sd": np.std(per_seed)}
        rows.append(row)
        log(f"1c ABMIL {comp}: n={n} | DX {A['DX', ihc]:.3f}->{A['DX', mol]:.3f} (d {row['delta_DX']:+.3f} "
            f"[{lo[1]:+.3f}, {hi[1]:+.3f}]) | TS {A['TS', ihc]:.3f}->{A['TS', mol]:.3f} (d {row['delta_TS']:+.3f} "
            f"[{lo[2]:+.3f}, {hi[2]:+.3f}]) | interaction {row['interaction']:+.3f} [{lo[0]:+.3f}, {hi[0]:+.3f}] "
            f"p={row['interaction_p_boot']:.3f} (per seed {np.mean(per_seed):+.3f} ± {np.std(per_seed):.3f})")
    for (marker, comp), g in oof[oof.comparison.str.startswith("1b")].groupby(["marker", "comparison"], sort=False):
        out = {"marker": marker, "comparison": comp, "model": "UNI ABMIL", "n": g.patient.nunique()}
        for s, gs in g.groupby("slide"):
            k = int(((gs.p >= 0.5).astype(int) == gs.y.astype(int)).sum())
            out.update({f"side_mol_{s}": k, f"side_ihc_{s}": len(gs) - k, f"p_{s}": binomtest(k, len(gs)).pvalue})
        rows.append(out)
        log(f"{comp} ABMIL: n={out['n']} | " + " | ".join(
            f"{s} molecular {out.get(f'side_mol_{s}')}:{out.get(f'side_ihc_{s}')} (p={out.get(f'p_{s}', np.nan):.3g})"
            for s in ("DX", "TS") if f"p_{s}" in out))
    pd.DataFrame(rows).to_csv(R("frozen_abmil_1c.csv"), index=False)

    # ---- 1e: seed-averaged P per pair in the frozen_1e_pairs.csv format, agreement as in frozen_1e.py ----
    e = raw[raw.exp == "1e"]
    if len(e):
        def icc21(x, y):
            M = np.column_stack([x, y])
            n, k = M.shape
            gm = M.mean()
            msr = k * ((M.mean(1) - gm) ** 2).sum() / (n - 1)
            msc = n * ((M.mean(0) - gm) ** 2).sum() / (k - 1)
            sse = ((M - M.mean(1, keepdims=True) - M.mean(0, keepdims=True) + gm) ** 2).sum()
            mse = sse / ((n - 1) * (k - 1))
            return (msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n)

        def logit(p):
            p = np.clip(p, 1e-4, 1 - 1e-4)
            return np.log(p / (1 - p))

        def stats(pt, pb, y):
            ct, cb = (pt >= 0.5).astype(int), (pb >= 0.5).astype(int)
            out = {"ICC_logit": icc21(logit(pt), logit(pb)), "pearson_P": np.corrcoef(pt, pb)[0, 1],
                   "kappa": cohen_kappa_score(ct, cb), "pct_flipped": 100 * (ct != cb).mean(),
                   "median_abs_dP": np.median(np.abs(pt - pb))}
            if len(set(y)) == 2:
                out.update(AUC_TS=roc_auc_score(y, pt), AUC_BS=roc_auc_score(y, pb))
            return out

        tum = man.set_index("slide_barcode").percent_tumor_cells
        pairs = man[man.pair_id != ""]
        tsb = pairs[pairs.section_type == "TS"].drop_duplicates("pair_id").set_index("pair_id").slide_barcode
        bsb = pairs[pairs.section_type == "BS"].drop_duplicates("pair_id").set_index("pair_id").slide_barcode
        per_pair, rows = [], []
        for lab, g in e.groupby("label", sort=False):
            w = g.pivot_table(index="pair_id", columns="section", values="p")  # mean over seeds
            meta = g.groupby("pair_id")[["patient", "y"]].first().loc[w.index]
            pt, pb, y = w.TS.values, w.BS.values, meta.y.astype(int).values
            est = stats(pt, pb, y)
            bsx = pd.DataFrame([stats(pt[i], pb[i], y[i]) for i in (rng.integers(0, len(w), len(w)) for _ in range(a.boot))])
            row = {"label": lab, "model": "UNI ABMIL", "n_pairs": len(w), "pos": int(y.sum())}
            for k, v in est.items():
                row[k] = v
                row[k + "_lo"], row[k + "_hi"] = np.nanpercentile(bsx[k], [2.5, 97.5])
            rows.append(row)
            per_pair.append(pd.DataFrame({"label": lab, "pair_id": w.index, "patient": meta.patient.values, "y": y,
                                          "P_TS": pt, "P_BS": pb, "tumour_TS": tum.reindex(tsb.loc[w.index]).values,
                                          "tumour_BS": tum.reindex(bsb.loc[w.index]).values}))
            log(f"1e ABMIL {lab}: {len(w)} pairs | ICC(logit P) {row['ICC_logit']:.3f} [{row['ICC_logit_lo']:.3f}, "
                f"{row['ICC_logit_hi']:.3f}] | r(P) {row['pearson_P']:.3f} | kappa {row['kappa']:.3f} "
                f"[{row['kappa_lo']:.3f}, {row['kappa_hi']:.3f}] | flipped {row['pct_flipped']:.1f}% "
                f"[{row['pct_flipped_lo']:.1f}, {row['pct_flipped_hi']:.1f}] | AUC TS {row.get('AUC_TS', np.nan):.3f} "
                f"BS {row.get('AUC_BS', np.nan):.3f}")
        pd.DataFrame(rows).to_csv(R("frozen_abmil_1e.csv"), index=False)
        pd.concat(per_pair).to_csv(R("frozen_abmil_1e_pairs.csv"), index=False)
    log("saved", R("frozen_abmil_1c.csv"), R("frozen_abmil_1c_oof.csv.gz"), R("frozen_abmil_1e.csv"),
        R("frozen_abmil_1e_pairs.csv"))
