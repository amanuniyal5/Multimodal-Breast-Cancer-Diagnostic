"""X2 comparison tables from the bench_mil.py outputs (docs/DESIGN_SPEC.md §5.2).

Reads <dir>/*_<name>_oof.csv.gz, averages predictions over seeds per patient and writes
  aucs.csv         method x task: AUC (DeLong CI), AUPRC, Brier, ECE; difference vs ABMIL and vs the mean-pool linear
                   probe on the same patients (DeLong); HER2 1+ vs 0 also with the ERBB2 mRNA comparator
  preference.csv   method x label pair: ΔAUC molecular − IHC (generalised DeLong) and 1b counts (binomial)
  runtime.csv      seconds per fit and parameters (GPU methods)
  repro.txt        ABMIL through this runner vs the fixed-split train_mil.py ABMIL (ER group, per seed)
Usage: python scripts/analyse_x2.py [--dir results/bench/x2] [--name uni_v1_20x]
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import binomtest, pearsonr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bench import stats, tasks  # noqa: E402
from bench.metrics import classification  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--dir", default="results/bench/x2")
ap.add_argument("--name", default="uni_v1_20x")
ap.add_argument("--ref-abmil", default="results/mil_uni_v1_20x_fixedsplits_abmil_oof.csv.gz")
a = ap.parse_args()
pd.set_option("display.width", 250)
O = pd.concat([pd.read_csv(f) for f in glob.glob(os.path.join(a.dir, f"*_{a.name}_oof.csv.gz"))])
ORDER = [m for m in ["mean", "meanstd", "meanmax", "knn25", "proto", "abmil", "clam", "transmil", "dsmil", "acmil"]
         if m in set(O.method)]
A1 = O[O.exp == "1a"].groupby(["method", "task", "group", "patient"], sort=False).agg(y=("y", "first"), p=("p", "mean"),
                                                                                        n_seeds=("seed", "nunique"))
A1 = A1.reset_index()
L = tasks.labels()


def block(m, t, g):
    return A1[(A1.method == m) & (A1.task == t) & (A1.group == g)].set_index("patient")


rows = []
for (t, g), _ in A1.groupby(["task", "group"], sort=False):
    for m in ORDER:
        b = block(m, t, g)
        if b.empty:
            continue
        y, p = b.y.astype(int).to_numpy(), b.p.to_numpy()
        auc, lo, hi = stats.auc_ci(y, p)
        r = {"method": m, "task": t, "group": g, "n": len(y), "n_seeds": int(b.n_seeds.min()), "auc": auc, "lo": lo,
             "hi": hi, **{k: v for k, v in classification(y, p).items() if k != "auc"}}
        for tag, ref in (("abmil", "abmil"), ("mean", "mean")):
            rb = block(ref, t, g)
            if m != ref and not rb.empty and set(rb.index) == set(b.index):
                d = stats.delong(y, rb.loc[b.index, "p"].to_numpy(), y, p)
                r.update({f"diff_vs_{tag}": d["diff"], f"lo_vs_{tag}": d["lo"], f"hi_vs_{tag}": d["hi"],
                          f"p_vs_{tag}": d["p"]})
        if t == "HER2low_1plus":
            r["auc_ERBB2_mRNA"] = stats.safe_auc(y, L.loc[b.index, "ERBB2_log2"].to_numpy())
        rows.append(r)
A = pd.DataFrame(rows)
A.to_csv(os.path.join(a.dir, "aucs.csv"), index=False)
print("AUC (seed-averaged):\n" + A.pivot(index="method", columns="task", values="auc").reindex(ORDER).round(3)
      .to_string())
print("\nvs ABMIL / vs mean pool (ΔAUC, p):\n" + A[["method", "task", "auc", "diff_vs_abmil", "p_vs_abmil",
                                                      "diff_vs_mean", "p_vs_mean"]].round(3).to_string(index=False))
print("\ncalibration:\n" + A[A.task.isin(["ER_IHC", "HER2_STRICT"])][["method", "task", "brier", "ece"]].round(3)
      .to_string(index=False))

pref = []
for m in ORDER:
    for g in A1[A1.method == m].group.unique():
        cols = g.split("|")
        for mol in cols[1:]:
            bi, bm = block(m, cols[0], g), block(m, mol, g)
            d = stats.delong(bi.y.astype(int), bi.p, bm.loc[bi.index, "y"].astype(int), bm.loc[bi.index, "p"])
            r = {"method": m, "pair": f"{cols[0]} vs {mol}", "n": len(bi), "auc_ihc": d["auc_ref"],
                 "auc_mol": d["auc_alt"], "delta": d["diff"], "lo": d["lo"], "hi": d["hi"], "p": d["p"]}
            b1 = O[(O.method == m) & (O.exp == "1b") & (O.task == r["pair"])]
            if len(b1):
                k = int(((b1.p >= 0.5).astype(int) == b1.y.astype(int)).sum())
                r.update(n_disc=len(b1), side_molecular=k, side_IHC=len(b1) - k, p_binom=binomtest(k, len(b1)).pvalue)
            pref.append(r)
P = pd.DataFrame(pref)
P.to_csv(os.path.join(a.dir, "preference.csv"), index=False)
print("\nlabel preference:\n" + P.round(3).to_string(index=False))

rt = [pd.read_csv(f) for f in glob.glob(os.path.join(a.dir, f"*_{a.name}_runtime.csv"))]
if rt:
    R = pd.concat(rt).groupby("method").agg(fits=("seconds", "size"), sec_per_fit=("seconds", "mean"),
                                           total_h=("seconds", lambda s: s.sum() / 3600), params=("params", "first"))
    R.to_csv(os.path.join(a.dir, "runtime.csv"))
    print("\nruntime:\n" + R.round(1).to_string())

if "abmil" in ORDER and os.path.exists(a.ref_abmil):
    ref = pd.read_csv(a.ref_abmil)
    lines = []
    for col in ["ER_IHC", "ESR1_RNApos", "ER_RPPApos"]:
        for s in sorted(O[(O.method == "abmil") & (O.task == col)].seed.unique()):
            x = O[(O.method == "abmil") & (O.exp == "1a") & (O.task == col) & (O.seed == s)].set_index("patient")
            r0 = ref[(ref.model == "abmil") & (ref.exp == "1a") & (ref.label == col) & (ref.seed == s)].set_index("patient")
            same = set(x.index) == set(r0.index)
            r0 = r0.loc[x.index] if same else r0
            a1, a0 = stats.safe_auc(x.y.astype(int), x.p), stats.safe_auc(r0.y.astype(int), r0.p)
            lines.append(f"{col} seed {s}: same patients {same} (n {len(x)} vs {len(r0)}); AUC here {a1:.4f} vs "
                         f"train_mil {a0:.4f} (diff {a1 - a0:+.4f}); r(p) "
                         f"{pearsonr(x.p, r0.p)[0] if same else float('nan'):.3f}")
    txt = "\n".join(lines)
    open(os.path.join(a.dir, "repro.txt"), "w").write(txt + "\n")
    print("\nABMIL reproduction:\n" + txt)
