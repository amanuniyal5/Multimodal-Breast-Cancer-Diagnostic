"""Tables for the GPU runs of X1 (ABMIL on a second encoder), X3 (ABMIL regression), X4b (multi-task), X4c (DANN) and
X4d (noisy labels, ABMIL). Predictions are averaged over seeds per patient; comparisons on the same patients (DeLong,
generalised to two labels where needed; paired bootstrap for correlations).

  x1_abmil/aucs.csv       H-optimus-0 ABMIL vs H-optimus-0 mean pool (X1) and vs UNI ABMIL (X2), per task; label preference
  x3/abmil_regression.csv ABMIL regression vs ridge (UNI): Pearson r; regression score as classifier vs ABMIL classifier
  x4/multitask.csv        multi-task vs single-task ABMIL (same patients), PAM50 vs the multinomial linear probe
  x4/dann.csv             DANN vs plain ABMIL: AUC and hospital-ID accuracy of the pooled embeddings (10 largest sites)
  x4/noisy_abmil.csv      noisy-label ABMIL variants vs plain BCE, scored against IHC, RNA and RPPA
Usage: python scripts/analyse_ext.py
"""
import os
import sys
import warnings

os.environ.setdefault("LOKY_MAX_CPU_COUNT", "4")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import binomtest  # noqa: E402
from sklearn.metrics import balanced_accuracy_score, roc_auc_score  # noqa: E402
from sklearn.model_selection import StratifiedKFold, cross_val_predict  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bench import stats, tasks  # noqa: E402
from common import linear_probe  # noqa: E402

warnings.filterwarnings("ignore")
pd.set_option("display.width", 250)
L = tasks.labels()
B = "results/bench"


def avg(df, keys=("variant", "task")):
    return df.groupby(list(keys) + ["patient"]).agg(y=("y", "first"), p=("p", "mean")).reset_index()


def show(name, D):
    D.to_csv(name, index=False)
    print(f"\n{name}:\n" + D.round(3).to_string(index=False))


def r_ci(y, p, boot=1000, seed=0):
    r = np.random.default_rng(seed)
    v = [np.corrcoef(y[i], p[i])[0, 1] for i in (r.integers(0, len(y), len(y)) for _ in range(boot))]
    return np.corrcoef(y, p)[0, 1], np.percentile(v, 2.5), np.percentile(v, 97.5)


X2A = None
if os.path.exists(f"{B}/x2/abmil_uni_v1_20x_oof.csv.gz"):
    X2A = pd.read_csv(f"{B}/x2/abmil_uni_v1_20x_oof.csv.gz")
    X2A = X2A[X2A.exp == "1a"].groupby(["task", "patient"]).agg(y=("y", "first"), p=("p", "mean"))

# ---------------- X1: ABMIL on H-optimus-0 ----------------
f = f"{B}/x1_abmil/abmil_hoptimus0_20x_oof.csv.gz"
if os.path.exists(f):
    O = pd.read_csv(f)
    A = O[O.exp == "1a"].groupby(["task", "patient"]).agg(y=("y", "first"), p=("p", "mean"))
    X1 = pd.read_csv(f"{B}/x1/oof.csv.gz")
    rows = []
    for t in A.index.get_level_values(0).unique():
        a_ = A.loc[t]
        y = a_.y.astype(int).to_numpy()
        auc, lo, hi = stats.auc_ci(y, a_.p.to_numpy())
        r = {"task": t, "n": len(y), "auc_hopt_abmil": auc, "lo": lo, "hi": hi}
        m = X1[(X1.encoder == "hoptimus0") & (X1.task == t)].set_index("patient")
        ids = a_.index.intersection(m.index)
        if len(ids) == len(a_):
            d = stats.delong(y, m.loc[a_.index, "p"].to_numpy(), y, a_.p.to_numpy())
            r.update(auc_hopt_mean=d["auc_ref"], diff_vs_mean=d["diff"], lo_vs_mean=d["lo"], hi_vs_mean=d["hi"],
                     p_vs_mean=d["p"])
        if X2A is not None and t in X2A.index.get_level_values(0):
            u = X2A.loc[t]
            ids = a_.index.intersection(u.index)
            d = stats.delong(a_.loc[ids, "y"].astype(int).to_numpy(), u.loc[ids, "p"].to_numpy(),
                             a_.loc[ids, "y"].astype(int).to_numpy(), a_.loc[ids, "p"].to_numpy())
            r.update(n_common_uni=len(ids), auc_uni_abmil=d["auc_ref"], diff_vs_uni_abmil=d["diff"],
                     lo_vs_uni=d["lo"], hi_vs_uni=d["hi"], p_vs_uni=d["p"])
        rows.append(r)
    for ihc, mol in (("ER_IHC", "ESR1_RNApos"), ("ER_IHC", "ER_RPPApos"), ("PR_IHC", "PGR_prevmatch"),
                     ("HER2_STRICT", "ERBB2_AMP")):
        bi, bm = A.loc[ihc], A.loc[mol].loc[A.loc[ihc].index]
        d = stats.delong(bi.y.astype(int), bi.p, bm.y.astype(int), bm.p)
        r = {"task": f"pref {mol} - {ihc}", "n": len(bi), "auc_hopt_abmil": d["diff"], "lo": d["lo"], "hi": d["hi"],
             "p_vs_mean": d["p"]}
        b1 = O[(O.exp == "1b") & (O.task == f"{ihc} vs {mol}")]
        if len(b1):
            k = int(((b1.p >= 0.5).astype(int) == b1.y.astype(int)).sum())
            r.update(n_common_uni=len(b1), auc_uni_abmil=k, diff_vs_uni_abmil=binomtest(k, len(b1)).pvalue)
        rows.append(r)
    show(f"{B}/x1_abmil/aucs.csv", pd.DataFrame(rows))
    print("  (preference rows: auc_hopt_abmil = delta, lo/hi = CI, p_vs_mean = p; 1b: n_common_uni = n discordant, "
          "auc_uni_abmil = molecular side, diff_vs_uni_abmil = binomial p)")

# ---------------- X3: ABMIL regression ----------------
f = f"{B}/x3/regress_uni_v1_20x_oof.csv.gz"
if os.path.exists(f):
    O = avg(pd.read_csv(f))
    R = pd.read_csv(f"{B}/x3/oof_reg.csv.gz")
    R = R[R.encoder == "uni_v1"].set_index(["gene", "patient"])
    rows = []
    for g in O.task.unique():
        a_ = O[O.task == g].set_index("patient")
        ids = a_.index.intersection(R.loc[g].index)
        y, pa, pr = a_.loc[ids, "y"].to_numpy(), a_.loc[ids, "p"].to_numpy(), R.loc[g].loc[ids, "pred"].to_numpy()
        ra, la, ha = r_ci(y, pa)
        rr, lr_, hr = r_ci(y, pr)
        rb = np.random.default_rng(1)
        dd = []
        for _ in range(1000):
            i = rb.integers(0, len(y), len(y))
            dd.append(np.corrcoef(y[i], pa[i])[0, 1] - np.corrcoef(y[i], pr[i])[0, 1])
        r = {"gene": g, "n": len(ids), "r_abmil": ra, "lo": la, "hi": ha, "r_ridge": rr, "ridge_lo": lr_,
             "ridge_hi": hr, "diff_r": ra - rr, "diff_lo": np.percentile(dd, 2.5), "diff_hi": np.percentile(dd, 97.5)}
        rows.append(r)
        pairs = {"ESR1_log2": ("ER_IHC", "ESR1_RNApos"), "PGR_log2": ("PR_IHC", "PGR_prevmatch"),
                 "ERBB2_log2": ("HER2_STRICT", "ERBB2_AMP")}[g]
        for lab in pairs:
            if X2A is None or lab not in X2A.index.get_level_values(0):
                continue
            c = X2A.loc[lab]
            ii = c.index.intersection(a_.index)
            yl = c.loc[ii, "y"].astype(int).to_numpy()
            d = stats.delong(yl, c.loc[ii, "p"].to_numpy(), yl, a_.loc[ii, "p"].to_numpy())
            rows.append({"gene": f"{g} as score for {lab}", "n": len(ii), "r_abmil": d["auc_alt"],
                         "r_ridge": d["auc_ref"], "diff_r": d["diff"], "diff_lo": d["lo"], "diff_hi": d["hi"],
                         "lo": d["p"]})
    show(f"{B}/x3/abmil_regression.csv", pd.DataFrame(rows))
    print("  (score rows: r_abmil = AUC of the ABMIL regression score, r_ridge = AUC of the ABMIL classifier (X2), "
          "lo = DeLong p)")

# ---------------- X4b: multi-task ----------------
fm, fs = f"{B}/x4/multitask_uni_v1_20x_oof.csv.gz", f"{B}/x4/single_uni_v1_20x_oof.csv.gz"
if os.path.exists(fm):
    M = pd.read_csv(fm)
    rows = []
    MB = avg(M[M.exp == "1a"])
    if os.path.exists(fs):
        S = avg(pd.read_csv(fs))
        for t in ["ER_IHC", "PR_IHC", "HER2_STRICT"]:
            a_, b_ = MB[MB.task == t].set_index("patient"), S[S.task == t].set_index("patient")
            ids = a_.index.intersection(b_.index)
            y = a_.loc[ids, "y"].astype(int).to_numpy()
            d = stats.delong(y, b_.loc[ids, "p"].to_numpy(), y, a_.loc[ids, "p"].to_numpy())
            rows.append({"task": t, "n": len(ids), "auc_single": d["auc_ref"], "auc_multitask": d["auc_alt"],
                         "diff": d["diff"], "lo": d["lo"], "hi": d["hi"], "p": d["p"]})
    MR = avg(M[M.exp == "reg"])
    for g in MR.task.unique():
        a_ = MR[MR.task == g]
        r, lo, hi = r_ci(a_.y.to_numpy(), a_.p.to_numpy())
        rows.append({"task": f"{g} (Pearson r)", "n": len(a_), "auc_multitask": r, "lo": lo, "hi": hi})
    P = M[M.exp == "pam50"]
    pc = [c for c in P.columns if c.startswith("p_") and c != "p"]
    PP = P.groupby("patient").agg({**{c: "mean" for c in pc}, "y": "first"})
    y = PP.y.astype(int).to_numpy()
    pr = PP[pc].to_numpy()
    rows.append({"task": "PAM50 macro AUC (one-vs-rest)", "n": len(y),
                 "auc_multitask": roc_auc_score(y, pr, multi_class="ovr", average="macro"),
                 "diff": balanced_accuracy_score(y, pr.argmax(1))})
    if os.path.exists(f"{B}/x4/pam50_lr.csv"):
        lr = pd.read_csv(f"{B}/x4/pam50_lr.csv").set_index("encoder").loc["uni_v1"]
        rows[-1].update(auc_single=lr.macro_auc_ovr, lo=lr.balanced_accuracy)
    show(f"{B}/x4/multitask.csv", pd.DataFrame(rows))
    print("  (PAM50 row: auc_single = UNI multinomial linear probe macro AUC; diff / lo = balanced accuracy multitask / probe)")

# ---------------- X4c: DANN ----------------
f = f"{B}/x4/dann_uni_v1_20x_oof.csv.gz"
if os.path.exists(f):
    D = avg(pd.read_csv(f))
    rows = []
    for t in D.task.unique():
        a_, b_ = D[(D.task == t) & (D.variant == "plain")].set_index("patient"), \
                 D[(D.task == t) & (D.variant == "dann")].set_index("patient")
        y = a_.y.astype(int).to_numpy()
        d = stats.delong(y, a_.p.to_numpy(), y, b_.loc[a_.index, "p"].to_numpy())
        rows.append({"task": t, "n": len(y), "auc_plain": d["auc_ref"], "auc_dann": d["auc_alt"], "diff": d["diff"],
                     "lo": d["lo"], "hi": d["hi"], "p": d["p"]})
    a_, b_ = D[(D.task == "ER_IHC") & (D.variant == "dann")].set_index("patient"), \
             D[(D.task == "ESR1_RNApos") & (D.variant == "dann")].set_index("patient")
    d = stats.delong(a_.y.astype(int), a_.p, b_.loc[a_.index, "y"].astype(int), b_.loc[a_.index, "p"])
    rows.append({"task": "pref ESR1_RNApos - ER_IHC (DANN)", "n": len(a_), "diff": d["diff"], "lo": d["lo"],
                 "hi": d["hi"], "p": d["p"]})
    E = np.load(f"{B}/x4/dann_uni_v1_20x_embeddings.npz", allow_pickle=True)
    for t in ["ER_IHC", "HER2_STRICT"]:
        pts = E[f"patients|{t}"].astype(str)
        site = L.loc[pts, "site"].to_numpy(dtype=str)
        top = pd.Series(site).value_counts().index[:10]
        m = np.isin(site, top)
        for v in ("plain", "dann"):
            pred = cross_val_predict(linear_probe(), E[f"{v}|{t}"][m], site[m],
                                     cv=StratifiedKFold(5, shuffle=True, random_state=0))
            acc = (pred == site[m]).mean()
            lo, hi = stats.wilson(int(round(acc * m.sum())), int(m.sum()))
            rows.append({"task": f"hospital ID from {v} ABMIL embedding ({t})", "n": int(m.sum()), "auc_plain": acc,
                         "lo": lo, "hi": hi, "p": pd.Series(site[m]).value_counts(normalize=True).max()})
    show(f"{B}/x4/dann.csv", pd.DataFrame(rows))
    print("  (hospital-ID rows: auc_plain = accuracy, p = chance)")

# ---------------- X4d: noisy labels (ABMIL) ----------------
f = f"{B}/x4/noisy_uni_v1_20x_oof.csv.gz"
if os.path.exists(f):
    N = pd.read_csv(f)
    N = N.groupby(["variant", "patient"]).agg(y=("y", "first"), y_rna=("y_rna", "first"), y_rppa=("y_rppa", "first"),
                                              p=("p", "mean")).reset_index()
    base = N[N.variant == "bce"].set_index("patient")
    rows = []
    for v in N.variant.unique():
        a_ = N[N.variant == v].set_index("patient").loc[base.index]
        for lab, col in (("ER_IHC", "y"), ("ESR1_RNApos", "y_rna"), ("ER_RPPApos", "y_rppa")):
            y = base[col].astype(int).to_numpy()
            d = stats.delong(y, base.p.to_numpy(), y, a_.p.to_numpy())
            rows.append({"variant": v, "scored_vs": lab, "n": len(y), "auc": d["auc_alt"], "diff_vs_bce": d["diff"],
                         "lo": d["lo"], "hi": d["hi"], "p": d["p"]})
    show(f"{B}/x4/noisy_abmil.csv", pd.DataFrame(rows))
print("ANALYSE_EXT_OK")
