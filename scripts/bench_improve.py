"""CPU parts of X3 (PAM50 gene panel) and X4 (c, d, e) on patient means (docs/DESIGN_SPEC.md §5.3-5.4).
Fixed site-grouped folds; inner folds from folds.csv (bench.folds.inner); 5 seeds averaged; mean-pool linear probe.

panel  X3 secondary: ridge (alpha per gene by inner CV, all 50 genes in one multi-output fit) from UNI, H-optimus-0,
       TITAN, CONCH means to log2 expression of the 50 PAM50 genes; out-of-fold Pearson r with patient-bootstrap CI,
       ranked; plus a multinomial PAM50 classifier on the means (reference for multi-task ABMIL X4b).
                                                                  -> x3/panel.csv, x3/oof_panel.csv.gz, x4/pam50_lr.csv
site   X4c(i): per-hospital feature standardisation (each hospital's patients z-scored with that hospital's own mean /
       SD of the unlabelled features; hospitals with < 5 patients use the cohort statistics). AUC per task vs plain,
       DeLong; RNA - IHC preference after the change; hospital-ID accuracy (10 largest sites, as X1) before / after.
                                                                  -> x4/site_std.csv, x4/site_std_hospital_id.csv
mag    X4e: UNI 20x + UNI 10x means concatenated (and + patch SD) vs UNI 20x.          -> x4/multimag.csv
noisy  X4d (linear heads): ER, train on IHC (plain), consensus (IHC = RNA only), soft (latent-class posterior from IHC /
       RNA / RPPA fitted on the training fold, as duplicated weighted rows) and label smoothing 0.1; scored against IHC,
       RNA and RPPA; UNI, H-optimus-0, TITAN.                                           -> x4/noisy_linear.csv
Usage: python scripts/bench_improve.py [--steps panel site mag noisy]
"""
import argparse
import os
import sys
import time
import warnings

os.environ.setdefault("LOKY_MAX_CPU_COUNT", "4")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from sklearn.linear_model import LogisticRegression, Ridge  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402
from sklearn.model_selection import StratifiedKFold, cross_val_predict  # noqa: E402
from sklearn.pipeline import make_pipeline  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bench import folds, registry, stats, tasks  # noqa: E402
from bench.latent import dawid_skene  # noqa: E402
from common import linear_probe, site_folds  # noqa: E402

warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser()
ap.add_argument("--steps", nargs="+", default=["panel", "site", "mag", "noisy"])
ap.add_argument("--boot", type=int, default=1000)
ap.add_argument("--jobs", type=int, default=4)
a = ap.parse_args()
for d in ("results/bench/x3", "results/bench/x4"):
    os.makedirs(d, exist_ok=True)
t0 = time.time()
pd.set_option("display.width", 250)


def log(*x):
    print(f"[{time.time() - t0:6.0f}s]", *x, flush=True)


L = tasks.labels()
FE = {}


def feats(e):
    if e not in FE:
        FE[e] = registry.load_means(e)
    return FE[e]


def oof(X, y, pts, seeds=range(5)):
    P = np.zeros(len(pts))
    for s in seeds:
        for tr, te in site_folds(pts, s):
            P[te] += linear_probe().fit(X[tr], y[tr]).predict_proba(X[te])[:, 1] / len(seeds)
    return P


TASKS = ["ER_IHC", "ESR1_RNApos", "PR_IHC", "PGR_prevmatch", "HER2_STRICT", "ERBB2_AMP"]
GROUP = {"ER_IHC": "ER", "ESR1_RNApos": "ER", "PR_IHC": "PR", "PGR_prevmatch": "PR", "HER2_STRICT": "HER2",
         "ERBB2_AMP": "HER2"}
GCOLS = {"ER": ["ER_IHC", "ESR1_RNApos", "ER_RPPApos"], "PR": ["PR_IHC", "PGR_prevmatch"],
         "HER2": ["HER2_STRICT", "ERBB2_AMP"]}

# ------------------------------- panel (X3) -------------------------------
if "panel" in a.steps:
    P = pd.read_csv("data/metadata/rna_panel.csv", index_col=0)
    ALPHAS = [10.0, 100.0, 1000.0, 10000.0, 100000.0]

    def ridge_panel(e, s):
        pts = P.index.intersection(feats(e).index).intersection(pd.read_csv("data/splits/folds.csv").patient).to_numpy()
        X, Y = feats(e).loc[pts].values, P.loc[pts].values
        out = np.zeros_like(Y)
        for tr, te in site_folds(pts, s):
            err = np.zeros((len(ALPHAS), Y.shape[1]))
            for itr, iva in folds.inner(pts[tr], s):
                sc = StandardScaler().fit(X[tr][itr])
                for k, al in enumerate(ALPHAS):
                    m = Ridge(alpha=al).fit(sc.transform(X[tr][itr]), Y[tr][itr])
                    err[k] += ((m.predict(sc.transform(X[tr][iva])) - Y[tr][iva]) ** 2).sum(0)
            best = err.argmin(0)
            sc = StandardScaler().fit(X[tr])
            for k in set(best):
                cols = np.where(best == k)[0]
                out[np.ix_(te, cols)] = Ridge(alpha=ALPHAS[k]).fit(sc.transform(X[tr]), Y[tr][:, cols]).predict(
                    sc.transform(X[te])).reshape(len(te), len(cols))
        return e, s, pts, out

    ENC = ["uni_v1", "hoptimus0", "titan", "conch_v1"]
    res = Parallel(n_jobs=a.jobs)(delayed(ridge_panel)(e, s) for e in ENC for s in range(5))
    rows, oo = [], []
    for e in ENC:
        parts = [r for r in res if r[0] == e]
        pts, pred = parts[0][2], np.mean([r[3] for r in parts], 0)
        Y = P.loc[pts].values
        rb = np.random.default_rng(0)
        B = [rb.integers(0, len(pts), len(pts)) for _ in range(a.boot)]
        for j, g in enumerate(P.columns):
            r = np.corrcoef(Y[:, j], pred[:, j])[0, 1]
            rbs = [np.corrcoef(Y[i, j], pred[i, j])[0, 1] for i in B]
            r2 = 1 - ((Y[:, j] - pred[:, j]) ** 2).sum() / ((Y[:, j] - Y[:, j].mean()) ** 2).sum()
            rows.append({"encoder": e, "gene": g, "n": len(pts), "pearson": r, "lo": np.percentile(rbs, 2.5),
                         "hi": np.percentile(rbs, 97.5), "r2": r2})
        oo.append(pd.DataFrame(pred, columns=P.columns).assign(encoder=e, patient=pts))
    PN = pd.DataFrame(rows)
    PN["rank_in_encoder"] = PN.groupby("encoder").pearson.rank(ascending=False)
    PN.to_csv("results/bench/x3/panel.csv", index=False)
    pd.concat(oo).to_csv("results/bench/x3/oof_panel.csv.gz", index=False)
    W = PN.pivot(index="gene", columns="encoder", values="pearson")
    W["mean"] = W.mean(1)
    log("PAM50 panel, out-of-fold Pearson r (sorted by mean):\n" + W.sort_values("mean", ascending=False).round(3).to_string())

    # multinomial PAM50 on the means (reference for multi-task ABMIL)
    PAM = ["LumA", "LumB", "Her2", "Basal", "Normal"]
    prow = []
    for e in ENC:
        d = L[["PAM50"]].dropna()
        pts = d.index.intersection(feats(e).index).to_numpy()
        y = d.loc[pts, "PAM50"].map({k: i for i, k in enumerate(PAM)}).to_numpy()
        X = feats(e).loc[pts].values
        pr = np.zeros((len(pts), 5))
        for s in range(5):
            for tr, te in site_folds(pts, s):
                m = make_pipeline(StandardScaler(), LogisticRegression(C=0.01, max_iter=3000, class_weight="balanced"))
                pr[te] += m.fit(X[tr], y[tr]).predict_proba(X[te]) / 5
        from sklearn.metrics import balanced_accuracy_score, f1_score
        prow.append({"encoder": e, "n": len(pts), "accuracy": (pr.argmax(1) == y).mean(),
                     "balanced_accuracy": balanced_accuracy_score(y, pr.argmax(1)),
                     "macro_f1": f1_score(y, pr.argmax(1), average="macro"),
                     "macro_auc_ovr": roc_auc_score(y, pr, multi_class="ovr", average="macro")})
        pd.DataFrame(pr, columns=[f"p_{k}" for k in PAM]).assign(patient=pts, y=y, encoder=e).to_csv(
            f"results/bench/x4/pam50_lr_{e}_oof.csv.gz", index=False)
    PL = pd.DataFrame(prow)
    PL.to_csv("results/bench/x4/pam50_lr.csv", index=False)
    log("PAM50 multinomial LR:\n" + PL.round(3).to_string(index=False))

# ------------------------------- site (X4c i) -------------------------------
if "site" in a.steps:
    def site_std(X, site, min_n=5):
        Z = np.empty_like(X)
        mu_all, sd_all = X.mean(0), X.std(0) + 1e-6
        for s in np.unique(site):
            m = site == s
            if m.sum() >= min_n:
                Z[m] = (X[m] - X[m].mean(0)) / (X[m].std(0) + 1e-6)
            else:
                Z[m] = (X[m] - mu_all) / sd_all
        return Z

    rows, hid = [], []
    for e in ["uni_v1", "hoptimus0", "titan", "conch_v1"]:
        F = feats(e)
        pts_all = F.index.intersection(L.index)
        Xs = pd.DataFrame(site_std(F.loc[pts_all].values, L.loc[pts_all, "site"].to_numpy(dtype=str)), index=pts_all)
        Pm = {}
        for t in TASKS:
            d = L[GCOLS[GROUP[t]]].dropna()
            pts = d.index.intersection(pts_all).to_numpy()
            y = d.loc[pts, t].astype(int).to_numpy()
            p0, p1 = oof(F.loc[pts].values, y, pts), oof(Xs.loc[pts].values, y, pts)
            Pm[t] = (pts, y, p0, p1)
            dl = stats.delong(y, p0, y, p1)
            rows.append({"encoder": e, "task": t, "n": len(pts), "auc_plain": dl["auc_ref"], "auc_site_std": dl["auc_alt"],
                         "diff": dl["diff"], "lo": dl["lo"], "hi": dl["hi"], "p": dl["p"]})
        for ihc, mol in (("ER_IHC", "ESR1_RNApos"), ("PR_IHC", "PGR_prevmatch"), ("HER2_STRICT", "ERBB2_AMP")):
            for k, tag in ((2, "plain"), (3, "site_std")):
                pi, pm = Pm[ihc], Pm[mol]
                dl = stats.delong(pi[1], pi[k], pm[1], pm[k])
                rows.append({"encoder": e, "task": f"pref {mol} - {ihc} ({tag})", "n": len(pi[0]), "diff": dl["diff"],
                             "lo": dl["lo"], "hi": dl["hi"], "p": dl["p"]})
        top = L.loc[pts_all, "site"].value_counts().index[:10]
        dd = L.loc[pts_all][L.loc[pts_all, "site"].isin(top)]
        for tag, X in (("plain", F.loc[dd.index].values), ("site_std", Xs.loc[dd.index].values)):
            pred = cross_val_predict(linear_probe(), X, dd.site.to_numpy(dtype=str),
                                     cv=StratifiedKFold(5, shuffle=True, random_state=0), n_jobs=a.jobs)
            acc = (pred == dd.site.to_numpy(dtype=str)).mean()
            lo, hi = stats.wilson(int(round(acc * len(dd))), len(dd))
            hid.append({"encoder": e, "features": tag, "n": len(dd), "accuracy": acc, "lo": lo, "hi": hi,
                        "chance": dd.site.value_counts(normalize=True).max()})
        log(f"site std {e} done")
    S, H = pd.DataFrame(rows), pd.DataFrame(hid)
    S.to_csv("results/bench/x4/site_std.csv", index=False)
    H.to_csv("results/bench/x4/site_std_hospital_id.csv", index=False)
    log("per-site standardisation:\n" + S.round(3).to_string(index=False))
    log("hospital ID:\n" + H.round(3).to_string(index=False))

# ------------------------------- mag (X4e) -------------------------------
if "mag" in a.steps:
    F20, F10 = registry.load_means("uni_v1", std=True), registry.load_means("uni_v1_10x", std=True)
    base = registry.load_means("uni_v1")
    pts_all = base.index.intersection(registry.load_means("uni_v1_10x").index)
    V = {"20x": base.loc[pts_all],
         "10x": registry.load_means("uni_v1_10x").loc[pts_all],
         "20x+10x": pd.concat([base.loc[pts_all], registry.load_means("uni_v1_10x").loc[pts_all].add_prefix("t")], axis=1)}
    if F20 is not None and F10 is not None:
        V["20x+10x (+SD)"] = pd.concat([F20.loc[pts_all], F10.loc[pts_all].add_prefix("t")], axis=1)
    else:
        V["10x+SD"] = F10.loc[pts_all] if F10 is not None else None
    V = {k: v for k, v in V.items() if v is not None}
    rows = []
    for t in TASKS:
        d = L[GCOLS[GROUP[t]]].dropna()
        pts = d.index.intersection(pts_all).to_numpy()
        y = d.loc[pts, t].astype(int).to_numpy()
        P = {k: oof(v.loc[pts].values, y, pts) for k, v in V.items()}
        for k in V:
            dl = stats.delong(y, P["20x"], y, P[k])
            rows.append({"task": t, "features": k, "n": len(pts), "auc": dl["auc_alt"], "diff_vs_20x": dl["diff"],
                         "lo": dl["lo"], "hi": dl["hi"], "p": dl["p"]})
    M = pd.DataFrame(rows)
    M.to_csv("results/bench/x4/multimag.csv", index=False)
    log("multi-magnification (UNI):\n" + M.round(3).to_string(index=False))

# ------------------------------- noisy (X4d, linear) -------------------------------
if "noisy" in a.steps:
    d = L[["ER_IHC", "ESR1_RNApos", "ER_RPPApos"]].dropna()

    def run(e, s):
        pts = d.index.intersection(feats(e).index).to_numpy()
        X = feats(e).loc[pts].values
        yi, yr, yp = (d.loc[pts, c].astype(int).to_numpy() for c in ("ER_IHC", "ESR1_RNApos", "ER_RPPApos"))
        out = {v: np.zeros(len(pts)) for v in ("plain", "consensus", "soft", "smooth")}
        for tr, te in site_folds(pts, s):
            out["plain"][te] = linear_probe().fit(X[tr], yi[tr]).predict_proba(X[te])[:, 1]
            k = yi[tr] == yr[tr]
            out["consensus"][te] = linear_probe().fit(X[tr][k], yi[tr][k]).predict_proba(X[te])[:, 1]
            for v, q in (("soft", dawid_skene(np.c_[yi[tr], yr[tr], yp[tr]])["posterior"]),
                         ("smooth", yi[tr] * 0.9 + 0.05)):
                X2 = np.r_[X[tr], X[tr]]
                y2 = np.r_[np.ones(len(tr), int), np.zeros(len(tr), int)]
                w = np.r_[q, 1 - q]
                m = make_pipeline(StandardScaler(), LogisticRegression(C=0.01, max_iter=3000))
                cw = np.where(y2 == 1, 0.5 / q.mean(), 0.5 / (1 - q.mean()))  # class balance as in the probe
                out[v][te] = m.fit(X2, y2, logisticregression__sample_weight=w * cw).predict_proba(X[te])[:, 1]
        return e, s, pts, yi, yr, yp, out

    res = Parallel(n_jobs=a.jobs)(delayed(run)(e, s) for e in ["uni_v1", "hoptimus0", "titan"] for s in range(5))
    rows = []
    for e in ["uni_v1", "hoptimus0", "titan"]:
        parts = [r for r in res if r[0] == e]
        pts, yi, yr, yp = parts[0][2:6]
        P = {v: np.mean([r[6][v] for r in parts], 0) for v in parts[0][6]}
        for v in P:
            for lab, y in (("ER_IHC", yi), ("ESR1_RNApos", yr), ("ER_RPPApos", yp)):
                dl = stats.delong(y, P["plain"], y, P[v])
                rows.append({"encoder": e, "variant": v, "scored_vs": lab, "n": len(pts), "auc": dl["auc_alt"],
                             "diff_vs_plain": dl["diff"], "lo": dl["lo"], "hi": dl["hi"], "p": dl["p"]})
    N = pd.DataFrame(rows)
    N.to_csv("results/bench/x4/noisy_linear.csv", index=False)
    log("noisy-label training (linear):\n" + N.round(3).to_string(index=False))
log("BENCH_IMPROVE_OK")
