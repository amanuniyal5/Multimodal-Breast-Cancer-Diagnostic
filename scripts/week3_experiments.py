"""Week 3 experiments C, D, E, H (docs/WEEK3_EXPERIMENTS.md) on patient means and saved out-of-fold predictions.
Fixed site-grouped folds (common.site_folds), mean-pool linear probe (C = 0.01, class-balanced), 5 seeds averaged.
Encoders (primary set, no TCGA-pretrained contrasts): UNI, H-optimus-0, TITAN, CONCH, UNI2-h, GigaPath, Hibou-L.

C  latent-class "truth" (Dawid-Skene, two classes, conditional independence):
   ER: IHC, RNA (GMM), RPPA (GMM); HER2: IHC/FISH (strict), GISTIC amplification, ERBB2 RNA, HER2 RPPA.
   Per-test sensitivity / specificity (patient bootstrap, 500), each H&E classifier's AUC against the latent truth vs
   against IHC (generalised DeLong, same patients), and accuracy of IHC vs the H&E model at matched prevalence.
   Sensitivity: the 4-rater model that also includes the UNI H&E call.          -> c_tests.csv, c_rescore.csv
D  H&E model as an IHC auditor, RNA as the oracle: score = P(opposite of the IHC call); target = IHC != RNA.
   AUC, precision / recall at the top 1, 2, 5, 10 % of the cohort, ER-low enrichment; baseline = IHC % stained.
   Then drop (or flip) the training patients the auditor flags (top 5 % of the training fold, flagged from inner-fold
   predictions only) and retrain: AUC vs IHC and vs RNA, DeLong against plain IHC training. -> d_audit.csv, d_retrain.csv
E  which label to train on: RNA-trained vs IHC-trained models, both scored against IHC (and RNA), same patients;
   ER, PR (prevalence-matched PGR), HER2 (amplification); learning curves n = 50 ... all (seed 0 folds, 5 repeats).
                                                                                  -> e_train_label.csv, e_curves.csv
H  leaderboard under each reference standard: 7 encoders (mean pool) + 5 MIL heads + 5 pooling baselines (UNI, X2)
   ranked by AUC under IHC / RNA / RPPA; Kendall tau between rankings, bootstrap P(top-1 changes) and
   P(rank of each model changes).                                                -> h_ranks.csv, h_tau.csv
Usage: python scripts/week3_experiments.py [--steps C D E H] [--out results/week3]
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
from scipy.stats import binomtest, kendalltau  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bench import folds, registry, stats, tasks  # noqa: E402
from bench.latent import dawid_skene  # noqa: E402
from common import linear_probe, site_folds  # noqa: E402

warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser()
ap.add_argument("--steps", nargs="+", default=["C", "D", "E", "H"])
ap.add_argument("--out", default="results/week3")
ap.add_argument("--x1", default="results/bench/x1/oof.csv.gz")
ap.add_argument("--x2", default="results/bench/x2")
ap.add_argument("--boot", type=int, default=500)
ap.add_argument("--jobs", type=int, default=4)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
t0 = time.time()
pd.set_option("display.width", 250)


def log(*x):
    print(f"[{time.time() - t0:6.0f}s]", *x, flush=True)


ENCS = ["uni_v1", "hoptimus0", "titan", "conch_v1", "uni_v2", "gigapath", "hibou_l"]
L = tasks.labels()
O1 = pd.read_csv(a.x1)
FE = {}


def feats(e):
    if e not in FE:
        FE[e] = registry.load_means(e)
    return FE[e]


def oof_lr(X, y, pts, seeds=range(5), w=None):
    """seed-averaged OOF probability of the linear probe (sample weights optional)"""
    P = np.zeros(len(pts))
    for s in seeds:
        for tr, te in site_folds(pts, s):
            m = linear_probe()
            kw = {} if w is None else {"logisticregression__sample_weight": w[tr]}
            P[te] += m.fit(X[tr], y[tr], **kw).predict_proba(X[te])[:, 1] / len(seeds)
    return P


# =============================== C: latent class ===============================
if "C" in a.steps:
    rows, resc = [], []
    SETS = {"ER": ["ER_IHC", "ESR1_RNApos", "ER_RPPApos"],
            "HER2": ["HER2_STRICT", "ERBB2_AMP", "ERBB2_RNApos", "HER2_RPPApos"]}
    for marker, cols in SETS.items():
        d = L[cols].dropna()
        R = d[cols].to_numpy(float)
        fit = dawid_skene(R)
        r = np.random.default_rng(0)
        B = [dawid_skene(R[r.integers(0, len(R), len(R))]) for _ in range(a.boot)]
        for j, c in enumerate(cols):
            se_b, sp_b = np.array([b["sens"][j] for b in B]), np.array([b["spec"][j] for b in B])
            rows.append({"marker": marker, "model": "tests only", "test": c, "n": len(d), "prevalence": fit["prevalence"],
                         "sens": fit["sens"][j], "sens_lo": np.percentile(se_b, 2.5), "sens_hi": np.percentile(se_b, 97.5),
                         "spec": fit["spec"][j], "spec_lo": np.percentile(sp_b, 2.5),
                         "spec_hi": np.percentile(sp_b, 97.5), "acc_vs_latent": np.mean((R[:, j] == 1) == (fit["posterior"] > 0.5))})
        post = pd.Series(fit["posterior"], index=d.index)
        post.rename("posterior").to_frame().assign(marker=marker).to_csv(os.path.join(a.out, f"c_posterior_{marker}.csv"))
        truth = (post > 0.5).astype(int)
        log(f"C {marker}: n {len(d)}, prevalence {fit['prevalence']:.3f}; uncertain (0.2-0.8) "
            f"{((post > 0.2) & (post < 0.8)).sum()}; sens {np.round(fit['sens'], 3)}, spec {np.round(fit['spec'], 3)}")
        ihc, mol = cols[0], cols[1]
        for e in ENCS:
            for trained in (ihc, mol):
                o = O1[(O1.encoder == e) & (O1.task == trained)].set_index("patient")
                pts = o.index.intersection(d.index)
                if len(pts) < 50:
                    continue
                yi, yt, p = d.loc[pts, ihc].astype(int), truth.loc[pts], o.loc[pts, "p"]
                dl = stats.delong(yi.values, p.values, yt.values, p.values)
                conf = (post.loc[pts] < 0.1) | (post.loc[pts] > 0.9)
                # matched-prevalence binary call of the H&E model vs the IHC test, both scored against the latent truth
                thr = np.quantile(p, 1 - yi.mean())
                acc_he = np.mean((p.values >= thr) == (yt.values == 1))
                acc_ihc = np.mean(yi.values == yt.values)
                rb = np.random.default_rng(1)
                db = []
                for _ in range(a.boot):
                    i = rb.integers(0, len(pts), len(pts))
                    db.append(np.mean((p.values[i] >= thr) == (yt.values[i] == 1)) - np.mean(yi.values[i] == yt.values[i]))
                resc.append({"marker": marker, "encoder": e, "trained_on": trained, "n": len(pts),
                             "auc_vs_ihc": dl["auc_ref"], "auc_vs_latent": dl["auc_alt"], "delta": dl["diff"],
                             "lo": dl["lo"], "hi": dl["hi"], "p": dl["p"],
                             "auc_vs_latent_confident": roc_auc_score(yt[conf], p[conf]), "n_confident": int(conf.sum()),
                             "acc_HE_matched": acc_he, "acc_IHC": acc_ihc, "acc_diff": acc_he - acc_ihc,
                             "acc_diff_lo": np.percentile(db, 2.5), "acc_diff_hi": np.percentile(db, 97.5)})
        # sensitivity: 4 raters incl. the UNI IHC-trained H&E call (prevalence-matched threshold)
        o = O1[(O1.encoder == "uni_v1") & (O1.task == ihc)].set_index("patient")
        pts = o.index.intersection(d.index)
        he = (o.loc[pts, "p"] >= np.quantile(o.loc[pts, "p"], 1 - d.loc[pts, ihc].mean())).astype(float)
        f4 = dawid_skene(np.c_[d.loc[pts, cols].to_numpy(float), he.to_numpy()])
        for j, c in enumerate(cols + ["H&E (UNI, IHC-trained)"]):
            rows.append({"marker": marker, "model": "tests + H&E", "test": c, "n": len(pts), "prevalence": f4["prevalence"],
                         "sens": f4["sens"][j], "spec": f4["spec"][j]})
    C = pd.DataFrame(rows)
    C.to_csv(os.path.join(a.out, "c_tests.csv"), index=False)
    RS = pd.DataFrame(resc)
    RS.to_csv(os.path.join(a.out, "c_rescore.csv"), index=False)
    log("C tests:\n" + C.round(3).to_string(index=False))
    log("C rescore:\n" + RS.round(3).to_string(index=False))

# =============================== D: label auditor ===============================
if "D" in a.steps:
    d = L[["ER_IHC", "ESR1_RNApos", "ER_pct_IHC", "PAM50", "site"]].dropna(subset=["ER_IHC", "ESR1_RNApos"])
    disc_all = (d.ER_IHC != d.ESR1_RNApos).astype(int)
    rows, retr, scores = [], [], {}
    common = d.index
    for e in ENCS:
        common = common.intersection(feats(e).index)
    d, disc_all = d.loc[common], disc_all.loc[common]
    pts = d.index.to_numpy()
    yi, yr = d.ER_IHC.astype(int).to_numpy(), d.ESR1_RNApos.astype(int).to_numpy()
    disc = (yi != yr).astype(int)
    erlow = ((d.ER_pct_IHC == "<10%") & (d.ER_IHC == 1)).to_numpy()
    log(f"D: {len(pts)} patients with IHC + RNA + all encoders; discordant {disc.sum()} ({disc.mean():.1%}); ER-low {erlow.sum()}")
    P = dict(zip(ENCS, Parallel(n_jobs=a.jobs)(delayed(oof_lr)(feats(e).loc[pts].values, yi, pts) for e in ENCS)))
    P["consensus (7 encoders)"] = np.mean([P[e] for e in ENCS], 0)
    P["consensus (UNI, H-optimus-0, TITAN)"] = np.mean([P[e] for e in ["uni_v1", "hoptimus0", "titan"]], 0)
    pct = d.ER_pct_IHC.map({"<10%": 5, "10-19%": 15, "20-29%": 25, "30-39%": 35, "40-49%": 45, "50-59%": 55,
                            "60-69%": 65, "70-79%": 75, "80-89%": 85, "90-99%": 95}).to_numpy(float)
    for name, p in P.items():
        sc = np.where(yi == 1, 1 - p, p)  # probability of the opposite of the IHC call
        scores[name] = sc
        auc, lo, hi = stats.auc_ci(disc, sc)
        r = {"auditor": name, "n": len(pts), "n_disc": int(disc.sum()), "base_rate": disc.mean(), "auc": auc, "lo": lo,
             "hi": hi, "auc_ihcpos": roc_auc_score(disc[yi == 1], sc[yi == 1]),
             "auc_ihcneg": roc_auc_score(disc[yi == 0], sc[yi == 0])}
        order = np.argsort(-sc)
        for k in (1, 2, 5, 10):
            top = order[:int(round(len(pts) * k / 100))]
            hit = disc[top].sum()
            r.update({f"prec@{k}%": hit / len(top), f"prec@{k}%_lo": stats.wilson(hit, len(top))[0],
                      f"prec@{k}%_hi": stats.wilson(hit, len(top))[1], f"recall@{k}%": hit / disc.sum(),
                      f"p_binom@{k}%": binomtest(int(hit), len(top), disc.mean(), alternative="greater").pvalue,
                      f"erlow@{k}%": int(erlow[top].sum())})
        rows.append(r)
    # baseline: IHC % stained (IHC+ patients with a recorded %): lower % = more suspicious
    m = (yi == 1) & ~np.isnan(pct)
    rows.append({"auditor": "IHC % stained (IHC+ with %)", "n": int(m.sum()), "n_disc": int(disc[m].sum()),
                 "base_rate": disc[m].mean(), "auc": roc_auc_score(disc[m], -pct[m])})
    for name in ["uni_v1", "consensus (7 encoders)"]:
        rows.append({"auditor": f"{name} on the same IHC+ with %", "n": int(m.sum()), "n_disc": int(disc[m].sum()),
                     "base_rate": disc[m].mean(), "auc": roc_auc_score(disc[m], scores[name][m])})
    DA = pd.DataFrame(rows)
    DA.to_csv(os.path.join(a.out, "d_audit.csv"), index=False)
    pd.DataFrame({"patient": pts, "ER_IHC": yi, "ESR1_RNApos": yr, "ER_low": erlow, "PAM50": d.PAM50.values,
                  **{f"score_{k}": v for k, v in scores.items()}}).to_csv(os.path.join(a.out, "d_scores.csv.gz"),
                                                                           index=False)
    log("D audit:\n" + DA.round(3).to_string(index=False))

    def retrain(e, s, variant, frac=0.05):
        X = feats(e).loc[pts].values
        out = np.zeros(len(pts))
        for tr, te in site_folds(pts, s):
            ytr = yi[tr].copy()
            keep = np.ones(len(tr), bool)
            if variant != "plain":
                pin = np.zeros(len(tr))  # inner-fold OOF on the training patients only
                for itr, iva in folds.inner(pts[tr], s):
                    pin[iva] = linear_probe().fit(X[tr][itr], ytr[itr]).predict_proba(X[tr][iva])[:, 1]
                sc = np.where(ytr == 1, 1 - pin, pin)
                flag = np.argsort(-sc)[:int(round(frac * len(tr)))]
                if variant == "drop":
                    keep[flag] = False
                else:
                    ytr[flag] = 1 - ytr[flag]
            out[te] = linear_probe().fit(X[tr][keep], ytr[keep]).predict_proba(X[te])[:, 1]
        return e, s, variant, out

    jobs = [(e, s, v) for e in ["uni_v1", "hoptimus0", "titan"] for s in range(5) for v in ("plain", "drop", "flip")]
    res = Parallel(n_jobs=a.jobs)(delayed(retrain)(*j) for j in jobs)
    PR = {}
    for e, s, v, p in res:
        PR.setdefault((e, v), []).append(p)
    for e in ["uni_v1", "hoptimus0", "titan"]:
        base = np.mean(PR[(e, "plain")], 0)
        for v in ("plain", "drop", "flip"):
            p = np.mean(PR[(e, v)], 0)
            for lab, y in (("ER_IHC", yi), ("ESR1_RNApos", yr)):
                dl = stats.delong(y, base, y, p)
                retr.append({"encoder": e, "variant": v, "scored_vs": lab, "n": len(pts), "auc": dl["auc_alt"],
                             "diff_vs_plain": dl["diff"], "lo": dl["lo"], "hi": dl["hi"], "p": dl["p"]})
    DR = pd.DataFrame(retr)
    DR.to_csv(os.path.join(a.out, "d_retrain.csv"), index=False)
    log("D retrain:\n" + DR.round(4).to_string(index=False))

# =============================== E: which label to train on ===============================
if "E" in a.steps:
    PAIRS = [("ER_IHC", "ESR1_RNApos"), ("PR_IHC", "PGR_prevmatch"), ("HER2_STRICT", "ERBB2_AMP")]
    rows, curves = [], []
    for ihc, mol in PAIRS:
        d = L[[ihc, mol]].dropna()
        for e in ["uni_v1", "hoptimus0", "titan", "conch_v1", "uni_v2"]:
            pts = d.index.intersection(feats(e).index).to_numpy()
            X = feats(e).loc[pts].values
            yi, ym = d.loc[pts, ihc].astype(int).values, d.loc[pts, mol].astype(int).values
            pi, pm = oof_lr(X, yi, pts), oof_lr(X, ym, pts)
            for lab, y in ((ihc, yi), (mol, ym)):
                dl = stats.delong(y, pi, y, pm)  # same label, two models: RNA-trained minus IHC-trained
                rows.append({"pair": f"{ihc} vs {mol}", "encoder": e, "scored_vs": lab, "n": len(pts),
                             "auc_ihc_trained": dl["auc_ref"], "auc_mol_trained": dl["auc_alt"], "diff": dl["diff"],
                             "lo": dl["lo"], "hi": dl["hi"], "p": dl["p"]})

            def curve(n, rep):
                out = []
                for tr, te in site_folds(pts, 0):
                    for lab, ytrain in (("ihc", yi), ("mol", ym)):
                        sub = folds.subsample(ytrain[tr], n, rep, 0)
                        p = linear_probe().fit(X[tr][sub], ytrain[tr][sub]).predict_proba(X[te])[:, 1]
                        out.append((lab, te, p))
                res = {}
                for lab, te, p in out:
                    res.setdefault(lab, np.zeros(len(pts)))[te] = p
                return n, rep, roc_auc_score(yi, res["ihc"]), roc_auc_score(yi, res["mol"])

            for n, rep, ai, am in Parallel(n_jobs=a.jobs)(delayed(curve)(n, r) for n in (50, 100, 200, 400, 10 ** 6)
                                                          for r in range(5 if n < 10 ** 6 else 1)):
                curves.append({"pair": f"{ihc} vs {mol}", "encoder": e, "n_train": n if n < 10 ** 6 else "all",
                               "rep": rep, "auc_vs_ihc_ihc_trained": ai, "auc_vs_ihc_mol_trained": am})
            log(f"E {ihc} vs {mol} {e} done")
    E = pd.DataFrame(rows)
    E.to_csv(os.path.join(a.out, "e_train_label.csv"), index=False)
    EC = pd.DataFrame(curves)
    EC.to_csv(os.path.join(a.out, "e_curves.csv"), index=False)
    log("E:\n" + E.round(3).to_string(index=False))
    log("E curves (mean over reps):\n" + EC.groupby(["pair", "encoder", "n_train"], sort=False)[
        ["auc_vs_ihc_ihc_trained", "auc_vs_ihc_mol_trained"]].mean().round(3).to_string())

# =============================== H: leaderboard under each reference standard ===============================
if "H" in a.steps:
    import glob
    GROUPS = {"ER": ["ER_IHC", "ESR1_RNApos", "ER_RPPApos"], "PR": ["PR_IHC", "PGR_prevmatch"],
              "HER2": ["HER2_STRICT", "ERBB2_AMP"]}
    X2 = pd.concat([pd.read_csv(f) for f in glob.glob(os.path.join(a.x2, "*_uni_v1_20x_oof.csv.gz"))])
    X2 = X2[X2.exp == "1a"].groupby(["method", "task", "patient"]).agg(y=("y", "first"), p=("p", "mean")).reset_index()
    rank_rows, tau_rows = [], []
    for g, cols in GROUPS.items():
        M = {}
        for e in ENCS:
            for c in cols:
                o = O1[(O1.encoder == e) & (O1.task == c)].set_index("patient")
                M.setdefault(f"{registry.label(e)} (mean)", {})[c] = o
        for m in X2.method.unique():
            for c in cols:
                M.setdefault(f"UNI {m}", {})[c] = X2[(X2.method == m) & (X2.task == c)].set_index("patient")
        pts = None
        for m in M:
            for c in cols:
                pts = M[m][c].index if pts is None else pts.intersection(M[m][c].index)
        pts = pts.to_numpy()
        names = list(M)
        Y = {c: M[names[0]][c].loc[pts, "y"].astype(int).to_numpy() for c in cols}
        PP = {(m, c): M[m][c].loc[pts, "p"].to_numpy() for m in names for c in cols}

        def table(idx):
            return pd.DataFrame({c: [roc_auc_score(Y[c][idx], PP[(m, c)][idx]) for m in names] for c in cols}, index=names)

        T = table(np.arange(len(pts)))
        R = T.rank(ascending=False)
        rb = np.random.default_rng(0)
        flips, taus, top1 = {c: [] for c in cols[1:]}, {c: [] for c in cols[1:]}, {c: [] for c in cols[1:]}
        for _ in range(a.boot):
            Tb = table(rb.integers(0, len(pts), len(pts)))
            for c in cols[1:]:
                taus[c].append(kendalltau(Tb[cols[0]], Tb[c])[0])
                top1[c].append(Tb[cols[0]].idxmax() != Tb[c].idxmax())
        for c in cols[1:]:
            tau = kendalltau(T[cols[0]], T[c])[0]
            tau_rows.append({"marker": g, "ref": cols[0], "alt": c, "n": len(pts), "models": len(names), "kendall_tau": tau,
                             "tau_lo": np.percentile(taus[c], 2.5), "tau_hi": np.percentile(taus[c], 97.5),
                             "top1_ref": T[cols[0]].idxmax(), "top1_alt": T[c].idxmax(),
                             "P_top1_differs": np.mean(top1[c]),
                             "max_rank_shift": int((R[cols[0]] - R[c]).abs().max())})
        for m in names:
            rank_rows.append({"marker": g, "model": m, **{f"auc_{c}": T.at[m, c] for c in cols},
                              **{f"rank_{c}": R.at[m, c] for c in cols}})
        log(f"H {g}: n {len(pts)}, {len(names)} models")
    HR, HT = pd.DataFrame(rank_rows), pd.DataFrame(tau_rows)
    HR.to_csv(os.path.join(a.out, "h_ranks.csv"), index=False)
    HT.to_csv(os.path.join(a.out, "h_tau.csv"), index=False)
    log("H ranks:\n" + HR.round(3).to_string(index=False))
    log("H tau:\n" + HT.round(3).to_string(index=False))
log("WEEK3_OK")
