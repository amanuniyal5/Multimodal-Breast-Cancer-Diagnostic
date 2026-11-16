"""X6 reliability (docs/DESIGN_SPEC.md §5.6) and experiment G (DX <-> frozen transfer, docs/WEEK3_EXPERIMENTS.md).

Predictor: deep ensemble = mean probability of the mean-pool linear probes on the top 3 encoders (UNI, H-optimus-0,
TITAN), fixed site-grouped folds, seed 0 folds for conformal, 5 seeds for the OOF probabilities.
Uncertainty scores: (i) ensemble predictive entropy; (ii) ensemble member SD; (iii, iv) ABMIL MC-dropout SD and
half-bag SD (results/bench/x6/reliab_uni_v1_20x_oof.csv.gz from Ada, if present); (v) |p - 0.5| of the UNI probe.
Conformal: split conformal per outer fold; calibration set = one site-grouped inner fold of the training patients
(fold_s1 column), model refitted on the remaining training patients; score 1 - p(true class); alpha = 0.10; Mondrian
(class-conditional) variant. Coverage per fold and per hospital; set size; singleton rate and accuracy; risk-coverage
(AURC) and selective AUC at 80 / 90 % coverage, vs IHC and vs RNA.
Links: Fisher odds ratio that an abstained case (set {+,-}) is IHC/RNA-discordant, ER-low, or a 1e TS/BS flip
(results/frozen_1e_pairs.csv); coverage on frozen TS slides with DX-calibrated sets.
G: UNI mean-pool models trained on DX and tested on TS (and back) vs within-slide-type, same patients, same folds:
AUC, Brier, ECE, calibration-in-the-large; RNA - IHC preference under transfer; 1b on discordant patients.
Outputs: results/bench/x6/{conformal.csv, conformal_by_site.csv, selective.csv, abstain_links.csv, uncertainty.csv,
oof_ensemble.csv.gz}, results/week3/g_transfer.csv, g_preference.csv
Usage: python scripts/bench_reliability.py [--frozen data/features/frozen_uni_v1_slides.npz]
"""
import argparse
import os
import sys
import time
import warnings

os.environ.setdefault("LOKY_MAX_CPU_COUNT", "4")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import binomtest, fisher_exact  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bench import registry, stats, tasks  # noqa: E402
from bench.metrics import classification, ece  # noqa: E402
from common import linear_probe, site_folds  # noqa: E402

warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser()
ap.add_argument("--frozen", default="data/features/frozen_uni_v1_slides.npz")
ap.add_argument("--alpha", type=float, default=0.10)
ap.add_argument("--out", default="results/bench/x6")
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
os.makedirs("results/week3", exist_ok=True)
t0 = time.time()
pd.set_option("display.width", 250)


def log(*x):
    print(f"[{time.time() - t0:6.0f}s]", *x, flush=True)


L = tasks.labels()
ENS = ["uni_v1", "hoptimus0", "titan"]
FE = {e: registry.load_means(e) for e in ENS}
TASKS = {"ER_IHC": "ESR1_RNApos", "PR_IHC": "PGR_prevmatch", "HER2_STRICT": "ERBB2_AMP"}
FOLDS = pd.read_csv("data/splits/folds.csv", index_col="patient")
pairs = pd.read_csv("results/frozen_1e_pairs.csv")
erlow = (L.ER_pct_IHC == "<10%") & (L.ER_IHC == 1)


def fitp(e, tr, te, y):
    return linear_probe().fit(FE[e].loc[tr].values, y).predict_proba(FE[e].loc[te].values)[:, 1]


def qhat(scores, alpha):
    n = len(scores)
    return np.quantile(scores, min(1.0, np.ceil((n + 1) * (1 - alpha)) / n), method="higher")


def aurc(conf, err):
    o = np.argsort(-conf)
    e = np.cumsum(err[o]) / np.arange(1, len(err) + 1)
    return e.mean()


conf_rows, site_rows, sel_rows, link_rows, unc_rows, oofs = [], [], [], [], [], []
MC = None
if os.path.exists(os.path.join(a.out, "reliab_uni_v1_20x_oof.csv.gz")):
    MC = pd.read_csv(os.path.join(a.out, "reliab_uni_v1_20x_oof.csv.gz")).groupby(["task", "patient"])[
        ["p", "mc_sd", "hb_sd", "mc_mean", "hb_mean"]].mean()
    log("ABMIL MC-dropout / half-bag results found")

for ihc, mol in TASKS.items():
    d = L[[ihc, mol, "site"]].dropna(subset=[ihc])
    pts = d.index
    for e in ENS:
        pts = pts.intersection(FE[e].index)
    pts = pts.to_numpy()
    y = d.loc[pts, ihc].astype(int).to_numpy()
    ym = d.loc[pts, mol].to_numpy(float)
    # ---- OOF ensemble over 5 seeds ----
    P = {e: np.zeros(len(pts)) for e in ENS}
    for s in range(5):
        for tr, te in site_folds(pts, s):
            for e in ENS:
                P[e][te] += fitp(e, pts[tr], pts[te], y[tr]) / 5
    pe = np.mean([P[e] for e in ENS], 0)
    sd = np.std([P[e] for e in ENS], 0)
    ent = -(pe * np.log(pe + 1e-12) + (1 - pe) * np.log(1 - pe + 1e-12))
    # ---- split conformal (seed 0 folds), standard and Mondrian ----
    pc = np.zeros(len(pts))
    sets = {"standard": np.zeros((len(pts), 2), bool), "mondrian": np.zeros((len(pts), 2), bool)}
    fold_of = np.zeros(len(pts), int)
    for k, (tr, te) in enumerate(site_folds(pts, 0)):
        fold_of[te] = k
        inner = FOLDS.loc[pts[tr], "fold_s1"].to_numpy()
        cal_fold = pd.Series(inner).value_counts().index[0]  # largest inner group = calibration set
        cal, prop = tr[inner == cal_fold], tr[inner != cal_fold]
        pcal = np.mean([fitp(e, pts[prop], pts[cal], y[prop]) for e in ENS], 0)
        pte = np.mean([fitp(e, pts[prop], pts[te], y[prop]) for e in ENS], 0)
        pc[te] = pte
        sc = np.where(y[cal] == 1, 1 - pcal, pcal)
        q = qhat(sc, a.alpha)
        sets["standard"][te] = np.c_[pte <= q, 1 - pte <= q]  # class 0 in set if 1 - p0 = p <= q; class 1 if 1 - p <= q
        for c in (0, 1):
            qc = qhat(sc[y[cal] == c], a.alpha)
            sets["mondrian"][te, c] = (pte if c == 0 else 1 - pte) <= qc
    for kind, S in sets.items():
        cover = S[np.arange(len(y)), y]
        size = S.sum(1)
        single = size == 1
        pred1 = S[:, 1] & single
        acc_single = (pred1[single].astype(int) == y[single]).mean()
        r = {"task": ihc, "kind": kind, "n": len(y), "alpha": a.alpha, "coverage": cover.mean(),
             "cov_lo": stats.wilson(int(cover.sum()), len(y))[0], "cov_hi": stats.wilson(int(cover.sum()), len(y))[1],
             "coverage_pos": cover[y == 1].mean(), "coverage_neg": cover[y == 0].mean(),
             "mean_set_size": size.mean(), "singleton_rate": single.mean(), "empty_rate": (size == 0).mean(),
             "abstain_rate": (size == 2).mean(), "acc_singletons_vs_ihc": acc_single,
             "fold_cov_min": min(cover[fold_of == k].mean() for k in np.unique(fold_of)),
             "fold_cov_max": max(cover[fold_of == k].mean() for k in np.unique(fold_of))}
        okm = ~np.isnan(ym)
        r["acc_singletons_vs_mol"] = (pred1[single & okm].astype(int) == ym[single & okm].astype(int)).mean()
        conf_rows.append(r)
        sites = d.loc[pts, "site"].to_numpy(dtype=str)
        for st in np.unique(sites):
            m = sites == st
            if m.sum() >= 15:
                site_rows.append({"task": ihc, "kind": kind, "site": st, "n": int(m.sum()), "coverage": cover[m].mean(),
                                  "abstain_rate": (size[m] == 2).mean()})
        # links: abstention vs discordance / ER-low / 1e flips
        ab = size == 2
        disc = okm & (ym != y)
        tests = [("IHC/RNA discordant", disc, okm)]
        if ihc == "ER_IHC":
            tests.append(("ER-low (IHC+ <10%)", erlow.reindex(pts).fillna(False).to_numpy(), np.ones(len(y), bool)))
        pr = pairs[pairs.label == ihc].set_index("patient")
        flip = ((pr.P_TS >= 0.5) != (pr.P_BS >= 0.5))
        inpair = np.isin(pts, flip.index)
        tests.append(("1e TS/BS flip", np.isin(pts, flip[flip].index), inpair))
        for name, flag, base in tests:
            t2 = np.array([[np.sum(ab & flag & base), np.sum(ab & ~flag & base)],
                           [np.sum(~ab & flag & base), np.sum(~ab & ~flag & base)]])
            orr, pv = fisher_exact(t2)
            link_rows.append({"task": ihc, "kind": kind, "event": name, "n": int(base.sum()),
                              "abstained": int(t2[0].sum()), "events": int(t2[:, 0].sum()),
                              "events_in_abstained": int(t2[0, 0]),
                              "rate_abstained": t2[0, 0] / max(t2[0].sum(), 1),
                              "rate_kept": t2[1, 0] / max(t2[1].sum(), 1), "odds_ratio": orr, "p_fisher": pv})
    # ---- selective prediction: risk-coverage for each uncertainty score ----
    U = {"ensemble |p-0.5|": np.abs(pe - 0.5), "ensemble entropy (neg)": -ent, "ensemble member SD (neg)": -sd,
         "UNI probe |p-0.5|": np.abs(P["uni_v1"] - 0.5)}
    if MC is not None and ihc in MC.index.get_level_values(0):
        m = MC.loc[ihc].reindex(pts)
        ok = m.p.notna().to_numpy()
        if ok.all():
            U["ABMIL |p-0.5|"] = np.abs(m.p.to_numpy() - 0.5)
            U["ABMIL MC-dropout SD (neg)"] = -m.mc_sd.to_numpy()
            U["ABMIL half-bag SD (neg)"] = -m.hb_sd.to_numpy()
    for name, conf in U.items():
        base_p = pe if not name.startswith("ABMIL") else MC.loc[ihc].reindex(pts).p.to_numpy()
        for lab, yy, ok in (("IHC", y, np.ones(len(y), bool)), ("RNA", ym, ~np.isnan(ym))):
            yy2, c2, p2 = yy[ok].astype(int), conf[ok], base_p[ok]
            err = ((p2 >= 0.5).astype(int) != yy2).astype(float)
            o = np.argsort(-c2)
            r = {"task": ihc, "score": name, "vs": lab, "n": len(yy2), "aurc": aurc(c2, err), "error_all": err.mean(),
                 "auc_all": roc_auc_score(yy2, p2)}
            for cov in (0.9, 0.8, 0.7):
                k = o[:int(round(cov * len(o)))]
                r[f"error@{int(cov * 100)}"] = err[k].mean()
                r[f"auc@{int(cov * 100)}"] = roc_auc_score(yy2[k], p2[k]) if len(set(yy2[k])) == 2 else np.nan
            # discordance among the 10 % least confident
            low = o[-int(round(0.1 * len(o))):]
            if lab == "IHC":
                dd = (~np.isnan(ym[ok])) & (ym[ok] != yy2)
                r["discordant_in_lowest10%"] = int(dd[low].sum())
                r["discordant_total"] = int(dd.sum())
                r["p_enrich"] = binomtest(int(dd[low].sum()), len(low), dd.mean(), alternative="greater").pvalue
            sel_rows.append(r)
    for name, conf in U.items():
        unc_rows.append({"task": ihc, "score": name, "auc_detect_discordance": roc_auc_score(
            ((~np.isnan(ym)) & (ym != y))[~np.isnan(ym)], -conf[~np.isnan(ym)])})
    oofs.append(pd.DataFrame({"task": ihc, "patient": pts, "y": y, "y_mol": ym, "p_ens": pe, "sd_ens": sd,
                              "p_conformal_model": pc, "set_standard_neg": sets["standard"][:, 0],
                              "set_standard_pos": sets["standard"][:, 1], "set_mondrian_neg": sets["mondrian"][:, 0],
                              "set_mondrian_pos": sets["mondrian"][:, 1],
                              **{f"p_{e}": P[e] for e in ENS}}))
    log(f"{ihc} done: n {len(y)}")

C = pd.DataFrame(conf_rows)
C.to_csv(os.path.join(a.out, "conformal.csv"), index=False)
pd.DataFrame(site_rows).to_csv(os.path.join(a.out, "conformal_by_site.csv"), index=False)
SL = pd.DataFrame(sel_rows)
SL.to_csv(os.path.join(a.out, "selective.csv"), index=False)
LK = pd.DataFrame(link_rows)
LK.to_csv(os.path.join(a.out, "abstain_links.csv"), index=False)
UN = pd.DataFrame(unc_rows)
UN.to_csv(os.path.join(a.out, "uncertainty.csv"), index=False)
pd.concat(oofs).to_csv(os.path.join(a.out, "oof_ensemble.csv.gz"), index=False)
log("conformal:\n" + C.round(3).to_string(index=False))
SB = pd.DataFrame(site_rows)
log("per-hospital coverage (sites with >= 15 patients): " + SB.groupby(["task", "kind"]).coverage.agg(
    ["min", "max", "count"]).round(3).to_string())
log("selective:\n" + SL.round(3).to_string(index=False))
log("abstention links:\n" + LK.round(3).to_string(index=False))
log("uncertainty as discordance detector:\n" + UN.round(3).to_string(index=False))

# ============================ G: DX <-> frozen transfer (UNI) ============================
if os.path.exists(a.frozen):
    z = np.load(a.frozen, allow_pickle=True)
    fr = pd.DataFrame(z["X"], index=z["slides"].astype(str))
    meta = pd.DataFrame({"patient": z["patients"].astype(str), "section": z["section"].astype(str),
                         "primary_TS": z["primary_TS"]}, index=fr.index)
    ts = meta[(meta.section == "TS") & meta.primary_TS]
    TS = fr.loc[ts.index].groupby(ts.patient.values).mean()
    DX = registry.load_means("uni_v1")
    GT = {"ER_IHC": "ESR1_RNApos", "PR_IHC": "PGR_prevmatch", "HER2_STRICT": "ERBB2_AMP"}
    rows, pref, fz = [], [], []
    for ihc, mol in GT.items():
        d = L[[ihc, mol]].dropna()
        pts = d.index.intersection(TS.index).intersection(DX.index).to_numpy()
        Y = {c: d.loc[pts, c].astype(int).to_numpy() for c in (ihc, mol)}
        PP = {}
        for c in (ihc, mol):
            for tr_name, Ftr in (("DX", DX), ("TS", TS)):
                for te_name, Fte in (("DX", DX), ("TS", TS)):
                    p = np.zeros(len(pts))
                    for s in range(5):
                        for tr, te in site_folds(pts, s):
                            p[te] += linear_probe().fit(Ftr.loc[pts[tr]].values, Y[c][tr]).predict_proba(
                                Fte.loc[pts[te]].values)[:, 1] / 5
                    PP[(c, tr_name, te_name)] = p
                    m = classification(Y[c], p)
                    rows.append({"label": c, "train": tr_name, "test": te_name, "n": len(pts), **m,
                                 "mean_p": p.mean(), "prevalence": Y[c].mean()})
                    fz.append(pd.DataFrame({"label": c, "train": tr_name, "test": te_name, "patient": pts,
                                            "y": Y[c], "p": p}))
        for tr_name, te_name in (("DX", "DX"), ("TS", "TS"), ("DX", "TS"), ("TS", "DX")):
            dl = stats.delong(Y[ihc], PP[(ihc, tr_name, te_name)], Y[mol], PP[(mol, tr_name, te_name)])
            pref.append({"pair": f"{mol} - {ihc}", "train": tr_name, "test": te_name, "n": len(pts), "delta": dl["diff"],
                         "lo": dl["lo"], "hi": dl["hi"], "p": dl["p"]})
        for c in (ihc, mol):  # transfer cost on the same test slides: DX->TS vs TS->TS, TS->DX vs DX->DX
            for ref, alt in ((("TS", "TS"), ("DX", "TS")), (("DX", "DX"), ("TS", "DX"))):
                dl = stats.delong(Y[c], PP[(c, *ref)], Y[c], PP[(c, *alt)])
                pref.append({"pair": f"transfer {c}", "train": f"{alt[0]} vs {ref[0]}", "test": alt[1], "n": len(pts),
                             "delta": dl["diff"], "lo": dl["lo"], "hi": dl["hi"], "p": dl["p"]})
        # 1b under transfer: train on concordant DX patients, score discordant patients' TS slides
        dd = L[[ihc, mol]].dropna()
        dd = dd[dd.index.isin(DX.index)]
        conc, disc = dd[dd[ihc] == dd[mol]], dd[(dd[ihc] != dd[mol]) & dd.index.isin(TS.index)]
        m = linear_probe().fit(DX.loc[conc.index].values, conc[ihc].astype(int))
        for te_name, F in (("DX", DX), ("TS", TS)):
            k = int(((m.predict_proba(F.loc[disc.index].values)[:, 1] >= 0.5).astype(int) == disc[mol].astype(int)).sum())
            pref.append({"pair": f"1b {ihc} vs {mol} (DX-trained)", "train": "DX", "test": te_name, "n": len(disc),
                         "delta": k / len(disc), "lo": np.nan, "hi": np.nan, "p": binomtest(k, len(disc)).pvalue,
                         "side_molecular": k, "side_IHC": len(disc) - k})
    G, GP = pd.DataFrame(rows), pd.DataFrame(pref)
    G.to_csv("results/week3/g_transfer.csv", index=False)
    GP.to_csv("results/week3/g_preference.csv", index=False)
    pd.concat(fz).to_csv("results/week3/g_oof.csv.gz", index=False)
    log("G transfer:\n" + G.round(3).to_string(index=False))
    log("G preference / transfer cost / 1b:\n" + GP.round(3).to_string(index=False))

    # X6 link: DX-calibrated conformal sets applied to frozen TS slides (ER, UNI probe only, seed 0 folds)
    d = L[["ER_IHC", "site"]].dropna()
    pts = d.index.intersection(DX.index).to_numpy()
    y = d.loc[pts, "ER_IHC"].astype(int).to_numpy()
    cov_dx, cov_ts, size_dx, size_ts = [], [], [], []
    for tr, te in site_folds(pts, 0):
        inner = FOLDS.loc[pts[tr], "fold_s1"].to_numpy()
        cal_fold = pd.Series(inner).value_counts().index[0]
        cal, prop = tr[inner == cal_fold], tr[inner != cal_fold]
        mdl = linear_probe().fit(DX.loc[pts[prop]].values, y[prop])
        pcal = mdl.predict_proba(DX.loc[pts[cal]].values)[:, 1]
        q = qhat(np.where(y[cal] == 1, 1 - pcal, pcal), a.alpha)
        for F, cov, size in ((DX, cov_dx, size_dx), (TS, cov_ts, size_ts)):
            tt = te[np.isin(pts[te], F.index)]
            tt = tt[np.isin(pts[tt], TS.index)]  # same patients for both slide types
            p = mdl.predict_proba(F.loc[pts[tt]].values)[:, 1]
            S = np.c_[p <= q, 1 - p <= q]
            cov.extend(S[np.arange(len(tt)), y[tt]])
            size.extend(S.sum(1))
    tsr = pd.DataFrame([{"slides": k, "n": len(c), "coverage": np.mean(c), "mean_set_size": np.mean(s),
                         "abstain_rate": np.mean(np.array(s) == 2)}
                        for k, c, s in (("DX", cov_dx, size_dx), ("TS (DX-calibrated)", cov_ts, size_ts))])
    tsr.to_csv(os.path.join(a.out, "conformal_frozen.csv"), index=False)
    log("conformal on frozen TS with DX calibration (ER, UNI probe):\n" + tsr.round(3).to_string(index=False))
else:
    log(f"{a.frozen} not found: G and the frozen conformal check skipped")
log("BENCH_RELIABILITY_OK")
