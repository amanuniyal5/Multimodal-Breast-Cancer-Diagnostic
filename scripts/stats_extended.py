"""Week 3 statistics beyond p-values, from saved predictions only (no training; runs on the laptop).

1a  Bayesian bootstrap (Dirichlet weights over patients, 2,000 draws) of ΔAUC = AUC(molecular) − AUC(IHC) for every
    model variant: posterior median, 95% credible interval, P(ΔAUC > 0). Calibration of each label's own model
    (seed-averaged out-of-fold P): Brier score, expected calibration error (10 equal-width bins), calibration intercept
    and slope (logistic regression of y on logit P; ideal 0 and 1).
1b  For each "which side" count (k of n discordant patients called on the molecular side): Wilson and Clopper–Pearson
    95% CIs, Cohen's h vs 50:50, Beta(1, 1) posterior (mean, 95% credible interval, P(θ > 0.5)), Bayes factors
    BF10 (θ ~ U(0, 1) vs θ = 0.5) and BF+0 (θ ~ U(0.5, 1) vs θ = 0.5). Paired DX vs TS calls on the same discordant
    patients: exact McNemar test.
1c  Equivalence (TOST) for the slide × label interaction with a stated margin (default ±0.05 AUC): bootstrap 90% CI
    (2,000 patient resamples), TOST p, and the smallest margin at which equivalence holds; permutation test
    (slide types swapped within patient, both labels together; 2,000 permutations).
1e  TS vs BS agreement beyond kappa: % flipped with Wilson CI, PABAK and Gwet's AC1 (robust to prevalence) with
    bootstrap CIs, exact McNemar test for a systematic TS/BS shift, Bland–Altman bias and 95% limits of agreement
    on logit P.

Frozen-slide inputs: results/frozen_1c_oof.csv.gz / frozen_1e_pairs.csv (UNI mean-pool) and, when present,
results/frozen_abmil_1c_oof.csv.gz / frozen_abmil_1e_pairs.csv (UNI ABMIL, scripts/frozen_abmil.py), same columns.

Usage: python scripts/stats_extended.py [--results results] [--boot 2000] [--margin 0.05]
Outputs: results/stats_ext_1a.csv, stats_ext_1b.csv, stats_ext_1c.csv, stats_ext_1e.csv
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy.special import betaln
from scipy.stats import beta, binomtest, norm
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

ap = argparse.ArgumentParser()
ap.add_argument("--results", default="results")
ap.add_argument("--boot", type=int, default=2000)
ap.add_argument("--margin", type=float, default=0.05, help="TOST equivalence margin for the 1c interaction (AUC)")
a = ap.parse_args()
R = lambda f: os.path.join(a.results, f)  # noqa: E731

# fixed-split variants, as in stats_delong_1b.py --abmil fixed
VARIANTS = [("TITAN linear", "titan_linear_fixedsplits_oof.csv.gz"),
            ("UNI mean", "mil_uni_v1_20x_fixedsplits_oof.csv.gz"),
            ("CONCH mean", "mil_conch_v1_20x_fixedsplits_oof.csv.gz"),
            ("UNI ABMIL", "mil_uni_v1_20x_fixedsplits_abmil_oof.csv.gz"),
            ("CONCH ABMIL", "mil_conch_v1_20x_fixedsplits_abmil_oof.csv.gz")]
COMPARISONS_1A = [("ER", "ER_IHC", "ESR1_RNApos"), ("ER", "ER_IHC", "ER_RPPApos"), ("PR", "PR_IHC", "PGR_RNApos"),
                  ("HER2", "HER2_STRICT", "ERBB2_AMP"), ("HER2", "HER2_STRICT", "HER2_RPPApos")]
# PR with the prevalence-matched PGR cut-off (under which PR shows no preference)
VARIANTS_PRPREV = [("TITAN linear", "titan_linear_prprev_oof.csv.gz"),
                   ("UNI ABMIL", "mil_uni_v1_20x_fixedsplits_prprev_oof.csv.gz"),
                   ("CONCH ABMIL", "mil_conch_v1_20x_fixedsplits_prprev_oof.csv.gz")]
FROZEN = [("UNI mean", "frozen_1c_oof.csv.gz", "frozen_1e_pairs.csv"),
          ("UNI ABMIL", "frozen_abmil_1c_oof.csv.gz", "frozen_abmil_1e_pairs.csv")]


def logit(p):
    p = np.clip(p, 1e-4, 1 - 1e-4)
    return np.log(p / (1 - p))


def wilson(k, n, z=1.96):
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return c - h, c + h


def prop_stats(k, n):
    """k of n on the molecular side: CIs, effect size, Bayesian posterior and Bayes factors vs 50:50"""
    w = wilson(k, n)
    cp = binomtest(k, n).proportion_ci(method="exact")
    post = beta(k + 1, n - k + 1)
    log_bf10 = betaln(k + 1, n - k + 1) - n * np.log(0.5)
    p_gt = post.sf(0.5)
    return {"k_molecular": k, "n": n, "prop": k / n, "wilson_lo": w[0], "wilson_hi": w[1], "exact_lo": cp.low,
            "exact_hi": cp.high, "binom_p": binomtest(k, n).pvalue,
            "cohens_h": 2 * np.arcsin(np.sqrt(k / n)) - 2 * np.arcsin(np.sqrt(0.5)),
            "post_mean": post.mean(), "cred_lo": post.ppf(0.025), "cred_hi": post.ppf(0.975), "P_theta_gt_half": p_gt,
            "BF10": np.exp(log_bf10), "BFplus0": np.exp(log_bf10) * 2 * p_gt}


def mcnemar(b, c):
    """exact McNemar: b, c = the two discordant cells"""
    return binomtest(int(b), int(b + c)).pvalue if b + c else 1.0


def calib(y, p):
    bins = np.minimum((p * 10).astype(int), 9)
    ece = sum(abs(y[bins == i].mean() - p[bins == i].mean()) * (bins == i).mean() for i in range(10) if (bins == i).any())
    lr = LogisticRegression(C=1e6, max_iter=1000).fit(logit(p).reshape(-1, 1), y)
    return {"brier": np.mean((p - y) ** 2), "ECE": ece, "calib_intercept": lr.intercept_[0], "calib_slope": lr.coef_[0, 0]}


rng = np.random.default_rng(0)

# ---------------- 1a: Bayesian bootstrap of ΔAUC, calibration ----------------
rows = []
for cutoff, var, f in [("GMM", v, f) for v, f in VARIANTS] + [("prevmatch", v, f) for v, f in VARIANTS_PRPREV]:
    if not os.path.exists(R(f)):
        print("missing", f)
        continue
    d = pd.read_csv(R(f))
    d = d[d.exp == "1a"]
    P = d.groupby(["label", "patient"]).p.mean()  # seed-averaged out-of-fold P
    Y = d.groupby(["label", "patient"]).y.first()
    for marker, ihc, mol in (COMPARISONS_1A if cutoff == "GMM" else [("PR", "PR_IHC", "PGR_RNApos")]):
        if ihc not in P.index.get_level_values(0) or mol not in P.index.get_level_values(0):
            continue
        pts = P[ihc].index.intersection(P[mol].index)
        yi, ym = Y[ihc][pts].astype(int).values, Y[mol][pts].astype(int).values
        pi, pm = P[ihc][pts].values, P[mol][pts].values
        delta = roc_auc_score(ym, pm) - roc_auc_score(yi, pi)
        W = rng.dirichlet(np.ones(len(pts)), a.boot)
        bb = np.array([roc_auc_score(ym, pm, sample_weight=w) - roc_auc_score(yi, pi, sample_weight=w) for w in W])
        row = {"cutoff": cutoff, "variant": var, "marker": marker, "comparison": f"{mol} - {ihc}", "n": len(pts), "delta_AUC": delta,
               "bayes_median": np.median(bb), "cred_lo": np.percentile(bb, 2.5), "cred_hi": np.percentile(bb, 97.5),
               "P_delta_gt_0": (bb > 0).mean()}
        for tag, y, p in (("ihc", yi, pi), ("mol", ym, pm)):
            row.update({f"{k}_{tag}": v for k, v in calib(y, p).items()})
        rows.append(row)
        print(f"1a {cutoff:9s} {var:12s} {mol} - {ihc}: n={len(pts)} dAUC {delta:+.3f}, Bayesian bootstrap 95% CrI "
              f"[{row['cred_lo']:+.3f}, {row['cred_hi']:+.3f}], P(d>0) = {row['P_delta_gt_0']:.3f} | Brier IHC "
              f"{row['brier_ihc']:.3f} / mol {row['brier_mol']:.3f}, slope {row['calib_slope_ihc']:.2f} / "
              f"{row['calib_slope_mol']:.2f}")
pd.DataFrame(rows).to_csv(R("stats_ext_1a.csv"), index=False)

# ---------------- 1b: proportions, Bayes, paired McNemar ----------------
rows = []
c1b = pd.read_csv(R("stats_combined_1b_final.csv"))
for _, r in c1b.iterrows():
    for var in ["consensus"] + [c[4:] for c in c1b.columns if c.startswith("mol_")]:
        k = r["consensus_mol" if var == "consensus" else "mol_" + var]
        if pd.notna(k) and r.n_discordant > 0:
            rows.append({"source": "FFPE (DX), 5 variants", "cutoff": r.cutoff, "pair": r.pair, "model": var,
                         **prop_stats(int(k), int(r.n_discordant))})
for model, f1c, _ in FROZEN:
    if not os.path.exists(R(f1c)):
        continue
    d = pd.read_csv(R(f1c))
    d = d[d.comparison.str.startswith("1b")]
    for comp, g in d.groupby("comparison"):
        w = g.pivot_table(index="patient", columns="slide", values="p")
        y = g.groupby("patient").y.first().loc[w.index].astype(int).values
        side = {s: ((w[s].values >= 0.5).astype(int) == y) for s in ("DX", "TS")}
        for s in ("DX", "TS"):
            rows.append({"source": f"frozen-slide patients, {model}", "cutoff": "GMM/prevmatch", "pair": comp[3:],
                         "model": f"{model} {s}", **prop_stats(int(side[s].sum()), len(y))})
        b, c = int((side["DX"] & ~side["TS"]).sum()), int((~side["DX"] & side["TS"]).sum())
        rows.append({"source": f"frozen-slide patients, {model}", "pair": comp[3:], "model": f"{model} DX vs TS",
                     "n": len(y), "mcnemar_DXonly": b, "mcnemar_TSonly": c, "mcnemar_p": mcnemar(b, c)})
out = pd.DataFrame(rows)
out.to_csv(R("stats_ext_1b.csv"), index=False)
for r in out.itertuples():
    if pd.notna(getattr(r, "prop", np.nan)) and (r.model == "consensus" or "frozen" in r.source):
        print(f"1b {r.source[:28]:28s} {r.pair:40s} {r.model:16s} {r.k_molecular:.0f}/{r.n:.0f} = {r.prop:.2f} "
              f"Wilson [{r.wilson_lo:.2f}, {r.wilson_hi:.2f}] h={r.cohens_h:+.2f} P(theta>.5)={r.P_theta_gt_half:.3f} "
              f"BF10={r.BF10:.2f} BF+0={r.BFplus0:.2f}")
    elif pd.notna(getattr(r, "mcnemar_p", np.nan)):
        print(f"1b {r.source[:28]:28s} {r.pair:40s} McNemar DX-only {r.mcnemar_DXonly:.0f} vs TS-only "
              f"{r.mcnemar_TSonly:.0f}: p = {r.mcnemar_p:.3f}")

# ---------------- 1c: TOST equivalence and permutation test of the interaction ----------------
rows = []
for model, f1c, _ in FROZEN:
    if not os.path.exists(R(f1c)):
        continue
    d = pd.read_csv(R(f1c))
    d = d[~d.comparison.str.startswith("1b")]
    for comp, g in d.groupby("comparison", sort=False):
        mol, ihc = comp.split(" vs ")
        w = g.pivot_table(index="patient", columns=["slide", "label"], values="p")
        yy = g.pivot_table(index="patient", columns="label", values="y").loc[w.index]
        Y = {lab: yy[lab].astype(int).values for lab in (ihc, mol)}
        P = {(s, lab): w[s, lab].values for s in ("DX", "TS") for lab in (ihc, mol)}
        n = len(w)

        def inter(i, PP=P):
            A = {k: roc_auc_score(Y[k[1]][i], PP[k][i]) for k in PP}
            return (A["TS", mol] - A["TS", ihc]) - (A["DX", mol] - A["DX", ihc])

        est = inter(np.arange(n))
        bs = []
        for _ in range(a.boot):
            i = rng.integers(0, n, n)
            if min(len(set(Y[ihc][i])), len(set(Y[mol][i]))) == 2:
                bs.append(inter(i))
        bs = np.array(bs)
        lo90, hi90 = np.percentile(bs, [5, 95])
        p_tost = max((bs <= -a.margin).mean(), (bs >= a.margin).mean())
        perm = []
        for _ in range(a.boot):  # swap DX and TS predictions within a random half of patients
            sw = rng.random(n) < 0.5
            PP = {(s, lab): np.where(sw, P["TS" if s == "DX" else "DX", lab], P[s, lab]) for s, lab in P}
            perm.append(inter(np.arange(n), PP))
        p_perm = (np.sum(np.abs(perm) >= abs(est)) + 1) / (len(perm) + 1)
        rows.append({"model": model, "comparison": comp, "n": n, "interaction": est, "CI95_lo": np.percentile(bs, 2.5),
                     "CI95_hi": np.percentile(bs, 97.5), "CI90_lo": lo90, "CI90_hi": hi90, "margin": a.margin,
                     "TOST_p": p_tost, "equivalent": bool(lo90 > -a.margin and hi90 < a.margin),
                     "smallest_equiv_margin": max(abs(lo90), abs(hi90)), "perm_p": p_perm})
        print(f"1c {model:9s} {comp:28s} n={n} interaction {est:+.3f} 90% CI [{lo90:+.3f}, {hi90:+.3f}] | TOST "
              f"(+/-{a.margin}) p = {p_tost:.4f} {'equivalent' if rows[-1]['equivalent'] else 'not shown equivalent'}"
              f" | smallest margin {rows[-1]['smallest_equiv_margin']:.3f} | permutation p = {p_perm:.3f}")
pd.DataFrame(rows).to_csv(R("stats_ext_1c.csv"), index=False)

# ---------------- 1e: agreement beyond kappa ----------------


def agree(ct, cb):
    pa = (ct == cb).mean()
    pi = (ct.mean() + cb.mean()) / 2
    pe = 2 * pi * (1 - pi)
    return {"PABAK": 2 * pa - 1, "AC1": (pa - pe) / (1 - pe)}


rows = []
for model, _, f1e in FROZEN:
    if not os.path.exists(R(f1e)):
        continue
    d = pd.read_csv(R(f1e))
    for lab, g in d.groupby("label"):
        ct, cb = (g.P_TS.values >= 0.5).astype(int), (g.P_BS.values >= 0.5).astype(int)
        n, flips = len(g), int((ct != cb).sum())
        est = agree(ct, cb)
        bs = pd.DataFrame([agree(ct[i], cb[i]) for i in (rng.integers(0, n, n) for _ in range(a.boot))])
        b, c = int((ct > cb).sum()), int((ct < cb).sum())
        dl = logit(g.P_TS.values) - logit(g.P_BS.values)
        sd = dl.std(ddof=1)
        w = wilson(flips, n)
        row = {"model": model, "label": lab, "n_pairs": n, "flipped": flips, "pct_flipped": 100 * flips / n,
               "flipped_wilson_lo": 100 * w[0], "flipped_wilson_hi": 100 * w[1],
               "TSpos_BSneg": b, "TSneg_BSpos": c, "mcnemar_p": mcnemar(b, c),
               "BA_bias_logit": dl.mean(), "BA_bias_lo": dl.mean() - 1.96 * sd / np.sqrt(n),
               "BA_bias_hi": dl.mean() + 1.96 * sd / np.sqrt(n), "BA_LoA_lo": dl.mean() - 1.96 * sd,
               "BA_LoA_hi": dl.mean() + 1.96 * sd}
        for k, v in est.items():
            row[k], (row[k + "_lo"], row[k + "_hi"]) = v, np.percentile(bs[k], [2.5, 97.5])
        rows.append(row)
        print(f"1e {model:9s} {lab:11s} {n} pairs: flipped {row['pct_flipped']:.1f}% (Wilson {row['flipped_wilson_lo']:.1f}-"
              f"{row['flipped_wilson_hi']:.1f}) | PABAK {row['PABAK']:.3f} AC1 {row['AC1']:.3f} "
              f"[{row['AC1_lo']:.3f}, {row['AC1_hi']:.3f}] | McNemar TS+/BS- {b} vs TS-/BS+ {c}: p = {row['mcnemar_p']:.3f}"
              f" | Bland-Altman logit bias {row['BA_bias_logit']:+.2f}, LoA [{row['BA_LoA_lo']:+.2f}, "
              f"{row['BA_LoA_hi']:+.2f}]")
pd.DataFrame(rows).to_csv(R("stats_ext_1e.csv"), index=False)
print("saved", ", ".join(R(f"stats_ext_{x}.csv") for x in ("1a", "1b", "1c", "1e")))
