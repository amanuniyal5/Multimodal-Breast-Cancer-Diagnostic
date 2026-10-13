"""Tier 2 statistics from saved out-of-fold predictions (no training): DeLong tests for 1a, combined evidence for 1b.

1a  For each model variant and each headline comparison (molecular label vs IHC label, same patients):
    - per seed: ΔAUC and DeLong p-value;
    - seed-averaged probabilities (each seed's prediction is out-of-fold, so their mean is too):
      ΔAUC, DeLong 95% CI and p-value, next to the paired bootstrap 95% CI (2,000 resamples) for comparison.
    The two AUCs use *different reference labels* on the same patients, so the DeLong covariance is computed from
    per-patient influence terms (DeLong 1988 structural components). With identical labels it reduces exactly to
    the classic DeLong test (checked in --selftest).
1b  The same discordant patients are scored by every model variant, so the variants are NOT independent and their
    counts cannot be summed into one binomial test. Combined evidence, at patient level (one call per patient):
    - consensus: mean P(positive) over variants, one call per patient, binomial test vs 50:50;
    - majority vote over variants, binomial test;
    - patients on which all variants agree, and on which side;
    - cluster bootstrap 95% CI (resampling patients) of the fraction of calls siding with the molecular label;
    - sign consistency: how many variants side with the molecular label in more than half of the patients.

Usage: python scripts/stats_delong_1b.py [--results results] [--boot 2000]
       python scripts/stats_delong_1b.py --selftest
Outputs: results/stats_delong_1a.csv, results/stats_combined_1b.csv, results/stats_combined_1b_patients.csv
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy.stats import binomtest, norm
from sklearn.metrics import roc_auc_score

ap = argparse.ArgumentParser()
ap.add_argument("--results", default="results")
ap.add_argument("--labels", default="data/metadata/labels_patient.csv")
ap.add_argument("--boot", type=int, default=2000)
ap.add_argument("--selftest", action="store_true")
a = ap.parse_args()

# (variant, file, model column value, folds). ABMIL predictions come from the pre-fixed-split folds (flagged).
VARIANTS = [
    ("TITAN linear", "titan_linear_fixedsplits_oof.csv.gz", "titan_linear", "fixed"),
    ("UNI mean", "mil_uni_v1_20x_fixedsplits_oof.csv.gz", "mean", "fixed"),
    ("CONCH mean", "mil_conch_v1_20x_fixedsplits_oof.csv.gz", "mean", "fixed"),
    ("UNI ABMIL", "mil_uni_v1_20x_oof.csv.gz", "abmil", "old"),
    ("CONCH ABMIL", "mil_conch_v1_20x_oof.csv.gz", "abmil", "old"),
]
# PR with the prevalence-matched PGR cut-off (the cut-off under which PR shows no preference)
VARIANTS_PRPREV = [
    ("TITAN linear", "titan_linear_prprev_oof.csv.gz", "titan_linear", "fixed"),
    ("UNI ABMIL", "mil_uni_v1_20x_prprev_oof.csv.gz", "abmil", "old"),
    ("CONCH ABMIL", "mil_conch_v1_20x_prprev_oof.csv.gz", "abmil", "old"),
]
COMPARISONS = [("ER", "ER_IHC", "ESR1_RNApos"), ("ER", "ER_IHC", "ER_RPPApos"), ("PR", "PR_IHC", "PGR_RNApos"),
               ("HER2", "HER2_STRICT", "ERBB2_AMP"), ("HER2", "HER2_STRICT", "HER2_RPPApos")]
PAIRS_1B = ["ER_IHC vs ESR1_RNApos", "ER_IHC vs ER_RPPApos", "PR_IHC vs PGR_RNApos", "HER2_STRICT vs ERBB2_AMP"]


def influence(y, p):
    """AUC and per-patient influence terms z with sum(z_k * z_l) = DeLong covariance (S10/m + S01/n form)."""
    y = np.asarray(y).astype(int)
    pos, neg = p[y == 1], p[y == 0]
    m, n = len(pos), len(neg)
    cmp = (pos[:, None] > neg[None, :]) + 0.5 * (pos[:, None] == neg[None, :])
    v10, v01 = cmp.mean(axis=1), cmp.mean(0)  # placement values
    auc = v10.mean()
    z = np.empty(len(y))
    z[y == 1] = (v10 - auc) / np.sqrt(m * (m - 1))
    z[y == 0] = (v01 - auc) / np.sqrt(n * (n - 1))
    return auc, z


def delong_paired(y_ref, p_ref, y_alt, p_alt):
    """ΔAUC (alt − ref) on the same patients, possibly with different labels; returns delta, se, p."""
    a_ref, z_ref = influence(y_ref, p_ref)
    a_alt, z_alt = influence(y_alt, p_alt)
    se = np.sqrt(np.sum((z_alt - z_ref) ** 2))
    d = a_alt - a_ref
    return a_ref, a_alt, d, se, 2 * norm.sf(abs(d) / se)


def classic_delong(y, p1, p2):
    """Textbook DeLong (1988) for two scores on one label, used only to check the generalised version."""
    y = np.asarray(y).astype(int)
    v10, v01, aucs = [], [], []
    for p in (p1, p2):
        pos, neg = p[y == 1], p[y == 0]
        c = (pos[:, None] > neg[None, :]) + 0.5 * (pos[:, None] == neg[None, :])
        v10.append(c.mean(axis=1)), v01.append(c.mean(0)), aucs.append(c.mean())
    s10, s01 = np.cov(np.vstack(v10)), np.cov(np.vstack(v01))
    m, n = (y == 1).sum(), (y == 0).sum()
    s = s10 / m + s01 / n
    var = s[0, 0] + s[1, 1] - 2 * s[0, 1]
    return aucs[1] - aucs[0], np.sqrt(var)


if a.selftest:
    r = np.random.default_rng(0)
    y = r.integers(0, 2, 400)
    p1, p2 = y + r.normal(0, 1.2, 400), y + r.normal(0, 1.0, 400)
    d1, se1 = classic_delong(y, p1, p2)
    _, _, d2, se2, _ = delong_paired(y, p1, y, p2)
    print(f"classic DeLong: delta {d1:.6f} se {se1:.6f} | generalised: delta {d2:.6f} se {se2:.6f}")
    assert abs(d1 - d2) < 1e-12 and abs(se1 - se2) < 1e-12
    # different labels: SE should agree with the bootstrap
    y2 = np.where(r.random(400) < 0.1, 1 - y, y)
    _, _, d3, se3, _ = delong_paired(y, p1, y2, p2)
    bs = []
    for _ in range(4000):
        i = r.integers(0, 400, 400)
        bs.append(roc_auc_score(y2[i], p2[i]) - roc_auc_score(y[i], p1[i]))
    print(f"different labels: generalised se {se3:.4f} | bootstrap sd {np.std(bs):.4f}")
    print("SELFTEST_OK")
    raise SystemExit


def load(fname, model):
    f = os.path.join(a.results, fname)
    if not os.path.exists(f):
        print(f"  (missing {f}: variant skipped)")
        return None
    d = pd.read_csv(f)
    return d[d.model == model]


# ---------------------------------------------------------------- 1a: DeLong next to the bootstrap
rows = []
rng = np.random.default_rng(0)
for cut, variants, comps in [("GMM", VARIANTS, COMPARISONS),
                             ("PR prevalence-matched", VARIANTS_PRPREV, [("PR", "PR_IHC", "PGR_RNApos")])]:
    for vname, fname, model, folds in variants:
        d = load(fname, model)
        if d is None:
            continue
        d = d[d.exp == "1a"]
        for marker, ref, alt in comps:
            w = {lab: d[d.label == lab].pivot(index="patient", columns="seed", values="p") for lab in (ref, alt)}
            if w[ref].empty or w[alt].empty:
                continue
            y = {lab: d[d.label == lab].groupby("patient").y.first() for lab in (ref, alt)}
            pts = w[ref].index.intersection(w[alt].index)
            assert len(pts) == len(w[ref]) == len(w[alt]), f"{vname} {ref}/{alt}: not paired"
            yr, ya = y[ref].loc[pts].values.astype(int), y[alt].loc[pts].values.astype(int)
            seeds = sorted(set(w[ref].columns) & set(w[alt].columns))
            per = [delong_paired(yr, w[ref].loc[pts, s].values, ya, w[alt].loc[pts, s].values) for s in seeds]
            pr, pa = w[ref].loc[pts].mean(axis=1).values, w[alt].loc[pts].mean(axis=1).values
            a_ref, a_alt, delta, se, p = delong_paired(yr, pr, ya, pa)
            bs = []
            for _ in range(a.boot):
                i = rng.integers(0, len(pts), len(pts))
                bs.append(roc_auc_score(ya[i], pa[i]) - roc_auc_score(yr[i], pr[i]))
            row = {"cutoff": cut, "variant": vname, "folds": folds, "marker": marker, "comparison": f"{alt} - {ref}",
                   "n": len(pts), "pos_ref": int(yr.sum()), "pos_alt": int(ya.sum()), "n_seeds": len(seeds),
                   "delta_seed_mean": np.mean([x[2] for x in per]), "delta_seed_sd": np.std([x[2] for x in per]),
                   "delong_p_seed_median": np.median([x[4] for x in per]),
                   "seeds_p_lt_0.05": int(sum(x[4] < 0.05 for x in per)),
                   "AUC_ref_avg": a_ref, "AUC_alt_avg": a_alt, "delta_avg": delta,
                   "delong_CI_lo": delta - 1.96 * se, "delong_CI_hi": delta + 1.96 * se, "delong_p": p,
                   "boot_CI_lo": np.percentile(bs, 2.5), "boot_CI_hi": np.percentile(bs, 97.5)}
            rows.append(row)
            print(f"1a [{cut}] {vname:12s} {alt:12s} - {ref:11s} n={len(pts):4d}  dAUC(avg) {delta:+.3f} "
                  f"DeLong CI [{row['delong_CI_lo']:+.3f}, {row['delong_CI_hi']:+.3f}] p={p:.2g} | "
                  f"boot CI [{row['boot_CI_lo']:+.3f}, {row['boot_CI_hi']:+.3f}] | "
                  f"per-seed dAUC {row['delta_seed_mean']:+.3f}+-{row['delta_seed_sd']:.3f}, "
                  f"p<0.05 in {row['seeds_p_lt_0.05']}/{len(seeds)} seeds", flush=True)
pd.DataFrame(rows).to_csv(os.path.join(a.results, "stats_delong_1a.csv"), index=False)

# ---------------------------------------------------------------- 1b: combined evidence over model variants
L = pd.read_csv(a.labels, index_col=0)
rows, pat_rows = [], []
for cut, variants, pairs in [("GMM", VARIANTS, PAIRS_1B), ("PR prevalence-matched", VARIANTS_PRPREV,
                                                           ["PR_IHC vs PGR_RNApos"])]:
    for pair in pairs:
        cols, ys = {}, {}
        for vname, fname, model, _ in variants:
            d = load(fname, model)
            if d is None:
                continue
            d = d[(d.exp == "1b") & (d.label == pair)].set_index("patient")
            if len(d):
                cols[vname], ys[vname] = d.p, d.y
        if not cols:
            continue
        P = pd.DataFrame(cols).dropna()
        y = pd.DataFrame(ys).loc[P.index].iloc[:, 0].astype(int)  # molecular label (same in every file)
        n, V = len(P), P.shape[1]
        side = (P.ge(0.5).astype(int).values == y.values[:, None])  # True = sides with the molecular label
        per_variant = side.sum(0)
        cons = int(((P.mean(axis=1) >= 0.5).astype(int) == y).sum())
        maj = int((side.sum(1) > V / 2).sum())
        ties = int((side.sum(1) == V / 2).sum())
        all_mol, all_ihc = int(side.all(1).sum()), int((~side).all(1).sum())
        frac = side.mean()
        bs = [side[rng.integers(0, n, n)].mean() for _ in range(a.boot)]
        row = {"cutoff": cut, "pair": pair, "n_discordant": n, "variants": "; ".join(P.columns), "n_variants": V,
               **{f"mol_{v}": int(k) for v, k in zip(P.columns, per_variant)},
               "variants_majority_molecular": int((per_variant > n / 2).sum()),
               "consensus_mol": cons, "consensus_p": binomtest(cons, n).pvalue,
               "majority_mol": maj, "majority_ties": ties,
               "majority_p": binomtest(maj, n - ties).pvalue if n - ties else np.nan,
               "all_agree_mol": all_mol, "all_agree_ihc": all_ihc,
               "frac_calls_mol": frac, "frac_CI_lo": np.percentile(bs, 2.5), "frac_CI_hi": np.percentile(bs, 97.5)}
        rows.append(row)
        print(f"1b [{cut}] {pair}: n={n}, {V} variants (" +
              ", ".join(f"{v} {k}" for v, k in zip(P.columns, per_variant)) +
              f"); consensus {cons}:{n - cons} p={row['consensus_p']:.2g}; majority {maj}:{n - ties - maj}"
              f" p={row['majority_p']:.2g}; all agree molecular {all_mol}, IHC {all_ihc}; "
              f"calls siding molecular {frac:.2f} [{row['frac_CI_lo']:.2f}, {row['frac_CI_hi']:.2f}]", flush=True)
        if pair == "ER_IHC vs ESR1_RNApos" and cut == "GMM":
            # subgroups: direction of discordance, ER-low (IHC+ with <10% cells stained)
            ihc = L.loc[P.index, "ER_IHC"].astype(int)
            low = L.loc[P.index, "ER_pct_IHC"].astype(str).str.contains("<10%") & (ihc == 1)
            for name, msk in [("IHC+/RNA-", ihc == 1), ("IHC-/RNA+", ihc == 0), ("ER-low (IHC+, <10%)", low)]:
                msk = msk.values
                c = int(((P.mean(axis=1) >= 0.5).astype(int) == y).values[msk].sum())
                k = int(msk.sum())
                rows.append({"cutoff": cut, "pair": f"{pair} | {name}", "n_discordant": k, "n_variants": V,
                             **{f"mol_{v}": int(s) for v, s in zip(P.columns, side[msk].sum(0))},
                             "consensus_mol": c, "consensus_p": binomtest(c, k).pvalue if k else np.nan,
                             "frac_calls_mol": side[msk].mean() if k else np.nan})
                print(f"   {name:20s} n={k:3d} consensus {c}:{k - c}", flush=True)
            pat_rows.append(P.assign(y_molecular=y, n_variants_molecular=side.sum(1),
                                     consensus_p=P.mean(axis=1)).reset_index().rename(columns={"index": "patient"}))
pd.DataFrame(rows).to_csv(os.path.join(a.results, "stats_combined_1b.csv"), index=False)
if pat_rows:
    pd.concat(pat_rows).to_csv(os.path.join(a.results, "stats_combined_1b_patients.csv"), index=False)
print("saved stats_delong_1a.csv, stats_combined_1b.csv, stats_combined_1b_patients.csv")
