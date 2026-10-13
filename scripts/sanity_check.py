"""Sanity checks on the patch-feature results (mean-pool + ABMIL, UNI + CONCH) and paired bootstrap CIs.

Checks (each prints PASS/FAIL):
  S1  every (model, label, seed) OOF block: one prediction per patient, no NaN, probabilities in [0, 1]
  S2  AUC recomputed from the OOF predictions == AUC_mean in the summary
  S3  1a is paired: within a group and seed, all labels are predicted for exactly the same patients
  S4  the y stored in the OOF file equals the label table (no label mix-up), except the shuffled control
  S5  negative control (shuffled ER) AUC within 0.40-0.60 for every model
  S6  determinism: mean-pool rows of the full run == the earlier mean-only run
  S7  1b side counts recomputed from the stored discordant-case probabilities == the summary
  S8  n per experiment matches the expected cohort sizes
Also writes results/paired_ci_patch.csv: paired bootstrap 95% CI of AUC(molecular) - AUC(IHC), per model, seed 0.

Usage: python scripts/sanity_check.py [--results results] [--labels data/metadata/labels_patient.csv]
"""
import argparse
import os
import warnings

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.mixture import GaussianMixture

warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser()
ap.add_argument("--results", default="results")
ap.add_argument("--labels", default="data/metadata/labels_patient.csv")
a = ap.parse_args()
R = a.results
fails = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name} {detail}")
    if not ok:
        fails.append(name)


L = pd.read_csv(a.labels, index_col=0)


def gmm_pos(x):
    x = x.dropna()
    g = GaussianMixture(2, random_state=0).fit(x.values.reshape(-1, 1))
    return pd.Series(g.predict(x.values.reshape(-1, 1)) == np.argmax(g.means_.ravel()), index=x.index).astype(float)


L["ER_RPPApos"] = gmm_pos(L.ER_RPPA).reindex(L.index)
L["HER2_RPPApos"] = gmm_pos(L.HER2_RPPA).reindex(L.index)
GROUPS = {"ER": ["ER_IHC", "ESR1_RNApos", "ER_RPPApos"], "PR": ["PR_IHC", "PGR_RNApos"],
          "HER2": ["HER2_STRICT", "ERBB2_AMP", "HER2_RPPApos"]}
EXPECTED_N = {"ER": 801, "PR": 986, "HER2": 478}

for enc in ["uni_v1_20x", "conch_v1_20x"]:
    print(f"\n========== {enc} ==========")
    S = pd.read_csv(os.path.join(R, f"mil_{enc}_summary.csv"))
    O = pd.read_csv(os.path.join(R, f"mil_{enc}_oof.csv.gz"))
    print("summary rows:", len(S), "| oof rows:", len(O), "| models:", sorted(S.model.unique()))
    o1a = O[O.exp == "1a"]
    # S1
    blk = o1a.groupby(["model", "label", "seed"])
    dup = blk.patient.apply(lambda s: s.duplicated().any()).any()
    nan = o1a.p.isna().any()
    rng = ((o1a.p < 0) | (o1a.p > 1)).any()
    check(f"S1 {enc} OOF integrity", not (dup or nan or rng), f"(blocks={blk.ngroups}, dup={dup}, nan={nan}, out_of_range={rng})")
    # S2
    re = blk.apply(lambda d: roc_auc_score(d.y.astype(int), d.p)).groupby(level=[0, 1]).mean()
    s1a = S[S.exp == "1a"].set_index(["model", "label"]).AUC_mean
    diff = (re - s1a.reindex(re.index)).abs().max()
    check(f"S2 {enc} AUC recomputed from OOF matches summary", diff < 1e-6, f"(max |diff| = {diff:.2e})")
    # S3
    paired = True
    for model in o1a.model.unique():
        for gname, cols in GROUPS.items():
            for seed in o1a.seed.unique():
                sets = [frozenset(o1a[(o1a.model == model) & (o1a.label == c) & (o1a.seed == seed)].patient) for c in cols]
                paired &= all(s == sets[0] for s in sets)
    check(f"S3 {enc} 1a paired (same patients across labels)", paired)
    # S4
    m = o1a.merge(L, left_on="patient", right_index=True, how="left")
    mism = sum((m.label == c).sum() and ((m[m.label == c].y.values != m[m.label == c][c].values).sum()) for c in sum(GROUPS.values(), []))
    check(f"S4 {enc} stored y equals label table", mism == 0, f"(mismatches = {mism})")
    # S5
    ctl = S[S.label == "ER_IHC_shuffled"].set_index("model").AUC_mean
    check(f"S5 {enc} shuffled control at chance", ctl.between(0.40, 0.60).all(), f"({ctl.round(3).to_dict()})")
    # S6
    mo = pd.read_csv(os.path.join(R, f"mil_{enc}_summary_meanonly.csv"))
    a1 = S[(S.model == "mean") & S.AUC_mean.notna()].set_index("label").AUC_mean
    a2 = mo[mo.AUC_mean.notna()].set_index("label").AUC_mean
    d6 = (a1 - a2.reindex(a1.index)).abs().max()
    check(f"S6 {enc} mean-pool reproducible across runs", d6 < 1e-6, f"(max |diff| = {d6:.2e})")
    # S7
    o1b = O[O.exp == "1b"]
    ok7 = True
    for (model, label), d in o1b.groupby(["model", "label"]):
        k = int(((d.p >= 0.5).astype(int) == d.y.astype(int)).sum())
        row = S[(S.model == model) & (S.exp == "1b") & (S.label == label)]
        ok7 &= (len(row) == 1) and int(row.side_molecular.iloc[0]) == k and int(row.n.iloc[0]) == len(d)
    check(f"S7 {enc} 1b counts recomputed from stored probabilities", ok7)
    # S8
    ns = S[S.exp == "1a"].groupby("label").n.first()
    ok8 = all(ns[c] == EXPECTED_N[g] for g, cols in GROUPS.items() for c in cols)
    nb = S[S.exp == "1b"].groupby("label").n.first().to_dict()
    ok8 &= nb == {"ER_IHC vs ESR1_RNApos": 54, "ER_IHC vs ER_RPPApos": 63, "PR_IHC vs PGR_RNApos": 169,
                  "HER2_STRICT vs ERBB2_AMP": 60}
    check(f"S8 {enc} cohort sizes as expected", ok8, f"(1a n: {ns.to_dict()}; 1b n: {nb})")

    # prevalence-matched PR runs
    P = pd.read_csv(os.path.join(R, f"mil_{enc}_prprev_summary.csv"))
    PO = pd.read_csv(os.path.join(R, f"mil_{enc}_prprev_oof.csv.gz"))
    check(f"S8 {enc} PR prevalence-matched run present", (P.exp == "1b").any() and int(P[P.exp == "1b"].n.iloc[0]) == 133,
          f"(models {sorted(P.model.unique())}, 1b n={P[P.exp == '1b'].n.tolist()}, oof rows {len(PO)})")

# paired bootstrap CIs (seed 0) for every model, group, molecular label
rows = []
r = np.random.default_rng(0)
for enc in ["uni_v1_20x", "conch_v1_20x"]:
    O = pd.read_csv(os.path.join(R, f"mil_{enc}_oof.csv.gz"))
    PO = pd.read_csv(os.path.join(R, f"mil_{enc}_prprev_oof.csv.gz"))
    for src, oo, groups in [("main", O, GROUPS), ("PR prevmatch", PO, {"PR": GROUPS["PR"]})]:
        oo = oo[(oo.exp == "1a") & (oo.seed == 0)]
        for model in oo.model.unique():
            for gname, cols in groups.items():
                base = oo[(oo.model == model) & (oo.label == cols[0])].set_index("patient")
                for c in cols[1:]:
                    other = oo[(oo.model == model) & (oo.label == c)].set_index("patient").loc[base.index]
                    y0, p0, y1, p1 = base.y.values.astype(int), base.p.values, other.y.values.astype(int), other.p.values
                    bs = []
                    for _ in range(2000):
                        i = r.integers(0, len(y0), len(y0))
                        if len(set(y0[i])) < 2 or len(set(y1[i])) < 2:
                            continue
                        bs.append(roc_auc_score(y1[i], p1[i]) - roc_auc_score(y0[i], p0[i]))
                    lo, hi = np.percentile(bs, [2.5, 97.5])
                    rows.append({"encoder": enc, "run": src, "model": model, "marker": gname, "comparison": f"{c} - {cols[0]}",
                                 "n": len(y0), "delta_seed0": roc_auc_score(y1, p1) - roc_auc_score(y0, p0),
                                 "CI_lo": lo, "CI_hi": hi, "CI_excludes_0": (lo > 0) or (hi < 0)})
CI = pd.DataFrame(rows)
CI.to_csv(os.path.join(R, "paired_ci_patch.csv"), index=False)
print("\nPaired bootstrap 95% CIs (seed 0):")
print(CI.round(3).to_string(index=False))
print(f"\n{len(fails)} check(s) failed: {fails}" if fails else "\nALL SANITY CHECKS PASSED")
