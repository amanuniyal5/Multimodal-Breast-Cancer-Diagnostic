"""Sanity checks on the patch-feature results (mean-pool + ABMIL, UNI + CONCH) and paired bootstrap CIs.

Checks (each prints PASS/FAIL):
  S1  every (model, label, seed) OOF block: one prediction per patient, no NaN, probabilities in [0, 1]
  S2  AUC recomputed from the OOF predictions == AUC_mean in the summary
  S3  1a is paired: within a group and seed, all labels are predicted for exactly the same patients
  S4  the y stored in the OOF file equals the label table (no label mix-up), except the shuffled control
  S5  negative control (shuffled ER) AUC within 0.40-0.60 for every model
  S6  determinism: mean-pool rows of the full run == the earlier mean-only run (CONCH: the mean-only run was made on
      the laptop with scikit-learn 1.5, so its folds differ; known and superseded by the fixed-split rerun, report 5b)
  S7  1b side counts recomputed from the stored discordant-case probabilities == the summary
  S8  n per experiment matches the expected cohort sizes
  S9-S15  frozen slides (1c, 1e): site-grouped splits, OOF pairing, labels, recomputed AUCs/interaction, 1b counts,
          1e pairs and agreement, extraction equivalence and cohort sizes
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
# ================= frozen slides (experiments 1c, 1e; scripts/frozen_1c.py, frozen_1e.py) =================
print("\n========== frozen slides (1c, 1e) ==========")
F = pd.read_csv("data/splits/folds.csv", index_col="patient")
# S9: the fixed splits are site-grouped (no tissue source site in two folds of the same seed)
ok9 = all((F.groupby("site")[c].nunique() == 1).all() for c in F.columns if c.startswith("fold_s"))
check("S9 fixed splits are site-grouped", ok9, f"({F.site.nunique()} sites, {len(F)} patients)")

c1 = pd.read_csv(os.path.join(R, "frozen_1c.csv"))
o1 = pd.read_csv(os.path.join(R, "frozen_1c_oof.csv.gz"))
q = L.PGR_log2.quantile(1 - L.PR_IHC.mean())
L["PGR_prevmatch"] = (L.PGR_log2 >= q).astype(float).where(L.PGR_log2.notna())
a1c = o1[~o1.label.str.startswith("1b")]
# S10: 1c OOF integrity and pairing: per marker, the 4 blocks (DX/TS x IHC/molecular) cover the same patients
blk = a1c.groupby(["comparison", "slide", "label"])
ok10 = not (blk.patient.apply(lambda s: s.duplicated().any()).any() or a1c.p.isna().any()
            or ((a1c.p < 0) | (a1c.p > 1)).any())
for mk, d in a1c.groupby("comparison"):
    sets = {k: frozenset(g.patient) for k, g in d.groupby(["slide", "label"])}
    ok10 &= len(sets) == 4 and len(set(sets.values())) == 1
check("S10 1c OOF integrity and pairing (DX/TS x IHC/molecular, same patients)", ok10, f"(blocks={blk.ngroups})")
# S11: stored y equals the label table
m = a1c.merge(L, left_on="patient", right_index=True, how="left")
mism = sum(int((m[m.label == c].y.values != m[m.label == c][c].values).sum()) for c in a1c.label.unique())
check("S11 1c stored y equals label table", mism == 0, f"(mismatches = {mism})")
# S12: AUCs and interaction recomputed from the stored predictions == summary
ok12, dmax = True, 0.0
for r_ in c1[c1.comparison.str.contains(" vs ") & ~c1.comparison.str.startswith("1b")].itertuples():
    mol, ihc = r_.comparison.split(" vs ")
    A = {}
    for src in ("DX", "TS"):
        for lab in (ihc, mol):
            d = a1c[(a1c.comparison == r_.comparison) & (a1c.slide == src) & (a1c.label == lab)]
            A[src, lab] = roc_auc_score(d.y.astype(int), d.p)
    inter = (A["TS", mol] - A["TS", ihc]) - (A["DX", mol] - A["DX", ihc])
    dmax = max(dmax, abs(inter - r_.interaction), abs(A["DX", ihc] - r_.AUC_DX_ihc), abs(A["TS", mol] - r_.AUC_TS_mol))
    ok12 &= dmax < 1e-9
check("S12 1c AUCs and interaction recomputed from OOF match summary", ok12, f"(max |diff| = {dmax:.1e})")
# S13: 1b-on-frozen counts recomputed
ok13 = True
for r_ in c1[c1.comparison.str.startswith("1b")].itertuples():
    for src in ("DX", "TS"):
        d = o1[(o1.label == r_.comparison) & (o1.slide == src)]
        k = int(((d.p >= 0.5).astype(int) == d.y.astype(int)).sum())
        ok13 &= len(d) == r_.n and k == getattr(r_, f"side_mol_{src}")
check("S13 1c 1b counts recomputed from stored probabilities", ok13)

# S14: 1e pairs are real pairs and the summary is reproducible
e1 = pd.read_csv(os.path.join(R, "frozen_1e.csv"))
ep = pd.read_csv(os.path.join(R, "frozen_1e_pairs.csv"))
man = pd.read_csv("data/metadata/frozen_manifest.csv")
man["pair_id"] = man.pair_id.fillna("")
pm = man[man.pair_id != ""].groupby("pair_id").agg(n_ts=("section_type", lambda s: (s == "TS").sum()),
                                                    n_bs=("section_type", lambda s: (s == "BS").sum()),
                                                    patient=("patient", "first"))
ok14 = ep.pair_id.isin(pm.index).all() and (pm.loc[ep.pair_id, "patient"].values == ep.patient.values).all()
ok14 &= ((pm.n_ts >= 1) & (pm.n_bs >= 1)).all() and (pm.index.str[:12] == pm.patient).all()
for lab, d in ep.groupby("label"):
    row = e1[e1.label == lab].iloc[0]
    ct, cb = (d.P_TS >= 0.5).astype(int), (d.P_BS >= 0.5).astype(int)
    ok14 &= len(d) == row.n_pairs and abs(100 * (ct != cb).mean() - row.pct_flipped) < 1e-9
    ok14 &= (d.y.values == L.loc[d.patient, lab].values).all()
check("S14 1e pairs = one TS + one BS of the same patient and portion; flips and labels reproducible", ok14,
      f"({ep.pair_id.nunique()} pairs)")
# S15: extraction equivalence and cohort sizes
eq = pd.read_csv(os.path.join(R, "equivalence_check.csv"))
ok15 = (eq.cos_same_slide > 0.999).all() and ((eq.P_ER_trace >= 0.5) == (eq.P_ER_ours >= 0.5)).all()
n1c = c1[c1.comparison == "ESR1_RNApos vs ER_IHC"].n.iloc[0]
ok15 &= len(eq) >= 40 and n1c == 947 and e1[e1.label == "ER_IHC"].n_pairs.iloc[0] == 368 and man.slide_barcode.nunique() == 1515
check("S15 extraction equivalence (cosine > 0.999, same ER calls) and frozen cohort sizes", ok15,
      f"(equivalence n={len(eq)}, min cos={eq.cos_same_slide.min():.5f}; 1c ER n={n1c}; manifest slides={man.slide_barcode.nunique()})")

print(f"\n{len(fails)} check(s) failed: {fails}" if fails else "\nALL SANITY CHECKS PASSED")
