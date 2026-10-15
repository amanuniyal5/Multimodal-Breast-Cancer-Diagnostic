"""Tier 5 final tables: every number from the fixed site-grouped splits and one environment (scikit-learn 1.9.1, Ada).

Reads (all in results/):
  stats_delong_1a_final.csv, stats_combined_1b_final.csv  <- scripts/stats_delong_1b.py --abmil fixed --tag _final
  extended_titan.csv (TITAN 1d + shuffled control, fixed splits, job 19737)
  mil_{uni,conch}_v1_20x_fixedsplits_summary.csv (mean-pool), mil_{uni,conch}_v1_20x_fixedsplits_abmil_summary.csv
  frozen_1c.csv, frozen_1e.csv
Writes the tables into docs/FINAL_RESULTS.md (between the FINAL TABLES markers) and results/final_tables.csv.

Usage: python scripts/final_tables.py [--results results] [--doc docs/FINAL_RESULTS.md]
"""
import argparse
import os

import pandas as pd
import sys

sys.stdout.reconfigure(encoding="utf-8")  # tables contain Δ, κ, −

ap = argparse.ArgumentParser()
ap.add_argument("--results", default="results")
ap.add_argument("--doc", default="docs/FINAL_RESULTS.md")
ap.add_argument("--dry-run", action="store_true", help="allow old-fold ABMIL inputs (format test only)")
a = ap.parse_args()
R = lambda f: os.path.join(a.results, f)  # noqa: E731
VAR = ["TITAN linear", "UNI mean", "CONCH mean", "UNI ABMIL", "CONCH ABMIL"]
out, md = [], []


def fp(p):
    return "< 0.001" if p < 0.001 else f"{p:.3f}"


# ---------- T1: 1a, delta AUC (molecular - IHC), seed-averaged predictions, DeLong 95% CI
d = pd.read_csv(R("stats_delong_1a_final.csv"))
assert a.dry_run or (d.folds == "fixed").all(), "T1 must use fixed-split predictions only"
md.append("### T1. Experiment 1a: ΔAUC (molecular − IHC label), same patients, DeLong 95% CI\n")
md.append("| Comparison | n | " + " | ".join(VAR) + " |")
md.append("|---|---|" + "---|" * len(VAR))
order = [("GMM", "ESR1_RNApos - ER_IHC", "**ER: RNA − IHC**"), ("GMM", "ER_RPPApos - ER_IHC", "ER: RPPA − IHC"),
         ("GMM", "PGR_RNApos - PR_IHC", "PR: RNA − IHC (GMM cut-off)"),
         ("PR prevalence-matched", "PGR_RNApos - PR_IHC", "PR: RNA − IHC (prevalence-matched)"),
         ("GMM", "ERBB2_AMP - HER2_STRICT", "HER2: amplification − strict"),
         ("GMM", "HER2_RPPApos - HER2_STRICT", "HER2: RPPA − strict")]
for cut, comp, name in order:
    sub = d[(d.cutoff == cut) & (d.comparison == comp)].set_index("variant")
    cells = []
    for v in VAR:
        if v in sub.index:
            r = sub.loc[v]
            star = "**" if r.delong_CI_lo > 0 or r.delong_CI_hi < 0 else ""
            cells.append(f"{star}{r.delta_avg:+.3f}{star} ({r.delong_CI_lo:+.3f}, {r.delong_CI_hi:+.3f})")
            out.append({"table": "T1", "row": name, "variant": v, "n": int(r.n), "estimate": r.delta_avg,
                        "lo": r.delong_CI_lo, "hi": r.delong_CI_hi, "p": r.delong_p})
        else:
            cells.append("–")
    n = int(sub.n.iloc[0]) if len(sub) else 0
    md.append(f"| {name} | {n} | " + " | ".join(cells) + " |")
md.append("\nBold: CI excludes 0. Seed-averaged out-of-fold probabilities; TITAN 5 seeds, mean-pool 5, ABMIL 3.\n")

# ---------- T2: 1b, which side on discordant patients
b = pd.read_csv(R("stats_combined_1b_final.csv"))
md.append("### T2. Experiment 1b: discordant patients, molecular : IHC side (per variant and patient-level consensus)\n")
md.append("| Discordant pair | n | " + " | ".join(VAR) + " | Consensus (p) | Calls siding molecular (95% CI) |")
md.append("|---|---|" + "---|" * (len(VAR) + 2))
names = {"ER_IHC vs ESR1_RNApos": "**ER: IHC vs RNA**", "ER_IHC vs ER_RPPApos": "ER: IHC vs RPPA",
         "PR_IHC vs PGR_RNApos": "PR: IHC vs RNA", "HER2_STRICT vs ERBB2_AMP": "**HER2: strict vs amplification**"}
for r in b.itertuples():
    if "|" in r.pair:
        continue
    name = names.get(r.pair, r.pair) + (" (prevalence-matched)" if r.cutoff != "GMM" else "")
    cells = []
    for v in VAR:
        k = b.loc[r.Index].get(f"mol_{v}")
        cells.append("–" if pd.isna(k) else f"{int(k)} : {r.n_discordant - int(k)}")
    md.append(f"| {name} | {r.n_discordant} | " + " | ".join(cells) +
              f" | {r.consensus_mol} : {r.n_discordant - r.consensus_mol} ({fp(r.consensus_p)}) | "
              f"{100 * r.frac_calls_mol:.0f}% ({100 * r.frac_CI_lo:.0f}–{100 * r.frac_CI_hi:.0f}%) |")
    out.append({"table": "T2", "row": name, "variant": "consensus", "n": r.n_discordant, "estimate": r.consensus_mol,
                "p": r.consensus_p})
md.append("\nConsensus = one call per patient from the mean probability over the variants; binomial test vs 50 : 50.\n")

# ---------- T3: 1d HER2-low and the shuffled-label control
e = pd.read_csv(R("extended_titan.csv"))
rows3 = {"HER2low_1plus": "HER2 0 vs 1+", "HER2_0_vs_low": "HER2 0 vs low", "ER_IHC_shuffled": "Shuffled ER (control)"}
md.append("### T3. Experiment 1d and negative control: AUC (mean ± SD over seeds)\n")
md.append("| Task | n | " + " | ".join(VAR) + " | ERBB2 mRNA alone |")
md.append("|---|---|" + "---|" * (len(VAR) + 1))
summ = {}
for enc, tag in (("uni", "UNI"), ("conch", "CONCH")):
    s = pd.read_csv(R(f"mil_{enc}_v1_20x_fixedsplits_summary.csv"))
    summ[f"{tag} mean"] = s[s.model == "mean"].set_index("label")
    f = R(f"mil_{enc}_v1_20x_fixedsplits_abmil_summary.csv")
    if os.path.exists(f):
        s = pd.read_csv(f)
        summ[f"{tag} ABMIL"] = s[s.model == "abmil"].set_index("label")
for lab, name in rows3.items():
    if lab == "ER_IHC_shuffled":
        t = e[(e.exp == "E1 repeated CV") & (e.label == lab)].iloc[0]
        tit, n, mrna = (t.AUC_mean, t.AUC_sd), int(t.n), None
    else:
        t = e[(e.exp == "E5 1d") & (e.label == lab)].iloc[0]
        tit, n, mrna = (t.AUC_HE_mean, t.AUC_HE_sd), int(t.n), t.AUC_ERBB2_mRNA
    cells = [f"{tit[0]:.3f} ± {tit[1]:.3f}"]
    for v in VAR[1:]:
        if v in summ and lab in summ[v].index:
            r = summ[v].loc[lab]
            cells.append(f"{r.AUC_mean:.3f} ± {r.AUC_sd:.3f}")
            out.append({"table": "T3", "row": name, "variant": v, "n": int(r.n), "estimate": r.AUC_mean})
        else:
            cells.append("–")
    md.append(f"| {name} | {n} | " + " | ".join(cells) + f" | {'–' if mrna is None else f'{mrna:.3f}'} |")
md.append("")

# ---------- T4: 1c
c = pd.read_csv(R("frozen_1c.csv"))
md.append("### T4. Experiment 1c: FFPE (DX) vs frozen (TS) slide, UNI mean-pool, same patients and folds\n")
md.append("| Label pair | n | DX ΔAUC (95% CI) | TS ΔAUC (95% CI) | Interaction TS − DX (95% CI) | p |")
md.append("|---|---|---|---|---|---|")
for r in c[~c.comparison.str.startswith("1b")].itertuples():
    md.append(f"| {r.marker}: {r.comparison} | {r.n} | {r.delta_DX:+.3f} ({r.delta_DX_lo:+.3f}, {r.delta_DX_hi:+.3f}) | "
              f"{r.delta_TS:+.3f} ({r.delta_TS_lo:+.3f}, {r.delta_TS_hi:+.3f}) | "
              f"{r.interaction:+.3f} ({r.interaction_lo:+.3f}, {r.interaction_hi:+.3f}) | {fp(r.interaction_p_boot)} |")
    out.append({"table": "T4", "row": r.comparison, "variant": "UNI mean", "n": r.n, "estimate": r.interaction,
                "lo": r.interaction_lo, "hi": r.interaction_hi, "p": r.interaction_p_boot})
md.append("\n| 1b on the same discordant patients | n | DX molecular : IHC (p) | TS molecular : IHC (p) |")
md.append("|---|---|---|---|")
for r in c[c.comparison.str.startswith("1b")].itertuples():
    md.append(f"| {r.comparison[3:]} | {r.n} | {int(r.side_mol_DX)} : {int(r.side_ihc_DX)} ({fp(r.p_DX)}) | "
              f"{int(r.side_mol_TS)} : {int(r.side_ihc_TS)} ({fp(r.p_TS)}) |")
md.append("")

# ---------- T5: 1e
t = pd.read_csv(R("frozen_1e.csv"))
md.append("### T5. Experiment 1e: agreement of two frozen sections (TS vs BS) of one tumour, ER model\n")
md.append("| ER model trained on | Pairs | ICC(2,1) logit P (95% CI) | κ (95% CI) | Calls flipped (95% CI) | AUC TS / BS |")
md.append("|---|---|---|---|---|---|")
for r in t[t.label.isin(["ER_IHC", "ESR1_RNApos"])].itertuples():
    md.append(f"| {'IHC' if r.label == 'ER_IHC' else 'RNA'} label | {int(r.n_pairs)} | {r.ICC_logit:.3f} "
              f"({r.ICC_logit_lo:.3f}, {r.ICC_logit_hi:.3f}) | {r.kappa:.3f} ({r.kappa_lo:.3f}, {r.kappa_hi:.3f}) | "
              f"{r.pct_flipped:.1f}% ({r.pct_flipped_lo:.1f}, {r.pct_flipped_hi:.1f}) | {r.AUC_TS:.3f} / {r.AUC_BS:.3f} |")
    out.append({"table": "T5", "row": r.label, "variant": "UNI mean", "n": int(r.n_pairs), "estimate": r.kappa,
                "lo": r.kappa_lo, "hi": r.kappa_hi})
md.append("")

pd.DataFrame(out).to_csv(R("final_tables.csv"), index=False)
doc = open(a.doc, encoding="utf-8").read()
START, END = "<!-- FINAL TABLES START -->", "<!-- FINAL TABLES END -->"
block = START + "\n\n*Generated by `scripts/final_tables.py`; do not edit by hand.*\n\n" + "\n".join(md) + "\n" + END
if START in doc:
    doc = doc[:doc.index(START)] + block + doc[doc.index(END) + len(END):]
else:
    doc = doc.rstrip() + "\n\n" + block + "\n"
open(a.doc, "w", encoding="utf-8", newline="\n").write(doc)
print("\n".join(md))
print(f"\nwrote {a.doc} and {R('final_tables.csv')}")
