"""Tier 3: describe the ER IHC-vs-RNA discordant patients (PAM50, ER % stained, histology, stage, PR/HER2, ESR1 level)
and how the model variants call them. Descriptive only (n = 54): no claims beyond counts.

Needs results/stats_combined_1b_patients.csv from scripts/stats_delong_1b.py (per-patient calls of every variant).

Usage: python scripts/describe_discordant.py [--esr1-cut 10.17]
Output: results/er_discordant_cases.csv, results/er_discordant_summary.txt
"""
import argparse
import io
import os

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--labels", default="data/metadata/labels_patient.csv")
ap.add_argument("--clinical", default="data/metadata/clinical_patient_brca.txt")
ap.add_argument("--calls", default="results/stats_combined_1b_patients.csv")
ap.add_argument("--esr1-cut", type=float, default=10.17, help="GMM cut-off printed by build_labels.py")
ap.add_argument("--out", default="results")
a = ap.parse_args()

L = pd.read_csv(a.labels, index_col=0)
c = pd.read_csv(a.clinical, sep="\t", skiprows=[1, 2], dtype=str).set_index("bcr_patient_barcode")
calls = pd.read_csv(a.calls).set_index("patient")
variants = [v for v in calls.columns if v not in ("y_molecular", "n_variants_molecular", "consensus_p")]

d = L.loc[calls.index, ["ER_IHC", "ESR1_RNApos", "ESR1_log2", "PR_IHC", "HER2_STRICT", "PAM50", "ER_pct_IHC", "site"]]
d["direction"] = np.where(d.ER_IHC == 1, "IHC+/RNA-", "IHC-/RNA+")
d["ER_low"] = (d.ER_IHC == 1) & d.ER_pct_IHC.astype(str).str.contains("<10%")
d["ESR1_minus_cut"] = d.ESR1_log2 - a.esr1_cut
hist = c.histological_type.reindex(d.index).replace({"Infiltrating Ductal Carcinoma": "Ductal",
                                                     "Infiltrating Lobular Carcinoma": "Lobular"})
d["histology"] = np.where(hist.isin(["Ductal", "Lobular"]), hist, "Other")
d["stage"] = c.ajcc_pathologic_tumor_stage.reindex(d.index).str.extract(r"Stage (I{1,3}V?|IV)")[0]
d["age"] = pd.to_numeric(c.age_at_diagnosis.reindex(d.index), errors="coerce")
d["P_ERpos_consensus"] = calls.consensus_p
d["n_variants"] = len(variants)
d["n_variants_side_RNA"] = calls.n_variants_molecular
d["consensus_side"] = np.where((calls.consensus_p >= 0.5).astype(int) == d.ESR1_RNApos.astype(int), "RNA", "IHC")
d = d.sort_values(["direction", "consensus_side", "PAM50"])
d.to_csv(os.path.join(a.out, "er_discordant_cases.csv"))

buf = io.StringIO()


def out(*x):
    print(*x)
    print(*x, file=buf)


out(f"ER IHC-vs-RNA discordant patients: n = {len(d)}; model variants: {', '.join(variants)}")
out(f"Consensus = mean P(ER+) over the {len(variants)} variants >= 0.5; side = label the call agrees with.\n")
for title, col in [("PAM50", "PAM50"), ("ER % cells stained (IHC)", "ER_pct_IHC"), ("Histology", "histology"),
                   ("Stage", "stage")]:
    t = pd.crosstab([d.direction, d[col].fillna("NA")], d.consensus_side)
    out(f"--- {title} x consensus side (counts)\n{t.to_string()}\n")
# context: PAM50 of concordant patients
conc = L[(L.ER_IHC == L.ESR1_RNApos) & L.ER_IHC.notna()]
out("--- PAM50 of concordant patients, for reference (row %)")
out((pd.crosstab(conc.ER_IHC.map({1: "ER+ concordant", 0: "ER- concordant"}), conc.PAM50, normalize="index")
     * 100).round(0).to_string() + "\n")
g = d.groupby(["direction", "consensus_side"])
out("--- numeric summaries (median)")
out(g[["ESR1_minus_cut", "P_ERpos_consensus", "age"]].median().round(2).to_string())
out(g.size().rename("n").to_string() + "\n")
out(f"PR IHC positive: IHC+/RNA- {int(d[d.direction == 'IHC+/RNA-'].PR_IHC.sum())}/"
    f"{int((d.direction == 'IHC+/RNA-').sum())}, IHC-/RNA+ {int(d[d.direction == 'IHC-/RNA+'].PR_IHC.sum())}/"
    f"{int((d.direction == 'IHC-/RNA+').sum())}")
out(f"HER2 strict positive (of known): {int(d.HER2_STRICT.sum())}/{int(d.HER2_STRICT.notna().sum())}")
out(f"ER-low (IHC+, <10%): n = {int(d.ER_low.sum())}; consensus sides with RNA in "
    f"{int((d.ER_low & (d.consensus_side == 'RNA')).sum())}; PAM50 {d[d.ER_low].PAM50.value_counts().to_dict()}")
out(f"Unanimous (all {len(variants)} variants same side): RNA {int((d.n_variants_side_RNA == len(variants)).sum())}, "
    f"IHC {int((d.n_variants_side_RNA == 0).sum())}, split {int(((d.n_variants_side_RNA > 0) & (d.n_variants_side_RNA < len(variants))).sum())}")
with open(os.path.join(a.out, "er_discordant_summary.txt"), "w", encoding="utf-8") as f:
    f.write(buf.getvalue())
print("saved er_discordant_cases.csv, er_discordant_summary.txt")
