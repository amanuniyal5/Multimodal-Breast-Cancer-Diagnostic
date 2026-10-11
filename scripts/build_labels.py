"""Build the patient-level label table used by every Option 1 experiment.

Input : data/metadata/{clinical_patient_brca.txt, molecular.csv, pam50.csv}
Output: data/metadata/labels_patient.csv  (+ prints the discordance audit)

Label definitions
  ER_IHC / PR_IHC    clinical IHC status (Positive=1, Negative=0, else NaN)
  HER2_STRICT        1 if IHC 3+ or FISH Positive; 0 if IHC 0/1+ (FISH not +) or IHC 2+ with FISH Negative
  HER2_CLIN          looser: her2_status_by_ihc Positive or FISH Positive (noisy, kept for reference only)
  *_RNApos           2-component Gaussian mixture on log2(RSEM+1) of ESR1 / PGR / ERBB2; high component = 1
  ERBB2_AMP          GISTIC ERBB2 == 2
  ER_RPPA, HER2_RPPA RPPA protein (continuous)
  PAM50              cBioPortal SUBTYPE
  site               tissue source site = characters 6-7 of the barcode (used for site-grouped CV)

Usage: python scripts/build_labels.py [--meta data/metadata]
"""
import argparse
import os

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.mixture import GaussianMixture

ap = argparse.ArgumentParser()
ap.add_argument("--meta", default="data/metadata")
a = ap.parse_args()


def P(f):
    return os.path.join(a.meta, f)


c = pd.read_csv(P("clinical_patient_brca.txt"), sep="\t", skiprows=[1, 2], dtype=str).set_index("bcr_patient_barcode")
assert len(c) > 1000, f"clinical table has {len(c)} rows - download was truncated, re-run fetch_metadata.py"
m = pd.read_csv(P("molecular.csv"))
m = m[m.sampleId.str.endswith("-01")]  # primary tumour only
mol = m.pivot_table(index="patientId", columns=["prof", "gene"], values="value", aggfunc="first")
pam = pd.read_csv(P("pam50.csv")).set_index("patientId").value.str.replace("BRCA_", "")


def gmm_label(x):
    x = x.dropna()
    v = np.log2(x + 1).values.reshape(-1, 1)
    g = GaussianMixture(2, random_state=0).fit(v)
    hi = np.argmax(g.means_.ravel())
    lab = pd.Series(g.predict(v) == hi, index=x.index)
    return lab, v[lab.values].min()


def ihc(s):
    return s.map({"Positive": 1, "Negative": 0})


df = pd.DataFrame(index=c.index)
df["ER_IHC"] = ihc(c.er_status_by_ihc)
df["PR_IHC"] = ihc(c.pr_status_by_ihc)
st, fi, sc = c.her2_status_by_ihc, c.her2_fish_status, c.her2_ihc_score
loose = pd.Series(np.nan, index=c.index)
loose[(st == "Positive") | (fi == "Positive")] = 1
loose[((st == "Negative") & (fi != "Positive")) | ((st == "Equivocal") & (fi == "Negative"))] = 0
df["HER2_CLIN"] = loose
strict = pd.Series(np.nan, index=c.index)
strict[(sc == "3+") | (fi == "Positive")] = 1
strict[((sc.isin(["0", "1+"])) & (fi != "Positive")) | ((sc == "2+") & (fi == "Negative"))] = 0
df["HER2_STRICT"] = strict
df["HER2_IHCscore"] = sc.where(sc.isin(["0", "1+", "2+", "3+"]))
for g in ["ESR1", "PGR", "ERBB2"]:
    lab, thr = gmm_label(mol[("rna_seq_v2_mrna", g)])
    df[g + "_RNApos"] = lab.astype(float)
    print(f"{g}: GMM cut-off log2 {thr:.2f}, positive {int(lab.sum())}/{len(lab)}")
df["ERBB2_AMP"] = (mol[("gistic", "ERBB2")] == 2).astype(float).where(mol[("gistic", "ERBB2")].notna())
for g in ["ESR1", "PGR", "ERBB2"]:
    df[g + "_log2"] = np.log2(mol[("rna_seq_v2_mrna", g)] + 1)
df["ER_RPPA"] = mol[("rppa", "ESR1")]
df["HER2_RPPA"] = mol[("rppa", "ERBB2")]
df["PAM50"] = pam
df["ER_pct_IHC"] = c.er_status_ihc_Percent_Positive
df["site"] = df.index.str[5:7]
df.to_csv(P("labels_patient.csv"))

print(f"\npatients: {len(df)}, sites: {df.site.nunique()}")
for x, y in [("ER_IHC", "ESR1_RNApos"), ("PR_IHC", "PGR_RNApos"), ("HER2_STRICT", "ERBB2_AMP")]:
    t = pd.crosstab(df[x], df[y])
    n = t.values.sum()
    d = n - np.trace(t.values)
    print(f"\n{x} vs {y}: n={n}, discordant={d} ({d / n:.1%})\n{t}")
print("\nER IHC+ & Basal:", int(((df.ER_IHC == 1) & (df.PAM50 == "Basal")).sum()),
      "| HER2 strict- & Her2-enriched:", int(((df.HER2_STRICT == 0) & (df.PAM50 == "Her2")).sum()))
for lab, score in [("ER_IHC", "ESR1_log2"), ("HER2_STRICT", "ERBB2_log2")]:
    d = df[[lab, score]].dropna()
    print(f"label-only AUC {score} -> {lab}: {roc_auc_score(d[lab], d[score]):.3f} (n={len(d)})")
