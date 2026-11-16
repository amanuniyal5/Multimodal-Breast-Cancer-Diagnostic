"""X3 secondary targets: log2(RSEM + 1) of the 50 PAM50 genes (Parker et al. 2009) for TCGA-BRCA primary tumours, from the
same cBioPortal study and RNA profile as fetch_metadata.py (brca_tcga_pan_can_atlas_2018, rna_seq_v2_mrna).

Output: data/metadata/rna_panel.csv (patient x gene, primary tumour sample "-01" only). Asserts >= 1,000 patients.
Usage: python scripts/fetch_rna_panel.py [--out data/metadata/rna_panel.csv]
"""
import argparse

import numpy as np
import pandas as pd
import requests

CBIO = "https://www.cbioportal.org/api"
STUDY = "brca_tcga_pan_can_atlas_2018"
PAM50 = ["ACTR3B", "ANLN", "BAG1", "BCL2", "BIRC5", "BLVRA", "CCNB1", "CCNE1", "CDC20", "CDC6", "CDH3", "CENPF", "CEP55",
         "CXXC5", "EGFR", "ERBB2", "ESR1", "EXO1", "FGFR4", "FOXA1", "FOXC1", "GPR160", "GRB7", "KIF2C", "KRT14",
         "KRT17", "KRT5", "MAPT", "MDM2", "MELK", "MIA", "MKI67", "MLPH", "MMP11", "MYBL2", "MYC", "NAT1", "NDC80",
         "NUF2", "ORC6", "PGR", "PHGDH", "PTTG1", "RRM2", "SFRP1", "SLC39A6", "TMEM45B", "TYMS", "UBE2C", "UBE2T"]

ap = argparse.ArgumentParser()
ap.add_argument("--out", default="data/metadata/rna_panel.csv")
a = ap.parse_args()
assert len(PAM50) == 50
g = requests.post(f"{CBIO}/genes/fetch", params={"geneIdType": "HUGO_GENE_SYMBOL"}, json=PAM50, timeout=120).json()
ids = {x["entrezGeneId"]: x["hugoGeneSymbol"] for x in g}
missing = set(PAM50) - set(ids.values())
assert not missing, f"genes not found: {missing}"
r = requests.post(f"{CBIO}/molecular-profiles/{STUDY}_rna_seq_v2_mrna/molecular-data/fetch",
                  params={"projection": "SUMMARY"}, json={"entrezGeneIds": list(ids), "sampleListId": f"{STUDY}_all"},
                  timeout=300)
d = pd.DataFrame(r.json())
d = d[d.sampleId.str.endswith("-01")]
d["gene"] = d.entrezGeneId.map(ids)
d["log2"] = np.log2(d.value.astype(float).clip(lower=0) + 1)
P = d.pivot_table(index="patientId", columns="gene", values="log2", aggfunc="mean")[PAM50]
P.index.name = "patient"
assert len(P) >= 1000, f"only {len(P)} patients: truncated download?"
P.round(4).to_csv(a.out)
print(f"{P.shape[0]} patients x {P.shape[1]} genes -> {a.out}; missing values {int(P.isna().sum().sum())}")
