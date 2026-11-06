"""Targets built with the project's frozen label rules (same definitions as train_mil.py / extended_titan.py)."""
import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture

BINARY = ["ER_IHC", "ESR1_RNApos", "ER_RPPApos", "PR_IHC", "PGR_prevmatch", "HER2_STRICT", "ERBB2_AMP",
          "HER2_0_vs_low", "HER2low_1plus"]
REGRESSION = ["ESR1_log2", "PGR_log2", "ERBB2_log2"]
# label pairs for the label-preference check (IHC/clinical label, molecular label)
PAIRS = [("ER_IHC", "ESR1_RNApos"), ("ER_IHC", "ER_RPPApos"), ("PR_IHC", "PGR_prevmatch"), ("HER2_STRICT", "ERBB2_AMP")]


def gmm_pos(x):
    x = x.dropna()
    g = GaussianMixture(2, random_state=0).fit(x.values.reshape(-1, 1))
    return pd.Series(g.predict(x.values.reshape(-1, 1)) == np.argmax(g.means_.ravel()), index=x.index).astype(float)


def labels(path="data/metadata/labels_patient.csv", clinical="data/metadata/clinical_patient_brca.txt"):
    L = pd.read_csv(path, index_col=0)
    c = pd.read_csv(clinical, sep="\t", skiprows=[1, 2], dtype=str).set_index("bcr_patient_barcode")
    q = L.PGR_log2.quantile(1 - L.PR_IHC.mean())  # prevalence-matched PGR cut-off (6.95)
    L["PGR_prevmatch"] = (L.PGR_log2 >= q).astype(float).where(L.PGR_log2.notna())
    L["ER_RPPApos"] = gmm_pos(L.ER_RPPA).reindex(L.index)
    L["HER2_RPPApos"] = gmm_pos(L.HER2_RPPA).reindex(L.index)
    L["HER2low_1plus"] = L.HER2_IHCscore.map({"0": 0, "1+": 1})
    fish = c.her2_fish_status.reindex(L.index)
    low = (L.HER2_IHCscore == "1+") | ((L.HER2_IHCscore == "2+") & (fish == "Negative"))
    L["HER2_0_vs_low"] = np.where(L.HER2_IHCscore == "0", 0.0, np.where(low, 1.0, np.nan))
    L["site"] = L.site.astype(str)
    return L
