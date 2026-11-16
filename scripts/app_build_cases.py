"""ReceptorCheck case bundles (docs/DESIGN_SPEC.md §7.2-7.4) from out-of-fold outputs only, plus X5d retrieval.

Per patient (seed-0 site-grouped folds; the patient's hospital is never in training):
  tasks      ER / PR / HER2: ensemble probability (UNI + H-optimus-0 + TITAN probes, X6), Platt-calibrated p_cal and the
             split-conformal set (alpha 0.10) from results/bench/x6/oof_ensemble.csv.gz; predicted ESR1 / PGR / ERBB2
             (mean of the UNI, H-optimus-0, TITAN ridge OOF predictions, X3) with its percentile among measured values
             of RNA-positive tumours; status (CONFIDENT / UNCERTAIN / OUT-OF-DISTRIBUTION); IHC consistency with tau
             chosen on the other folds' patients (precision >= 0.5 vs RNA-defined discordance, else no flags).
  ood        Mahalanobis distance of the UNI patient mean in the PCA-50 space of the training patients (standardised);
             flag above the training patients' 99th percentile.
  concepts   CONCH zero-shot concept shares (results/conch_zeroshot_scores.csv) as z-scores vs the cohort.
  neighbours X5d: 5 most similar training-fold patients (cosine on standardised UNI means) with their labels.
  frozen     DX-trained UNI probe applied to the patient's frozen TS slide (experiment G), if available.
  assets     X5b map panels (results/bench/x5/maps/panels) copied to app/cases/<patient>/.
Demo set: all ER IHC/RNA-discordant patients, all ER-low, the map patients and 150 random concordant patients (seed 0).
Also writes results/bench/x5/retrieval.csv (neighbour-vote accuracy) and app/cases/index.csv.
Usage: python scripts/app_build_cases.py [--all]
"""
import argparse
import glob
import os
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "app"))
from bench import registry, stats, tasks  # noqa: E402
from common import site_folds  # noqa: E402
from receptorcheck_core import MARKERS, consistency, status, write  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--out", default="app/cases")
ap.add_argument("--all", action="store_true", help="bundle every patient, not only the demo set")
ap.add_argument("--n-random", type=int, default=150)
a = ap.parse_args()
L = tasks.labels()
ENS = pd.read_csv("results/bench/x6/oof_ensemble.csv.gz")
REG = pd.read_csv("results/bench/x3/oof_reg.csv.gz")
REG = REG[REG.encoder.isin(["uni_v1", "hoptimus0", "titan"])].groupby(["gene", "patient"]).pred.mean()
CON = pd.read_csv("results/conch_zeroshot_scores.csv", index_col="patient").filter(like="c_")
CONZ = (CON - CON.mean()) / CON.std()
U = registry.load_means("uni_v1")
G = pd.read_csv("results/week3/g_oof.csv.gz") if os.path.exists("results/week3/g_oof.csv.gz") else None
sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()

# ---------------- OOD and neighbours (seed-0 folds over all patients with UNI means) ----------------
pts = U.index.intersection(L.index).to_numpy()
OOD, NB = {}, {}
ret = []
for tr, te in site_folds(pts, 0):
    sc = StandardScaler().fit(U.loc[pts[tr]].values)
    Xtr, Xte = sc.transform(U.loc[pts[tr]].values), sc.transform(U.loc[pts[te]].values)
    pca = PCA(50, random_state=0).fit(Xtr)
    ztr = pca.transform(Xtr) / np.sqrt(pca.explained_variance_)
    zte = pca.transform(Xte) / np.sqrt(pca.explained_variance_)
    dtr, dte = np.linalg.norm(ztr, axis=1), np.linalg.norm(zte, axis=1)
    thr = np.percentile(dtr, 99)
    for p, d in zip(pts[te], dte):
        OOD[p] = {"score": round(float(d), 2), "threshold": round(float(thr), 2), "flag": bool(d > thr)}
    ntr = Xtr / np.linalg.norm(Xtr, axis=1, keepdims=True)
    nte = Xte / np.linalg.norm(Xte, axis=1, keepdims=True)
    S = nte @ ntr.T
    top = np.argsort(-S, axis=1)[:, :5]
    for i, p in enumerate(pts[te]):
        nb = []
        for j in top[i]:
            q = pts[tr][j]
            nb.append({"patient": q, "sim": round(float(S[i, j]), 3), "ER_IHC": None if pd.isna(L.at[q, "ER_IHC"]) else int(L.at[q, "ER_IHC"]),
                       "ESR1_RNApos": None if pd.isna(L.at[q, "ESR1_RNApos"]) else int(L.at[q, "ESR1_RNApos"]),
                       "PAM50": None if pd.isna(L.at[q, "PAM50"]) else str(L.at[q, "PAM50"])})
        NB[p] = nb
        for lab in ("ER_IHC", "ESR1_RNApos"):
            if pd.notna(L.at[p, lab]):
                v = [n[lab] for n in nb if n[lab] is not None]
                ret.append({"label": lab, "patient": p, "y": int(L.at[p, lab]), "vote": np.mean(v) if v else np.nan})
        if pd.notna(L.at[p, "PAM50"]) and nb[0]["PAM50"] is not None:
            ret.append({"label": "PAM50 top-1", "patient": p, "y": 1, "vote": float(nb[0]["PAM50"] == L.at[p, "PAM50"])})
R = pd.DataFrame(ret)
rrows = []
for lab, g in R.groupby("label"):
    g = g.dropna()
    if lab == "PAM50 top-1":
        rrows.append({"label": lab, "n": len(g), "accuracy": g.vote.mean(),
                      "chance": (L.PAM50.value_counts(normalize=True) ** 2).sum()})
    else:
        acc = ((g.vote >= 0.5).astype(int) == g.y).mean()
        rrows.append({"label": lab, "n": len(g), "accuracy": acc, "auc_vote": stats.safe_auc(g.y, g.vote),
                      "chance": max(g.y.mean(), 1 - g.y.mean())})
os.makedirs("results/bench/x5", exist_ok=True)
pd.DataFrame(rrows).to_csv("results/bench/x5/retrieval.csv", index=False)
print("X5d retrieval (5 nearest training-fold patients):\n" + pd.DataFrame(rrows).round(3).to_string(index=False))

# ---------------- tasks, tau (cross-fitted over folds) ----------------
fold0 = {}
for k, (tr, te) in enumerate(site_folds(pts, 0)):
    for p in pts[te]:
        fold0[p] = k
TASK = {}
for m, (ihc, mol, gene) in MARKERS.items():
    e = ENS[ENS.task == ihc].set_index("patient")
    e = e[e.index.isin(fold0)]
    sets = e.set_standard_neg.astype(bool).to_numpy(), e.set_standard_pos.astype(bool).to_numpy()
    y = e.y.astype(int).to_numpy()
    score = np.where(y == 1, 1 - e.p_cal, e.p_cal)
    contra = np.where(y == 1, sets[0] & ~sets[1], sets[1] & ~sets[0])
    disc = (L.loc[e.index, mol] != L.loc[e.index, ihc]) & L.loc[e.index, mol].notna()
    f = e.index.map(fold0).to_numpy()
    taus = {}
    for k in np.unique(f):
        o = f != k
        cand = np.sort(np.unique(score[o & contra]))
        tau = 1.01
        for t in cand:
            sel = o & contra & (score >= t)
            if sel.sum() >= 3 and disc.to_numpy()[sel].mean() >= 0.5:
                tau = float(t)
                break
        taus[k] = tau
    TASK[m] = (e, taus)
    print(f"{m}: tau per fold {dict((k, round(v, 3)) for k, v in taus.items())}")

# measured expression of RNA-positive tumours, for percentiles
REF = {gene: np.sort(L.loc[L[mol] == 1, f"{gene}_log2"].dropna().to_numpy()) for _, (_, mol, gene) in MARKERS.items()}

# ---------------- demo set ----------------
er = L[["ER_IHC", "ESR1_RNApos"]].dropna()
disc = er[er.ER_IHC != er.ESR1_RNApos].index
erlow = L[(L.ER_pct_IHC == "<10%") & (L.ER_IHC == 1)].index
maps = pd.read_csv("data/metadata/maps_patients.csv").patient
conc = er[er.ER_IHC == er.ESR1_RNApos].index.intersection(pts)
rnd = pd.Index(np.random.default_rng(0).choice(conc, a.n_random, replace=False))
demo = pts if a.all else pd.Index(disc).union(erlow).union(maps).union(rnd).intersection(pts)
group = {p: "random concordant" for p in rnd}
group.update({p: "ER-low (IHC+ <10%)" for p in erlow})
group.update({p: "ER IHC/RNA discordant" for p in disc})
for p in maps:
    group[p] = group.get(p, "map reference") + " · maps"

# ---------------- bundles ----------------
if os.path.isdir(a.out):
    for f in glob.glob(os.path.join(a.out, "*.json")):
        os.remove(f)
rows = []
for p in demo:
    case = {"patient": p, "slide_type": "DX", "hospital_seen_in_training": False, "site": str(L.at[p, "site"]),
            "group": group.get(p, "other"),
            "model": {"encoder": "UNI + H-optimus-0 + TITAN (ensemble)", "head": "mean-pool linear probe, Platt",
                      "version": sha},
            "ood": OOD[p], "tasks": {}, "concepts": {}, "neighbours": NB[p], "assets": {}, "frozen": None,
            "truth_for_evaluation_only": {}}
    for m, (ihc, mol, gene) in MARKERS.items():
        e, taus = TASK[m]
        if p not in e.index:
            continue
        r = e.loc[p]
        st = (["neg"] if r.set_standard_neg else []) + (["pos"] if r.set_standard_pos else [])
        t = {"p": round(float(r.p_ens), 3), "p_cal": round(float(r.p_cal), 3), "set": st,
             "tau": round(taus[fold0[p]], 3)}
        if (f"{gene}_log2", p) in REG.index:
            v = float(REG.loc[(f"{gene}_log2", p)])
            t["expr"] = {"gene": gene, "pred_log2": round(v, 2),
                         "pct_vs_RNApos": int(round(100 * np.searchsorted(REF[gene], v) / len(REF[gene])))}
        t["status"] = status(t, OOD[p]["flag"])
        t["ihc"] = None if pd.isna(L.at[p, ihc]) else int(L.at[p, ihc])
        t["consistency"] = consistency(t, t["ihc"], m)
        if m == "HER2":
            t["her2_low"] = "not assessable from H&E (HER2 0 vs 1+ AUC 0.42-0.54 across models)"
        case["tasks"][m] = t
        case["truth_for_evaluation_only"][mol] = None if pd.isna(L.at[p, mol]) else int(L.at[p, mol])
    case["truth_for_evaluation_only"].update({"PAM50": None if pd.isna(L.at[p, "PAM50"]) else str(L.at[p, "PAM50"]),
                                              "ER_pct_IHC": str(L.at[p, "ER_pct_IHC"])})
    if p in CONZ.index:
        case["concepts"] = {c[2:]: round(float(v), 2) for c, v in CONZ.loc[p].items()}
    if G is not None:
        g = G[(G.patient == p) & (G.label == "ER_IHC") & (G.train == "DX") & (G.test == "TS")]
        if len(g):
            case["frozen"] = {"ER": {"p": round(float(g.p.iloc[0]), 3), "note": "DX-trained UNI probe on the frozen TS slide"}}
    panels = sorted(glob.glob(os.path.join("results/bench/x5/maps/panels", f"{p}_*.png")))
    if panels:
        dst = os.path.join(a.out, p)
        os.makedirs(dst, exist_ok=True)
        for f in panels:
            shutil.copy(f, dst)
            case["assets"][os.path.basename(f)[len(p) + 1:-4]] = f"{p}/{os.path.basename(f)}"
    write(case, a.out)
    rows.append({"patient": p, "group": case["group"], "site": case["site"],
                 "ER_status": case["tasks"].get("ER", {}).get("status"),
                 "ER_consistency": case["tasks"].get("ER", {}).get("consistency", {}).get("state"),
                 "ood": OOD[p]["flag"], "frozen": case["frozen"] is not None, "maps": bool(panels)})
I = pd.DataFrame(rows)
I.to_csv(os.path.join(a.out, "index.csv"), index=False)
print(f"{len(I)} bundles -> {a.out}; groups:\n{I.group.value_counts().to_string()}\nER status:\n"
      f"{I.ER_status.value_counts().to_string()}\nER consistency:\n{I.ER_consistency.value_counts().to_string()}")
print("APP_BUILD_OK")
