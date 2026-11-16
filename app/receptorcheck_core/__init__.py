"""ReceptorCheck core logic (docs/DESIGN_SPEC.md §7.3), shared by the case builder (scripts/app_build_cases.py), the
retrospective evaluation (scripts/app_evaluate.py) and the Streamlit UI (app/receptorcheck.py).

Research use only: every prediction in the bundles is out-of-fold and site-grouped (the patient's hospital was never
seen in training).
"""
import json
import os

MARKERS = {"ER": ("ER_IHC", "ESR1_RNApos", "ESR1"), "PR": ("PR_IHC", "PGR_prevmatch", "PGR"),
           "HER2": ("HER2_STRICT", "ERBB2_AMP", "ERBB2")}
# markers that passed the retrospective acceptance checks (scripts/app_evaluate.py, results/app_eval/acceptance.csv):
# HER2 did not (coverage 0.94 vs 0.90 target, FLAG precision not above the base rate), so it never flags
VALIDATED = {"ER": True, "PR": True, "HER2": False}
UNCERTAIN_BAND = 0.10  # |p_cal - 0.5| below this = within the section-to-section variability (1e flips)


def status(task, ood_flag):
    """CONFIDENT / UNCERTAIN / OUT-OF-DISTRIBUTION for one marker of one case"""
    if ood_flag:
        return "OUT-OF-DISTRIBUTION"
    if len(task["set"]) != 1 or abs(task["p_cal"] - 0.5) < UNCERTAIN_BAND:
        return "UNCERTAIN"
    return "CONFIDENT"


def consistency(task, ihc, marker=None):
    """IHC consistency check: score = calibrated probability of the opposite of the entered IHC call; FLAG FOR RETEST
    when the conformal set is a single class contradicting IHC and score >= tau (tau chosen on other folds for
    precision >= 0.5 against RNA-defined discordance)."""
    if ihc is None:
        return {"state": "NO IHC ENTERED", "score": None}
    if marker is not None and not VALIDATED.get(marker, True):
        return {"state": "NOT VALIDATED", "score": None}
    ihc = int(ihc)
    score = 1 - task["p_cal"] if ihc == 1 else task["p_cal"]
    contra = task["set"] == (["neg"] if ihc == 1 else ["pos"])
    flag = contra and score >= task["tau"]
    return {"state": "FLAG FOR RETEST" if flag else "CONSISTENT", "score": round(float(score), 3),
            "contradicting_singleton": bool(contra)}


def write(case, folder):
    os.makedirs(folder, exist_ok=True)
    with open(os.path.join(folder, f"{case['patient']}.json"), "w") as f:
        json.dump(case, f, indent=1)


def read(patient, folder):
    with open(os.path.join(folder, f"{patient}.json")) as f:
        return json.load(f)
