"""ReceptorCheck: research demo of an H&E receptor-status second reader (docs/DESIGN_SPEC.md §7.5).

Run:  streamlit run app/receptorcheck.py      (bundles: python scripts/app_build_cases.py)
Every number comes from out-of-fold, site-grouped predictions (the patient's hospital was never seen in training).
RESEARCH USE ONLY: not a medical device, not for clinical decisions.
"""
import json
import os
import sys
import time

import pandas as pd
import streamlit as st

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from receptorcheck_core import MARKERS, consistency, status  # noqa: E402

CASES = os.path.join(HERE, "cases")
COLOR = {"CONFIDENT": "#1b7f3b", "UNCERTAIN": "#b26b00", "OUT-OF-DISTRIBUTION": "#b00020"}
st.set_page_config(page_title="ReceptorCheck (research demo)", layout="wide")
st.markdown("<div style='background:#b00020;color:white;padding:6px 12px;border-radius:4px;font-weight:600'>"
            "RESEARCH USE ONLY · not a medical device · do not use for patient care</div>", unsafe_allow_html=True)
st.title("ReceptorCheck")
st.caption("ER / PR / HER2 from an H&E slide, with a reliability status and an IHC consistency check. "
           "Out-of-fold predictions on TCGA-BRCA; the patient's hospital was not in the training data.")

idx_path = os.path.join(CASES, "index.csv")
if not os.path.exists(idx_path):
    st.error("No case bundles found. Run `python scripts/app_build_cases.py` first.")
    st.stop()
IDX = pd.read_csv(idx_path)

# ---------------- sidebar ----------------
with st.sidebar:
    st.header("Case")
    filt = st.radio("Show", ["all", "ER IHC/RNA discordant", "ER-low", "frozen slide available", "with maps",
                             "flagged by the check"], index=0)
    sel = IDX
    if filt == "ER IHC/RNA discordant":
        sel = IDX[IDX.group.str.contains("discordant")]
    elif filt == "ER-low":
        sel = IDX[IDX.group.str.contains("ER-low")]
    elif filt == "frozen slide available":
        sel = IDX[IDX.frozen]
    elif filt == "with maps":
        sel = IDX[IDX.maps]
    elif filt == "flagged by the check":
        sel = IDX[IDX.ER_consistency == "FLAG FOR RETEST"]
    patient = st.selectbox(f"Patient ({len(sel)})", sel.patient.tolist())
    t0 = time.time()
    with open(os.path.join(CASES, f"{patient}.json")) as f:
        case = json.load(f)
    load_s = time.time() - t0
    slide = st.radio("Slide type", ["DX (FFPE)", "TS (frozen)"], index=0,
                     help="Frozen mode shows the FFPE-trained model applied to the frozen section (experiment G).")
    st.subheader("Lab IHC result (editable)")
    ihc_in = {}
    for m in MARKERS:
        if m in case["tasks"]:
            v = case["tasks"][m]["ihc"]
            ch = st.selectbox(f"{m} IHC", ["not entered", "positive", "negative"],
                              index=0 if v is None else (1 if v == 1 else 2), key=f"ihc_{patient}_{m}")
            ihc_in[m] = None if ch == "not entered" else int(ch == "positive")
    evaluation = st.toggle("Evaluation mode (show hidden molecular truth)", value=False)
    st.caption(f"bundle loaded in {load_s * 1000:.0f} ms · site {case['site']} · group: {case['group']} · "
               f"model {case['model']['version']}")

# ---------------- row 1: marker cards ----------------
ood = case["ood"]
if ood["flag"]:
    st.warning(f"OUT-OF-DISTRIBUTION: this slide's features are far from the training cohort "
               f"(Mahalanobis {ood['score']} > {ood['threshold']}). Do not use the estimates below.")
cols = st.columns(3)
for col, m in zip(cols, MARKERS):
    with col:
        if m not in case["tasks"]:
            st.info(f"{m}: no estimate")
            continue
        t = dict(case["tasks"][m])
        if slide.startswith("TS"):
            if case.get("frozen") and m in case["frozen"]:
                t["p_cal"] = case["frozen"][m]["p"]
                t["set"] = ["pos"] if t["p_cal"] >= 0.6 else ["neg"] if t["p_cal"] <= 0.4 else ["neg", "pos"]
            else:
                st.info(f"{m}: no frozen-slide estimate for this case")
                continue
        s = status(t, ood["flag"])
        call = "positive" if t["set"] == ["pos"] else "negative" if t["set"] == ["neg"] else "undetermined"
        st.markdown(f"<div style='border-left:6px solid {COLOR[s]};padding:4px 10px'>"
                    f"<h3 style='margin:0'>{m}: {call}</h3>"
                    f"<b>P({m}+) = {t['p_cal']:.2f}</b> · prediction set {{{', '.join(t['set']) or 'empty'}}}<br>"
                    f"<span style='color:{COLOR[s]};font-weight:600'>{s}</span></div>", unsafe_allow_html=True)
        if "expr" in t and slide.startswith("DX"):
            e = t["expr"]
            st.caption(f"predicted {e['gene']} log2 {e['pred_log2']} · {e['pct_vs_RNApos']}th percentile of "
                       f"measured {e['gene']} in RNA-positive tumours")
        if m == "HER2":
            st.caption("HER2 from H&E is weak (AUC about 0.6-0.7) and did not pass the acceptance checks; "
                       "HER2-low (0 vs 1+/2+) is not assessable from H&E.")

# ---------------- row 2: IHC consistency ----------------
st.subheader("IHC consistency check")
ccols = st.columns(3)
for col, m in zip(ccols, MARKERS):
    with col:
        if m not in case["tasks"]:
            continue
        t = case["tasks"][m]
        c = consistency(t, ihc_in.get(m), m)
        if c["state"] == "FLAG FOR RETEST":
            st.error(f"{m}: FLAG FOR RETEST · discordance score {c['score']} (threshold {t['tau']})")
        elif c["state"] == "CONSISTENT":
            st.success(f"{m}: consistent with IHC · score {c['score']}")
        elif c["state"] == "NOT VALIDATED":
            st.info(f"{m}: check not offered (failed the retrospective acceptance checks)")
        else:
            st.info(f"{m}: enter the IHC result to run the check")
conc = case.get("concepts", {})
nb = case.get("neighbours", [])
if conc:
    top2 = sorted(conc.items(), key=lambda kv: -kv[1])[:2]
    pam = pd.Series([n["PAM50"] for n in nb if n["PAM50"]]).value_counts()
    st.caption(f"Why: highest concepts vs the cohort: {top2[0][0]} (z {top2[0][1]:+.1f}), {top2[1][0]} "
               f"(z {top2[1][1]:+.1f}); most similar reference patients are mostly "
               f"{pam.index[0] if len(pam) else 'unknown'} ({pam.iloc[0] if len(pam) else 0} of {len(nb)}).")

# ---------------- row 3: evidence ----------------
tabs = st.tabs(["H&E", "Contribution map", "Attention map", "Concepts", "Similar cases"])
assets = case.get("assets", {})


def show(keys, caption):
    k = [x for x in keys if x in assets]
    if not k:
        st.info("Maps are available only for the map patients (filter 'with maps').")
        return
    for x in k:
        st.image(os.path.join(CASES, assets[x]), caption=f"{caption}: {x.replace('_', ' ')}", use_container_width=True)


with tabs[0]:
    show(["he"], "H&E thumbnail")
with tabs[1]:
    show(["contribution"], "Exact per-patch contribution of the linear probe (red = towards ER+, blue = towards ER-)")
with tabs[2]:
    show([k for k in assets if k.endswith("attention")], "MIL attention (percentile within the slide)")
with tabs[3]:
    if conc:
        s = pd.Series(conc).sort_values()
        st.bar_chart(s, horizontal=True)
        st.caption("CONCH zero-shot concept shares, z-scores against the TCGA cohort (0 = cohort average).")
    show([k for k in assets if k.startswith("conch_concept")], "CONCH concept map")
with tabs[4]:
    st.dataframe(pd.DataFrame(nb), hide_index=True, use_container_width=True)
    st.caption("Five most similar patients from the training hospitals (cosine similarity of UNI slide features).")

# ---------------- evaluation mode and footer ----------------
if evaluation:
    st.subheader("Hidden truth (evaluation mode)")
    st.json(case["truth_for_evaluation_only"])
st.divider()
st.caption("Scope: estimates within ±0.10 of 0.5 are within section-to-section variability (two frozen sections of one "
           "tumour flip 8-11 % of calls) and are shown as UNCERTAIN. HER2-low is not assessable. Prediction sets are "
           "split-conformal with 90 % target coverage. Out-of-distribution slides are rejected. "
           "Research use only.")
