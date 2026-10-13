# Which Truth Does the Slide See?

**Auditing what pathology foundation models learn when IHC, mRNA and PAM50 disagree in breast cancer.**

Course project for *Topics in AI for Healthcare* (IIIT Hyderabad), built on the hospital project
**Kamineni – Multimodal Breast Cancer Diagnostic & Outcome Assistant**. Public data: TCGA-BRCA.

Team: Aman Uniyal · Aman Srivastava

> **New here (human or Claude Code)? Read in this order:**
> [`docs/HANDOFF.md`](docs/HANDOFF.md) → [`docs/PROJECT_PLAN.md`](docs/PROJECT_PLAN.md) → [`comms/LOG.md`](comms/LOG.md).
> Claude Code sessions should also read [`CLAUDE.md`](CLAUDE.md).

---

## The question

ER, PR and HER2 status decide breast-cancer treatment and are measured clinically by **IHC**. IHC disagrees
with mRNA-level status in a sizeable minority of patients, and for HER2 0 vs 1+ (HER2-low) pathologists agree
on only about a quarter of cases. H&E AI models are trained on whichever label a dataset provides, so nobody
knows **which biological truth they actually track**.

**Research question:** when IHC, RNA and PAM50 disagree, which one does an H&E foundation-model classifier
follow, and does the answer change between the FFPE diagnostic slide and the frozen slide cut from the
sequenced tissue?

Novelty type: a **new, stricter evaluation** of existing models (no new architecture).

| ID | Experiment | Status |
|---|---|---|
| 1a | IHC- vs RNA-defined label: which is more predictable from H&E? (same patients, same folds) | Preliminary done |
| 1b | Train on concordant cases; on IHC≠RNA cases, which side does the model take? | Preliminary done |
| 1c | FFPE-DX vs frozen slide × IHC vs RNA label (interaction) | Needs frozen-slide features |
| 1d | HER2 0 vs 1+: separable from H&E vs from ERBB2 mRNA? | First number only |
| 1e | Within-patient agreement: top (TS) vs bottom (BS) frozen slide | Needs frozen-slide features |

## Headline preliminary results (TITAN slide embeddings + logistic regression, site-grouped CV)

- **1a, ER:** RNA-defined ER is more predictable from H&E than IHC-defined ER: AUC 0.930 vs 0.882,
  Δ +0.048 (95% CI +0.023 to +0.074), same 989 patients.
- **1b, ER:** a model trained on concordant cases sides with **RNA in 39 of 54** discordant cases (p = 0.0015);
  36 of 45 after dropping borderline cases.
- **Leakage:** HER2 0 vs 1+ looks predictable with random folds (AUC 0.63) but falls to chance with
  site-grouped folds (0.43). The slide embedding identifies the hospital 82% of the time.

Full report with all encoders, ABMIL and sanity checks: [`docs/baseline_report.md`](docs/baseline_report.md) (running notes: [`docs/BASELINE_RESULTS.md`](docs/BASELINE_RESULTS.md)).

## Quick start (laptop, CPU, ~10 minutes)

```bash
pip install -r requirements.txt
python scripts/fetch_metadata.py       # optional: data/metadata/ is already committed
python scripts/build_labels.py         # -> data/metadata/labels_patient.csv + discordance audit
python scripts/fetch_titan.py          # -> data/features/titan_brca_20x.parquet (6.8 MB, not committed)
python scripts/baseline_titan.py       # baselines B2-B4 + preliminary 1b -> results/baseline_results.csv
python scripts/paired_ihc_vs_rna.py    # preliminary 1a + 1b robustness
```

Run everything from the repo root. Tested with Python 3.12 on Windows.

## Repository layout

```
README.md               this file
CLAUDE.md               rules + context for Claude Code sessions working in this repo
requirements.txt
docs/
  HANDOFF.md            full context to date: history, decisions, what is done, what is next
  PROJECT_PLAN.md       experiments, pipeline, timeline, owners, compute budget
  DATASET.md            dataset study: sources, label counts, discordance, slides, feasibility
  PIPELINE_AND_MODELS.md  foundation-model choice and step-by-step pipeline (Ada + laptop)
  baseline_report.md    consolidated baseline report: all encoders, ABMIL, sanity checks, interpretation
  BASELINE_RESULTS.md   running log of baseline numbers and how they were produced
  NOVELTY_RESEARCH.md   literature audit that led to this project (why the old pitches were dropped)
  archive/              original individual pitches (B1, B2) and the four-proposal research doc
comms/
  README.md             how the two of us (and our Claude Code sessions) communicate
  LOG.md                append-only message log, newest at the bottom
scripts/                runnable code (fetch metadata, build labels, baselines)
data/metadata/          small open-access TCGA metadata tables (committed, ~3 MB)
data/features/          large feature files (git-ignored)
results/                small result tables (committed)
slides/                 Beamer source of the project-update deck (compile on Overleaf)
.github/                issue and pull-request templates
```

## Data and licences

- TCGA-BRCA clinical, molecular and slide metadata: open access (GDC, cBioPortal).
- Precomputed features: [IUCompPath/TRACE](https://huggingface.co/datasets/IUCompPath/TRACE) (licence not stated).
  UNI, CONCH and TITAN weights are CC-BY-NC-ND. **Academic use only; never commit features or slides.**
