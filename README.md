# Which Truth Does the Slide See?

**Auditing what pathology foundation models learn when IHC, mRNA and protein disagree in breast cancer.**

Course project for *Topics in AI for Healthcare* (IIIT Hyderabad), built on the hospital project
**Kamineni – Multimodal Breast Cancer Diagnostic & Outcome Assistant**. Public data: TCGA-BRCA.

Team: Aman Uniyal · Aman Srivastava

> **New here (human or Claude Code)? Read in this order:**
> [`docs/HANDOFF.md`](docs/HANDOFF.md) → [`docs/PROJECT_PLAN.md`](docs/PROJECT_PLAN.md) → [`comms/LOG.md`](comms/LOG.md).
> Claude Code sessions should also read [`CLAUDE.md`](CLAUDE.md).
> **Final results:** [`docs/FINAL_RESULTS.md`](docs/FINAL_RESULTS.md) · **Report:** [`docs/FINAL_REPORT.md`](docs/FINAL_REPORT.md) (Word: `docs/FINAL_REPORT.docx`).

---

## The question

ER, PR and HER2 status decide breast-cancer treatment and are measured clinically by **IHC**, which disagrees with
RNA- or DNA-level status in 5-17% of TCGA-BRCA tumours. H&E AI models are trained on whichever label a dataset provides,
so it is unclear **which biological truth they actually track**.

**Research question:** when IHC, RNA and protein disagree, which one does an H&E foundation-model classifier follow,
and does the answer change between the FFPE diagnostic slide and the frozen slide cut from the sequenced tissue?

Novelty type: a **new, stricter evaluation** of existing models (same patients, slides and folds, different reference
standards); no new architecture.

| ID | Experiment | Status |
|---|---|---|
| 1a | IHC- vs RNA-defined label: which is more predictable from H&E? | Done (5 models) |
| 1b | Train on concordant cases; on IHC≠RNA cases, which side does the model take? | Done |
| 1c | FFPE (DX) vs frozen (TS) slide × IHC vs RNA label (interaction) | Done (1,515 frozen slides) |
| 1d | HER2 0 vs 1+ / 0 vs low: visible in H&E? | Done |
| 1e | Two frozen sections (TS vs BS) of one tumour: do predictions agree? | Done (373 pairs) |

## Headline results (fixed site-grouped splits; 5 models: TITAN, UNI and CONCH with linear or ABMIL heads)

- **ER, 1a:** the RNA-defined ER label is more predictable from H&E than the IHC label in all 5 models
  (ΔAUC +0.02 to +0.06; DeLong CI excludes 0 in 4 of 5; primary model UNI mean-pool +0.036, Holm-adjusted p = 0.042).
- **ER, 1b:** on 54 IHC-vs-RNA discordant patients every model leans to RNA (consensus 34 : 20, p = 0.076; frozen slides
  39 : 13, p = 0.0004); calls follow the PAM50 subtype (Basal → ER−, Luminal A → ER+; HER2-enriched is the exception).
- **ER, 1c:** the RNA advantage is identical on the FFPE and the frozen slide (+0.032 each; interaction 0.000,
  95% CI −0.024 to +0.024; n = 947): not a tissue-sampling artefact.
- **1e:** two frozen sections of one tumour agree well (ICC 0.87-0.89) but ~1 in 10 ER calls flips.
- **HER2** tracks DNA amplification (consensus 42 : 18, p = 0.003); **PR** shows no reliable preference;
  **HER2-low** is at chance from H&E (0.42-0.54) while ERBB2 mRNA alone reaches 0.65-0.69.
- Site leakage is real: embeddings identify the hospital 82% of the time, and HER2 0 vs 1+ looks learnable with random
  folds (0.63) but not with site-grouped folds.

## Where every result comes from

| Result | File(s) | Script → job |
|---|---|---|
| Labels, cohort | `data/metadata/labels_patient.csv` | `scripts/fetch_metadata.py`, `scripts/build_labels.py` |
| Fixed site-grouped splits | `data/splits/folds.csv` | `scripts/make_splits.py` |
| TITAN baselines, robustness | `results/baseline_results.csv`, `results/extended_titan.csv`, `results/robustness_mean.csv` | `scripts/baseline_titan.py`, `extended_titan.py`, `robustness_mean.py` → `jobs/rerun_linear.sbatch` |
| DX patch features (UNI, CONCH) | Ada `data/features/` | `scripts/build_bags.py` → `jobs/uni_features.sbatch` |
| Mean-pool and ABMIL, 1a/1b/1d (fixed splits) | `results/mil_{uni,conch}_v1_20x_fixedsplits*_{summary,oof}` | `scripts/train_mil.py` → `jobs/train_mil.sbatch`, `jobs/rerun_linear.sbatch` |
| DeLong 1a, combined 1b | `results/stats_delong_1a_final.csv`, `results/stats_combined_1b_final.csv` | `scripts/titan_oof.py`, `scripts/stats_delong_1b.py --abmil fixed --tag _final` |
| Discordant-case description | `results/er_discordant_*` | `scripts/describe_discordant.py` |
| Frozen-slide features | Ada `data/features/frozen_uni_v1/` | `scripts/make_frozen_manifest.py`, `make_idc_objects.py`, `extract_frozen.py` → `jobs/extract_frozen.sbatch` |
| Extraction equivalence (ours vs TRACE) | `results/equivalence_check.csv` | `scripts/equivalence_check.py` → `jobs/equivalence_check.sbatch` |
| 1c, 1e | `results/frozen_1c*`, `results/frozen_1e*` | `scripts/frozen_1c.py`, `frozen_1e.py` → `jobs/frozen_analysis.sbatch` |
| Attention maps and overlays | `results/attention/`, `results/attention_v2/`, `results/attention_overlay/` | `scripts/attention_maps.py`, `attention_overlay.py` → `jobs/attention_maps.sbatch`, `attention_overlay.sbatch` |
| Final tables T1-T5, sanity checks S1-S15 | `docs/FINAL_RESULTS.md`, `results/final_tables.csv` | `scripts/final_tables.py`, `scripts/sanity_check.py` → `jobs/final_tables.sbatch` |

Ada maintenance jobs: `setup_env`, `trident_setup`, `fix_torch` (restore the pinned PyTorch), `check_env`,
`requeue_frozen.sh`, `list_hf_archive` (one-off). Ada usage and gotchas: [`docs/ADA.md`](docs/ADA.md), [`docs/HANDOFF.md`](docs/HANDOFF.md) section 5.

## Quick start (laptop, CPU, ~10 minutes)

```bash
pip install -r requirements.txt
python scripts/build_labels.py         # -> data/metadata/labels_patient.csv + discordance audit
python scripts/fetch_titan.py          # -> data/features/titan_brca_20x.parquet (6.8 MB, not committed)
python scripts/extended_titan.py       # TITAN 1a/1b/1d on the fixed splits -> results/extended_titan.csv
python scripts/sanity_check.py         # S1-S15 on the committed result files
python scripts/final_tables.py         # rebuild T1-T5 in docs/FINAL_RESULTS.md from the committed *_final files
```

Run everything from the repo root. Patch-feature models and frozen slides need a GPU cluster (see `docs/ADA.md`).

## Repository layout

```
README.md               this file
CLAUDE.md               rules + context for Claude Code sessions working in this repo
docs/
  FINAL_RESULTS.md      final tables (fixed splits), primary/secondary claims, multiple-testing note
  FINAL_REPORT.md/.docx final report
  HANDOFF.md            full context: history, decisions, current state, gotchas
  PROJECT_PLAN.md       tiered plan with checklists and guardrails
  baseline_report.md    detailed results with methods, sanity checks and interpretation (Sections 4-5g)
  DATASET.md            dataset study; PIPELINE_AND_MODELS.md, ADA.md, NOVELTY_RESEARCH.md, archive/
comms/LOG.md            append-only team log
scripts/                all analysis code;  jobs/  SLURM jobs for the Ada cluster
data/metadata/          small TCGA metadata tables, frozen-slide manifest, IDC slide listing (committed)
data/splits/folds.csv   fixed site-grouped folds used by every experiment
results/                result tables, predictions, attention figures (committed; features are not)
slides/                 Beamer decks (compile on Overleaf)
```

## Data and licences

- TCGA-BRCA clinical, molecular and slide metadata: open access (GDC, cBioPortal). Frozen slides: the NCI Imaging Data
  Commons copy of TCGA (DICOM, public cloud).
- Precomputed features: [IUCompPath/TRACE](https://huggingface.co/datasets/IUCompPath/TRACE) (licence not stated).
  UNI, CONCH and TITAN weights are CC-BY-NC-ND. **Academic use only; never commit features or slides.**
