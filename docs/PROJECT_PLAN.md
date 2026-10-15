# Project plan: "Which truth does the slide see?"

This plan is organised in **tiers**, not weeks. A tier is done when every box is ticked **and** its exit gate is met.
Work on a higher tier only when the tiers below it have passed their gates (optional Tier 6 last).
Tick boxes in this file as we go and log every change of scope in `comms/LOG.md`.

Legend: `[x]` done · `[~]` in progress (with name and date) · `[ ]` to do

**No fixed owners.** Either of us picks up whatever is needed next; see 0.6 for how we avoid doing the same task twice.

---

## 0. Guardrails: read before starting any task

### 0.1 The question (fixed; do not change without both agreeing in `comms/LOG.md`)

> When IHC, RNA and protein measurements of ER/PR/HER2 disagree, which one does an H&E foundation-model classifier
> follow, and does the answer change between the FFPE diagnostic slide and the frozen slide cut from the sequenced tissue?

**Novelty claim:** a new, stricter evaluation of existing models (same patients, same slides, same folds, different
reference standards). We do **not** propose a new architecture.

### 0.2 The five experiments (the whole project)

| ID | Question | Test | Data |
|---|---|---|---|
| 1a | Which label (IHC / RNA / RPPA) is more predictable from H&E? | Paired ΔAUC, same patients and folds; bootstrap CI | DX slides |
| 1b | On discordant patients, which side does a concordant-trained model take? | Binomial test vs 50:50 | DX slides |
| 1c | Does the answer change between FFPE (DX) and frozen (TS) slides? | 2 × 2 interaction (slide type × label source) | DX + TS slides |
| 1d | Is HER2 0 vs 1+ (and 0 vs low) visible in H&E at all? | AUC vs an ERBB2-mRNA-only comparator | DX slides |
| 1e | Do two frozen slides of one tumour agree? | ICC, Cohen's κ, % flipped calls | TS + BS pairs |

### 0.3 Frozen settings (use these everywhere; changing one means rerunning everything that depends on it)

| Setting | Value |
|---|---|
| Cohort | TCGA-BRCA; labels in `data/metadata/labels_patient.csv` (built by `scripts/build_labels.py`) |
| Labels | ER/PR IHC (clinical); RNA = GMM on log2 RSEM of ESR1/PGR/ERBB2; **HER2 strict** (IHC 3+ or FISH+ vs 0/1+ or 2+/FISH−), never `HER2_CLIN`; amplification = GISTIC 2; RPPA = GMM on protein |
| RNA cut-off robustness | Every RNA-label claim must hold under GMM, **prevalence-matched** and confident-only cut-offs |
| Splits | **`data/splits/folds.csv`** (site-grouped, 5 folds, seeds 0–4) via `common.site_folds()`. Never create folds any other way |
| Encoders | UNI v1 (main), CONCH v1 (replication), TITAN (fast baseline); all frozen |
| Classifiers | Logistic regression on mean-pooled / slide features (C = 0.01, class-balanced); ABMIL (settings in `scripts/train_mil.py`) |
| Repeats | ≥ 3 seeds for ABMIL, 5 seeds for linear models; report mean ± SD |
| Statistics | Paired bootstrap 95% CI (2,000) for ΔAUC; binomial test for 1b; report n for every number |
| Environment | scikit-learn 1.9.1 (pinned); heavy runs on Ada (`docs/ADA.md`) |

### 0.4 Out of scope (do not start these)

- New architectures or fusion models (the original B1/B2 pitches).
- Survival or outcome prediction.
- Radiology (MRI) or multimodal tri-fusion.
- Fine-tuning the encoder, except optional 1f in Tier 6.
- Hyper-parameter hunting to "win" AUC. Absolute AUC is not the point; the comparison between labels is.

### 0.5 Decision rules

- A finding goes in the report only if it is **consistent in direction across ≥ 2 encoders** and survives the cut-off schemes.
- If a result flips with the cut-off or the encoder, report it as "not robust" (the PR result is the example).
- Small subgroups (n < 20) are reported as exploratory leads, never as conclusions.
- Any change to labels, splits or settings: write it in `comms/LOG.md` first, then rerun affected tiers.

### 0.6 How we share the work (no fixed owners)

1. `git pull`, then read the **Status snapshot** below and the latest `comms/LOG.md` entries.
2. Pick the first open `[ ]` item in the **lowest tier that has not passed its gate** (Tier 6 only after Tiers 0–5).
3. Before starting, change it to `[~] in progress (Name, YYYY-MM-DD)`, commit and push, so the other person sees it.
4. When done: tick it `[x]`, add the result or file to the item, update the Status snapshot, add a `comms/LOG.md` entry, push.
5. Do not take an item marked `[~]` by the other person. If it looks stale (no commit for 2+ days), ask in `comms/LOG.md` first.
6. Ada jobs: each of us can submit under our own account (separate 40 CPU / 4 GPU quotas). Note the job IDs in the log.

---

## Tier 0: Foundations ✅

**Goal:** clean data, labels, splits and infrastructure that every later tier relies on.

- [x] Novelty audit; Option 1 chosen; B1/B2 dropped (`docs/NOVELTY_RESEARCH.md`)
- [x] Metadata pipeline: `fetch_metadata.py`, `build_labels.py` (1,097 patients, 4 label sources)
- [x] Dataset study and feasibility of all five experiments (`docs/DATASET.md`)
- [x] Repo, `CLAUDE.md`, `comms/LOG.md`, issue/PR templates
- [x] Ada access, venv, smoke test, storage plan (`docs/ADA.md`)
- [x] UNI v1 and CONCH v1 patch features for 1,055 patients (`data/features/` on Ada)
- [x] Fixed site-grouped splits `data/splits/folds.csv`; scikit-learn pinned
- [x] Both members have Ada key access and can run `jobs/check_env.sbatch` under their own account (Srivastava: job 19986)
- [ ] Both members have UNI/CONCH access on Hugging Face (IIIT e-mail)

**Exit gate:** labels, features and splits exist and are documented; both members can run jobs on Ada.

## Tier 1: Core baselines (1a, 1b, 1d on FFPE slides) ✅

**Goal:** answer 1a, 1b and 1d with three encoders and two classifier types.

- [x] Preliminary checks B0–B4 (label audit, label-only AUC, shuffled control, site-only, site identification)
- [x] TITAN linear: 1a, 1b, 1d, controls (5 seeds)
- [x] UNI and CONCH mean-pool: 1a, 1b, 1d, controls (3 seeds)
- [x] UNI and CONCH ABMIL: 1a, 1b, 1d, controls (3 seeds)
- [x] Sanity checks S1–S8 (`scripts/sanity_check.py`)
- [x] Baseline report (`docs/baseline_report.md`)

**Exit gate (met):** consistent ER result across 5 model variants; all sanity checks pass or are explained.

## Tier 2: Robustness and consistency ✅

**Goal:** make every Tier 1 claim reproducible and robust before interpreting it.

- [x] RNA cut-off sensitivity for ER and PR (GMM / prevalence-matched / confident-only)
- [x] PR prevalence-matched ABMIL re-check (UNI, CONCH): PR preference is a cut-off artefact
- [x] ER discordance breakdown (direction; ER-low vs higher % stained)
- [x] Rerun all linear baselines on the fixed splits, one environment (Ada job 19737): all conclusions unchanged
- [x] Update `docs/baseline_report.md` with the fixed-split numbers (Section 5b); ABMIL used the older folds
- [x] DeLong test (paired AUC) alongside the bootstrap for headline 1a comparisons (`scripts/stats_delong_1b.py`, `results/stats_delong_1a.csv`, job 19989; write-up in report pending)
- [x] Combined evidence across encoders for 1b (e.g. pooled binomial / sign consistency), clearly labelled (`results/stats_combined_1b.csv`: ER consensus 34:20, p = 0.076; HER2 42:18, p = 0.003; write-up pending)
- [x] Decide: rerun ABMIL on fixed splits? **No**: fixed-split linear numbers change no conclusion (revisit only for final tables)

**Exit gate:** every headline number comes from the fixed splits (or is explicitly flagged), and has n, mean ± SD and a CI or p-value.

## Tier 3: Interpretation of the FFPE results 🔄

**Goal:** explain why the slide follows RNA for ER (and amplification for HER2).

- [x] ABMIL attention heatmaps for 54 ER-discordant + 10 reference patients (Ada job 19722; `results/attention/`)
- [~] Overlay heatmaps on slide thumbnails for ~6 example patients (download those DX slides) to confirm attention sits on tumour (Srivastava, 2026-10-15)
- [~] Review the heatmaps: do ER-low cases look like ER-negative morphology? (Srivastava, 2026-10-15)
- [x] Attention statistics by group (entropy, fraction of patches holding 50% of attention; report Section 5c)
- [x] Describe discordant patients: PAM50, ER % stained, grade/histology (from the clinical table) (`scripts/describe_discordant.py`, report 5d: call follows intrinsic subtype; HER2-enriched IHC+/RNA− are the exception)
- [~] ER-low analysis written up as an exploratory lead (n = 11) (Srivastava, 2026-10-15)

**Exit gate:** a short, evidence-backed explanation section in the report, with 3–4 example heatmaps.

## Tier 4: Frozen slides (experiments 1c and 1e) ✅

**Goal:** test whether slide location (frozen next to the sequenced tissue vs FFPE block) changes which truth is tracked, and how reproducible single-slide predictions are.

- [x] UNI weights downloaded on Ada (needs Hugging Face approval) (Srivastava: loaded in job 21299)
- [x] Build GDC manifest of tumour TS/BS slides (1,517 slides, 430 GB) from `data/metadata/gdc_slides.csv` (`data/metadata/frozen_manifest.csv`: 1,516 slides, 373 pairs, 30 batches)
- [x] Time one TS slide through TRIDENT + UNI on an RTX 2080 Ti; log seconds per slide (job 21299, IDC source: download 30 s, segment 17 s, encode 15 s per slide; download overlaps encoding)
- [x] Extraction-equivalence check: our TRIDENT + UNI on ~50 DX slides vs TRACE features (cosine similarity, same predictions) (job 21300: 46 slides, cosine 1.0000 vs 0.33 for other patients, max |diff| ≤ 0.002, ER P(+) r = 1.000, 46/46 same call; 4 JPEG-2000 slides unreadable; `results/equivalence_check.csv`)
- [x] Batch download + encode all TS/BS slides (job array, delete raw slides after each batch); keep mean + capped bags (array 20080: 1,515 slides from the IDC copy; `~/project/data/features/frozen_uni_v1/` on Srivastava's Ada account)
- [x] Add % tumour cells per frozen slide (biospecimen table) as a covariate (columns in `data/metadata/frozen_manifest.csv`; 1,513 of 1,516 slides)
- [x] **1c:** 2 x 2 (DX vs TS x IHC vs RNA) with mean-pool (ABMIL optional), fixed splits, interaction test (job 22323, report 5e: ER interaction 0.000 [-0.024, +0.024], n = 947)
- [x] **1e:** TS vs BS agreement on 374 pairs: ICC of P(ER+), kappa of calls, % flipped (373 pairs: ICC 0.87-0.89, kappa 0.73/0.81, 8-11% flipped; report 5e)

**Exit gate:** 1c interaction estimate with CI and 1e agreement metrics, both with sanity checks like S1–S8.

## Tier 5: Statistics, consolidation and robustness of the whole story ⬜

**Goal:** one consistent set of final tables.

- [~] Final tables for 1a–1e from one environment and the fixed splits (Srivastava, 2026-10-15)
- [~] Multiple-testing note (how many tests, which are primary) (Srivastava, 2026-10-15)
- [~] Primary vs secondary claims written down (primary: ER 1a/1b; secondary: HER2, PR null, HER2-low, ER-low) (Srivastava, 2026-10-15)
- [~] Re-run `sanity_check.py` (extend it to the frozen-slide outputs) (Srivastava, 2026-10-15)

**Exit gate:** final results tables frozen; no number in the report lacks a script and file that produce it.

## Tier 6: Optional extensions (only after Tiers 0–5 pass their gates) ⬜

- [ ] External check on CPTAC-BRCA (134 patients: H&E, RNA, protein; no survival) for 1a/1b
- [ ] 1f: LoRA fine-tuning of UNI on IHC vs RNA labels (does adaptation change which truth is learned? ~30–60 GPU-h)
- [ ] Kamineni hospital data as an external test, if provided

## Tier 7: Deliverables ⬜

- [ ] Slides: fix PR bullet ("no reliable preference"), Ada CONCH numbers, remove/keep "attention maps" consistently
- [ ] Project guide / docs kept in sync with the final numbers
- [ ] Final report (problem, data, methods, results 1a–1e, robustness, interpretation, limitations)
- [ ] Final presentation and rehearsal; Q&A sheet (see the project guide)
- [ ] Repo clean-up: README up to date, every result reproducible from `scripts/` + `jobs/`

---

## Status snapshot (update when ticking boxes)

| Tier | Status | Blocking item |
|---|---|---|
| 0 Foundations | Done except Uniyal HF/CONCH item | Srivastava: Ada + UNI done, CONCH access pending |
| 1 Core baselines | Done | — |
| 2 Robustness | Done (14 Oct, report 5d) | — |
| 3 Interpretation | In progress | Thumbnail overlay + heatmap review (attention `.npz` on Uniyal's Ada account); ER-low write-up |
| 4 Frozen slides | Done (15 Oct, report 5e) | Sanity checks for frozen outputs moved to Tier 5 |
| 5 Consolidation | Not started (unblocked) | Final tables; multiple-testing note; extend `sanity_check.py` to 1c/1e |
| 6 Optional | Not started | Tiers 0–5 |
| 7 Deliverables | Partly (update deck exists) | Final numbers |

## Key files

| What | Where |
|---|---|
| Context and decisions | `docs/HANDOFF.md` |
| Baseline results and sanity checks | `docs/baseline_report.md`, `results/` |
| Dataset facts | `docs/DATASET.md` |
| Models and pipeline | `docs/PIPELINE_AND_MODELS.md` |
| Ada usage | `docs/ADA.md`, `jobs/` |
| Team log | `comms/LOG.md` |
