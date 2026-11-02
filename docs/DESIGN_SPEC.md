# Design and technical specification: AI extensions and the *ReceptorCheck* application

*Written 3 Nov 2026. Implements the proposal in [`AI_EXTENSIONS.md`](AI_EXTENSIONS.md) (reviewed and accepted as is).
It does not change the finished results or the Week 4 plan (extended statistics; experiments C, D, G); those feed into
the application as described in §8. Section numbers X1–X6 refer to the proposal.*

---

## 1. Goals and non-goals

**Goals**
1. Turn the project into an AI benchmark:
   - encoders (X1);
   - aggregation methods (X2);
   - training targets, including regression (X3);
   - improvement methods (X4);
   - interpretability (X5);
   - reliability (X6).
   All of these run on the same patients, the same fixed site-grouped folds and the same statistics as the finished work.
2. Make every new model pass the same stress tests:
   - label preference (IHC vs molecular);
   - hospital signature;
   - FFPE → frozen transfer, where features exist.
3. Build *ReceptorCheck*, a working demo of an H&E receptor-status second reader that uses the best components from
   1–2.

**Non-goals**
- No new foundation model. No full fine-tuning; LoRA is a stretch only.
- No new confirmatory claims: P1–P3 and their Holm family stay as they are. Everything here is labelled secondary or
  engineering.
- No clinical deployment. The application is a retrospective research demo on TCGA.

---

## 2. System overview

```
             +----------------------------------------------------------------------------------+
  DATA       | labels_patient.csv · folds.csv · molecular.csv · pam50.csv · (new) rna_panel.csv  |
             +----------------------------------------------------------------------------------+
                                  |
  FEATURES   TRACE parquet (UNI, UNI2, GigaPath, CONCH, CTransPath, CHIEF, TITAN)   IDC slides → TRIDENT
             (streamed on /scratch → build_bags.py)                                 (ResNet-50, optional DINOv2)
                                  |                                                  |
                                  v                                                  v
             +--------------------- feature store  data/features/<enc>_20x_* --------------------+
             |  <enc>_mean.npz (mean, std)   <enc>_capK.npy (fp16 bags)   <enc>_index.csv          |
             |  <enc>_capK_coords.npy (new)  frozen_<enc>/<slide>.npz (TS/BS)                      |
             +----------------------------------------------------------------------------------+
                                  |
  HARNESS    scripts/bench/  (registry · tasks · folds · heads · metrics · stats · OOF writer)
                                  |
     +-------------+--------------+--------------+---------------+--------------+
     | X1 encoders | X2 MIL heads | X3 regression| X4 improve    | X6 reliability |
     +-------------+--------------+--------------+---------------+--------------+
                                  |                       X5 interpretability (reads OOF + bags + coords)
                                  v
  OUTPUTS    results/bench/<exp>/*_oof.csv.gz · *_summary.csv · figures/
                                  |
  APP        scripts/app_build_cases.py → app/cases/<patient>.json + PNGs → app/receptorcheck.py (Streamlit)
```

Design rules (unchanged from the project):
- Scripts run from the repo root, take paths as arguments and print their numbers.
- Folds come only from `common.site_folds()`.
- Heavy jobs run on Ada with scikit-learn 1.9.1.
- Nothing under `data/features/` is committed.

---

## 3. Data and feature store

### 3.1 Existing inputs (unchanged)

| File | Content | Used for |
|---|---|---|
| `data/metadata/labels_patient.csv` | ER_IHC, PR_IHC, HER2_STRICT, HER2_IHCscore, ESR1/PGR/ERBB2 RNApos and log2, ERBB2_AMP, ER/HER2 RPPA, PAM50, ER_pct_IHC, site | All tasks |
| `data/splits/folds.csv` | patient, site, fold_s0 … fold_s4 | Outer CV (and inner CV, §4.3) |
| `data/metadata/frozen_manifest.csv` | TS/BS slides, pairs | Frozen transfer, 1e reliability |
| `data/metadata/idc_slides.csv` | IDC objects for DX/TS/BS | ResNet-50 / DINOv2 extraction, overlays |

### 3.2 New inputs

| File | Built by | Content |
|---|---|---|
| `data/metadata/rna_panel.csv` | `scripts/fetch_rna_panel.py` (same cBioPortal API as `fetch_metadata.py`) | log2(RSEM+1) of the 50 PAM50 genes + *ESR1/PGR/ERBB2*, patient × gene, primary tumour only. Assert ≥ 1,000 patients |
| `data/metadata/dx_manifest.csv` | `scripts/make_dx_manifest.py` | One DX slide per row from `idc_slides.csv`. Same columns as `frozen_manifest.csv` where applicable; excludes the 48 JPEG-2000 DX slides |

### 3.3 Feature-store format (one convention for every encoder)

`data/features/<enc>_20x_mean.npz`
- Keys: `X` (patients × d float32, mean over patches); new key `S` (patients × d, per-patch std); `patients`.
- Existing files gain `S` when rebuilt. Readers must tolerate its absence.

`data/features/<enc>_20x_cap<K>.npy`
- Shape: patients × K × d, fp16, rows in `index.csv` order.
- Patches are randomly subsampled to K with seed 0 when a slide has more.

`data/features/<enc>_20x_index.csv`
- Columns: `patient, row, n_kept, n_total`.

`data/features/<enc>_20x_cap<K>_coords.npy` **(new)**
- Shape: patients × K × 2, int32, level-0 pixel coordinates of each kept patch, padded with −1.
- Needed for the feature maps (X5b). Written by `build_bags.py --keep-coords`; TRACE has a `coords` column.

`data/features/frozen_<enc>/<slide_barcode>.npz`
- As now: keys `mean`, `bag`, `coords`, `n_patches`.

**Encoder registry** `configs/encoders.json`: one entry per encoder.

```json
{"uni_v1":   {"dim": 1024, "patch_px": 256, "source": "trace", "trace_dir": "uni_v1/20x", "cap": 4096, "tcga_pretrained": false},
 "uni2_h":   {"dim": 1536, "patch_px": 256, "source": "trace", "trace_dir": "<verify>", "cap": 2048, "tcga_pretrained": false},
 "gigapath": {"dim": 1536, "patch_px": 256, "source": "trace", "trace_dir": "<verify>", "cap": 2048, "tcga_pretrained": false},
 "conch_v1": {"dim": 512,  "patch_px": 512, "source": "trace", "trace_dir": "conch_v1/20x", "cap": 4096, "tcga_pretrained": false},
 "ctranspath": {"dim": 768, "patch_px": 256, "source": "trace", "trace_dir": "<verify>", "cap": 2048, "tcga_pretrained": true},
 "chief":    {"dim": 768,  "patch_px": 256, "source": "trace", "trace_dir": "<verify>", "cap": 0, "tcga_pretrained": true},
 "resnet50_imagenet": {"dim": 1024, "patch_px": 256, "source": "trident", "trident_name": "resnet50", "cap": 2048, "tcga_pretrained": false},
 "dinov2_vitl_natural": {"dim": 1024, "patch_px": 256, "source": "custom", "timm": "vit_large_patch14_dinov2.lvd142m", "cap": 0, "tcga_pretrained": false},
 "titan":    {"dim": 768,  "source": "trace_slide", "cap": 0, "tcga_pretrained": false}}
```

- `<verify>`: the TRACE folder names are read from the dataset listing on day 1, before any download.
- `cap = 0` means means only, no bags.
- `tcga_pretrained = true` marks an encoder as a leakage contrast. It is shown in grey in every table and never chosen
  as "best".

### 3.4 Storage budget (home quota 50 GB, shared)

Bag size = 1,060 × K × d × 2 bytes.

| Item | Size | Decision |
|---|---|---|
| Existing UNI cap4096 / CONCH cap4096 / frozen UNI | ~8.9 / ~4.4 / ~13 GB | keep |
| UNI2-h, GigaPath cap2048 | ~6.7 GB each | keep only if X2 uses them; else means only |
| ResNet-50 cap2048 | ~4.4 GB | keep (needed for X5b maps comparison) |
| CTransPath cap2048 | ~3.3 GB | optional |
| All means (+std) | < 50 MB total | keep |
| Raw TRACE parquet (47–93 GB per encoder) | – | **never in home**: streamed on node-local `/scratch` inside the job, deleted at job end |

**Gate.** Before each feature-building job, run `du -sh ~` and check the room the job needs. Stop and ask if
> 40 GB is used.

---

## 4. Shared evaluation harness (`scripts/bench/`)

A small package that every experiment imports, so the experiments differ only in the piece under test.

| Module | Responsibility | Key functions |
|---|---|---|
| `registry.py` | Read `configs/encoders.json`; load means / bags / coords | `load_means(enc) -> DataFrame`, `load_bags(enc) -> (np.memmap, index)`, `load_coords(enc)` |
| `tasks.py` | Build every target from `labels_patient.csv` (+ `rna_panel.csv`) using the frozen rules (GMM cut-offs, prevalence-matched PGR, HER2_STRICT) | `task(name) -> Series`; `TASKS` list (§4.1) |
| `cohort.py` | Benchmark cohort = patients with features for **all** encoders in a comparison and a non-missing target | `common_patients(encs, task)` |
| `folds.py` | Outer folds via `common.site_folds`; inner folds (§4.3); label-efficiency subsamples (§5.1) | `outer(pts, seed)`, `inner(train_pts, seed)`, `subsample(train_pts, n, rep)` |
| `heads.py` | Linear probe, ridge, MLP, MIL heads (wraps `mil_models.py`) | `fit_predict(head, Xtr, ytr, Xte, cfg, seed)` |
| `metrics.py` | AUC, AUPRC, balanced accuracy, Brier, ECE (15 bins); Pearson, Spearman, R² | `score(task_type, y, p)` |
| `stats.py` | Paired bootstrap of metric differences (2,000 resamples over patients), DeLong for AUC pairs (reuses `stats_delong_1b.py`), Wilson CI | `paired_diff(y, p1, p2)`, `delong(y1, y2, p)` |
| `io.py` | OOF and summary writers with a fixed schema (§4.4) | `write_oof`, `write_summary` |

### 4.1 Tasks

| Task | Type | Definition | Note |
|---|---|---|---|
| `ER_IHC`, `ESR1_RNApos`, `ER_RPPApos` | binary | as in the finished work | label-preference axis |
| `PR_IHC`, `PGR_prevmatch` | binary | prevalence-matched PGR cut-off (6.95) | GMM version as sensitivity |
| `HER2_STRICT`, `ERBB2_AMP` | binary | as in the finished work | |
| `HER2_0_vs_low`, `HER2_0_vs_1plus` | binary | as in 1d | expected at chance; checks for leakage |
| `PAM50` | 5-class | LumA, LumB, Her2, Basal, Normal | macro-AUC (one-vs-rest), macro-F1 |
| `ESR1_log2`, `PGR_log2`, `ERBB2_log2` | regression | log2 RSEM | X3 |
| `PANEL_<gene>` | regression | 50 PAM50 genes | X3, optional |
| `SITE10` | 10-class | the 10 largest hospitals | representation check only (random stratified CV within those sites, as in the finished site-ID test) |

### 4.2 Training protocol (defaults; experiment sections list only what differs)

- **Outer CV:** 5 site-grouped folds. Linear and ridge heads use 5 seeds (fold_s0–s4); neural heads use 3 seeds
  (s0–s2). Predictions are averaged over seeds per patient, as in the finished work.
- **Linear probe:** `StandardScaler` + logistic regression, C = 0.01, class-balanced (the project default). A tuned-C
  variant (C ∈ {1e-3, 1e-2, 1e-1, 1}, inner-CV AUC) is reported as a sensitivity analysis only.
- **Neural MIL defaults:**

  | Setting | Value |
  |---|---|
  | Optimiser | AdamW, lr 2e-4, weight decay 1e-4 |
  | Epochs | 15 |
  | Training patches | 2,048 random patches per bag per step (all at test) |
  | Precision | fp16 autocast (2080 Ti) |
  | Loss | BCE with `pos_weight` = neg/pos |

  These match `train_mil.py`, so ABMIL results are directly comparable.
- **Hyper-parameters** (for new methods only, small grids) are chosen on inner folds of seed 0. They are then frozen
  for all seeds and outer folds and written to `results/bench/<exp>/hparams.json`.

### 4.3 Inner folds (no new random splits)

- For an outer fold at seed *s*, the training patients are re-split by their `fold_s{(s+1) mod 5}` column. That gives
  up to 5 inner folds; groups with no training patients are dropped.
- These are site-grouped by construction, deterministic and library-independent. This keeps the rule "never create
  folds another way".
- `scripts/sanity_check.py` gains S16: no inner-validation patient is in the outer test fold, and no hospital appears
  in both the inner train and inner validation sets.

### 4.4 Output schema

`results/bench/<exp>/<enc>__<head>__<task>_oof.csv.gz`
- Columns: `patient, seed, fold, y, pred` (probability, or the predicted value for regression; PAM50 gives 5
  probability columns).

`results/bench/<exp>/summary.csv`
- One row per (enc, head, task, variant).
- Columns: `exp, enc, head, task, variant, n, n_pos, metric, value, sd_seeds, ci_lo, ci_hi, ref_enc, ref_head, diff, diff_ci_lo, diff_ci_hi, diff_p`.
- The reference is UNI + mean-pool LR unless stated otherwise.

### 4.5 Statistics and reporting

- **Every number:** n, a 95% CI (patient-level bootstrap with 2,000 resamples; DeLong where two AUCs share patients),
  and SD over seeds.
- **Method vs reference:** paired bootstrap of the difference on the same patients. Comparisons are reported as a full
  table, never selected by significance.
- **Label-preference check per method:** ΔAUC(ESR1_RNApos − ER_IHC) with DeLong CI, plus 1b consensus counts on the
  54 discordant patients.
- **Status:** all X-experiments are secondary or engineering results. They are not added to the Holm family. A method
  "improves" only if the CI of its difference vs the reference excludes 0 **and** the direction holds for ≥ 2 encoders
  or both heads.

---

## 5. Experiment specifications

### 5.1 X1: encoder benchmark and "do we need pathology foundation models?"

| Item | Specification |
|---|---|
| Encoders | `resnet50_imagenet`, `uni_v1`, `uni2_h`, `gigapath`, `conch_v1`, `titan`; contrasts `ctranspath`, `chief`; optional `dinov2_vitl_natural` |
| Feature building | TRACE: `jobs/trace_features.sbatch ENC=<name>` (CPU, `/scratch`, streams parquet → `build_bags.py --keep-coords --cap <K>`). ResNet-50 / DINOv2: `jobs/extract_dx.sbatch` running `extract_frozen.py --manifest data/metadata/dx_manifest.csv --encoder resnet50` (new choices `resnet50`, `dinov2_natural`; TRIDENT `resnet50` = ImageNet ResNet-50 truncated after stage 3, 1,024-d, as in CLAM) |
| Equivalence check | ResNet-50 on the same 46 DX slides used before vs a second run (determinism); for TRACE encoders, 5 slides re-extracted with TRIDENT vs TRACE (cosine) where weights are available |
| Heads | Mean-pool LR (primary); mean + std LR; ABMIL for the top 2 non-contrast encoders by inner-CV AUC |
| Tasks | All of §4.1 except the regression panel |
| Label efficiency | Training subsets n ∈ {25, 50, 100, 200, 400, all}, stratified by label, 5 repeats per n per outer fold (`folds.subsample`, seed = 1000·s + rep), always tested on the full outer test fold. Mean-pool LR only (cheap). Output: AUC vs n with 95% band per encoder, and the **n at which each encoder reaches 0.85 ER AUC** (bootstrap CI) |
| Extra outputs | `SITE10` accuracy per encoder; ER ΔAUC (RNA − IHC) per encoder; HER2-low AUC per encoder |
| Script | `scripts/bench_encoders.py --encoders ... --heads ... --tasks ... --out results/bench/x1` |
| Compute | TRACE streaming: ~2–4 h per encoder (download-bound, CPU job). ResNet-50 DX extraction: ~1,055 slides at an estimated 60–120 s per slide → 18–35 GPU-hours, split into 3 parallel array batches (fits the 40-CPU cap). Heads: < 2 h CPU; ABMIL 2 × ~4 h GPU |
| Decision it feeds | Encoder for the application (§8): best ER/PR/HER2 AUC under site-grouped CV, tie-break on label efficiency then on lower `SITE10` accuracy |
| Expected result | ResNet-50 clearly below the pathology foundation models, with the gap largest at small n; UNI2/GigaPath ≈ UNI (±0.02); contrasts possibly inflated |

### 5.2 X2: aggregation methods beyond ABMIL

Implemented in `scripts/mil_models.py`: one file, plain PyTorch, a shared interface `forward(x) -> (logit, aux)`.
`aux` carries attention or instance scores for X5. Methods are re-implemented from the papers, kept small, with no
external MIL packages.

| Method | Architecture (defaults) | Extra loss |
|---|---|---|
| ABMIL (reference) | existing class from `train_mil.py` (h 256, att 128, dropout 0.25) | – |
| Mean + std, top-k (k = 64) pooling | Pooled vector → LR / 1-hidden-layer MLP (256) | – |
| CLAM-SB | Gated attention + instance classifier on top-8 / bottom-8 attended patches | Instance loss (smooth-SVM), bag weight 0.7 |
| TransMIL | Linear 1024→512, 2 Nyström self-attention layers (8 heads, 256 landmarks), PPEG positional conv, CLS token | – |
| DSMIL | Instance classifier → critical instance → attention by distance to it | 0.5 × (bag + max-instance) loss |
| ACMIL | 5 attention branches, stochastic top-10 attention masking (p 0.6), diversity loss | Branch-diversity loss |
| kNN / prototype | Cosine kNN on patient means (k ∈ {10, 25}) / class prototypes | – |

| Item | Specification |
|---|---|
| Encoders | UNI (all methods); CONCH and the X1 winner for the best 2 methods (replication) |
| Tasks | ER_IHC, ESR1_RNApos, PR_IHC, PGR_prevmatch, HER2_STRICT, ERBB2_AMP |
| Tuning | Inner folds, seed 0, ≤ 4 configs per method (lr {1e-4, 2e-4} × dropout {0.1, 0.25}); epochs fixed at 15 (20 for TransMIL) |
| Script / job | `scripts/bench_mil.py --method <m> --enc <e>`; `jobs/bench_mil.sbatch` (1 GPU, `-c 6`, ~2–4 h per method per encoder) |
| Memory | 2,048 patches × 1,024-d fp16 per step; TransMIL Nyström is linear in the number of patches (fits 11 GB) |
| Expected result | Gains ≤ 0.02 AUC for ER (a global-morphology task); possibly larger for HER2 (focal). A null result is a finding |

### 5.3 X3: regression on mRNA expression

| Item | Specification |
|---|---|
| Targets | `ESR1_log2`, `PGR_log2`, `ERBB2_log2` (primary); `PANEL_*` 50 genes (secondary) |
| Heads | Ridge on means (α ∈ {1, 10, 100, 1000} by inner CV); ABMIL with a linear regression head and MSE loss on the z-scored target (training-fold mean / SD) |
| Metrics | Pearson r and Spearman ρ with bootstrap CI, R² (out-of-fold) |
| Regression → classification | Use the predicted value as a score: AUC against ER_IHC and ESR1_RNApos (and the PR, HER2 equivalents). Compare with the binary classifier trained on each label (paired DeLong) |
| Questions answered | (a) Does training on the continuous signal beat training on a binary label? (b) Does the predicted *PGR* remove the PR cut-off artefact (AUC vs PR_IHC independent of the cut-off)? (c) Which PAM50 genes are predictable from H&E (ranked r with CI)? |
| Script | `scripts/bench_regression.py`; panel genes via `fetch_rna_panel.py` |
| Compute | Ridge: minutes, CPU. ABMIL regression: 3 targets × ~1.5 h GPU |

### 5.4 X4: improvement methods

| # | Method | Specification | Compared with |
|---|---|---|---|
| a | Early fusion | Concatenate per-encoder z-scored means (UNI + CONCH + TITAN; then + X1 winner), LR with C tuned on inner folds | Best single encoder |
| a | Late fusion / ensemble | Average of logits from the single-encoder LR and ABMIL models | Same |
| b | Multi-task ABMIL | Shared `fc` + attention. Heads: ER, PR, HER2_STRICT (BCE), PAM50 (5-class CE), ESR1/PGR/ERBB2 (MSE on z-scores). Missing targets masked; equal loss weights (sensitivity: uncertainty weighting, Kendall et al. 2018) | Single-task ABMIL per task |
| c | Site-invariant features | (i) Per-site feature standardisation: each hospital's slides are z-scored with that hospital's own mean/SD, using unlabeled features only, which a hospital can do for its own slides. (ii) Domain-adversarial ABMIL: gradient-reversal site classifier on training sites, λ ramped 0 → 0.1 | Same head without it; outcomes are AUC **and** `SITE10` accuracy of the adjusted features |
| d | Noisy-label training | (i) Train only on IHC = RNA patients (consensus). (ii) Soft labels = latent-class posterior P(ER+) from experiment C, **fitted on the training fold only**. (iii) Generalised cross-entropy (q = 0.7). (iv) Label smoothing 0.1. Scored against ER_IHC **and** ESR1_RNApos | Plain BCE on IHC |
| e | Multi-magnification | Concatenate UNI 10x and 20x means (TRACE UNI 10x, 16 GB, streamed) | UNI 20x |
| f | Calibration | Temperature scaling fitted on inner folds; ECE and Brier before / after; reliability diagrams | Uncalibrated |
| g | *LoRA on UNI (stretch)* | Rank 8 on qkv of the last 8 blocks; 512 patches per bag per step; gradient checkpointing | Frozen UNI ABMIL. Only if X1–X3 finish with ≥ 3 days left |

Scripts: `scripts/bench_fusion.py` (a, e), `scripts/bench_mil.py --multitask | --dann | --loss gce|smooth|soft`
(b, c, d) and `scripts/calibrate.py` (f).

### 5.5 X5: interpretability and representation analysis

| # | Analysis | Specification | Output |
|---|---|---|---|
| a | t-SNE / UMAP | Patient means per encoder, PCA → 50 → t-SNE (perplexity 30, 3 seeds) and UMAP (15 neighbours, min_dist 0.1). Colour by ER (IHC, RNA), PAM50, hospital; the 54 discordant patients as crosses | `figures/x5a_<enc>_{er,pam50,site}.png` |
| a | Quantities behind the maps | In the **full** feature space (not 2-D): 10-NN label purity for ER_IHC, ESR1_RNApos and PAM50; 10-NN hospital purity; silhouette by PAM50; for each discordant patient, the fraction of its 10 neighbours on the RNA side (vs the IHC side), tested against 0.5 | `results/bench/x5/embedding_metrics.csv` |
| b | Exact linear contribution map | For the scaled linear probe, logit = b + (1/N) Σᵢ cᵢ with cᵢ = Σⱼ (wⱼ/σⱼ)(xᵢⱼ − μⱼ). Each patch's cᵢ is drawn at its coordinates (red = towards ER+, blue = ER−). Unit test: the mean of the cᵢ plus b equals the model logit (to 1e-4) | Heatmap PNG per slide per encoder |
| b | Attention / instance maps | From `aux` of each X2 method (ABMIL, CLAM, DSMIL, ACMIL) | PNG per slide |
| b | Concept maps | Per-patch CONCH concept shares (from `conch_zeroshot.py`) for tubules, pleomorphism, necrosis, apocrine | PNG per slide |
| b | PCA-RGB map | PCA fitted on a 200k-patch cohort sample per encoder; the first 3 components → RGB at each patch | PNG per slide, UNI vs ResNet-50 vs CONCH |
| b | Map quantities | (i) Spearman correlation between the attention and contribution maps. (ii) Contribution share from patches in the top CONCH "tumour" concepts vs stroma / adipose. (iii) Deletion curve: remove top-k contributing patches, re-score, AUC drop vs random removal | `results/bench/x5/map_metrics.csv` |
| c | CKA | Linear CKA between encoders on patient means (same patients), bootstrap CI | Heatmap + CSV |
| d | Retrieval explanations | For a query patient: top-5 most similar *training-fold* patients (cosine on means) with their labels; top-8 most ER+ / ER− contributing patches across the training fold (gallery) | Used in the app |
| e | Concept probing | Ridge from each encoder's means to each of the 14 CONCH concept shares; out-of-fold R² per concept × encoder | Heatmap + CSV |
| f | *Sparse autoencoder (stretch)* | On 2M UNI patches: 8,192 latents, L1 1e-3, 5 epochs; latents vs ER contribution; top-activating patches per latent | Exploratory figure |

- **Slides for the maps:** the 10 patients in `results/attention_overlay/`, which include ER-low Basal and HER2-enriched
  exceptions, plus 4 concordant references.
- **Thumbnails:** come from IDC, as `attention_overlay.py` does now.
- **Script:** `scripts/maps.py` (needs bags + coords, so it runs on Ada).

### 5.6 X6: reliability

| Item | Specification |
|---|---|
| Uncertainty scores | (i) Predictive entropy of the deep ensemble (seeds × the top 3 encoders). (ii) MC-dropout (20 passes) in ABMIL. (iii) Patch-resampling SD: 20 random half-bags per slide |
| Conformal prediction | Split conformal per outer fold. Calibration set = one inner fold (site-grouped) of the training patients. Score 1 − p(true class). α = 0.10 → prediction sets {+}, {−}, {+, −}. Mondrian (class-conditional) variant as sensitivity |
| Metrics | Empirical coverage per outer fold and per hospital; mean set size; singleton rate; accuracy on singletons; risk–coverage curve and AURC; selective AUC at 80% / 90% coverage |
| Checks against the finished results | Odds ratio (Fisher) that an abstained case is an IHC/RNA-discordant patient; that it is one of the 1e TS/BS flips (frozen features); coverage on frozen TS slides with DX-calibrated sets (links to G) |
| Script | `scripts/bench_reliability.py` (CPU on saved OOF + ensemble; MC-dropout and half-bags on Ada) |

---

## 6. Code organisation

| Path | Status | Content |
|---|---|---|
| `configs/encoders.json` | new | Encoder registry (§3.3) |
| `scripts/bench/` (`registry, tasks, cohort, folds, heads, metrics, stats, io`) | new | Shared harness (§4) |
| `scripts/mil_models.py` | new | ABMIL (copied from `train_mil.py`, which stays unchanged), CLAM-SB, TransMIL, DSMIL, ACMIL, multi-task and DANN wrappers |
| `scripts/bench_encoders.py`, `bench_mil.py`, `bench_regression.py`, `bench_fusion.py`, `calibrate.py`, `bench_reliability.py`, `maps.py`, `embeddings.py` (X5a, c, e) | new | One per experiment |
| `scripts/fetch_rna_panel.py`, `make_dx_manifest.py` | new | Inputs |
| `scripts/build_bags.py` | extend | `--keep-coords`, `std` output, `--encoder-dir` for TRACE folder names |
| `scripts/extract_frozen.py` | extend | `--encoder resnet50 / dinov2_natural`; `--manifest` already generic |
| `scripts/sanity_check.py` | extend | S16 inner folds; S17 feature rows ↔ index ↔ labels aligned per encoder; S18 OOF complete (every patient predicted exactly once per seed); S19 contribution-map identity |
| `jobs/trace_features.sbatch`, `extract_dx.sbatch`, `bench_mil.sbatch`, `bench_cpu.sbatch`, `maps.sbatch` | new | Ada jobs |
| `app/` | new | Application (§7) |
| `results/bench/x1 … x6/`, `results/bench/figures/` | new | Outputs (small CSVs and PNGs committed; no features) |

**Backward compatibility.** Finished scripts are not modified except for the opt-in flags listed above. The ABMIL copy
in `mil_models.py`, run through the harness, must reproduce the fixed-split UNI ABMIL ER
OOF predictions. Check: rerun one seed and match `mil_uni_v1_20x_fixedsplits_abmil_oof.csv.gz` within GPU
non-determinism (AUC ±0.003).

---

## 7. Application: *ReceptorCheck*

### 7.1 What it does

For one patient's H&E slide (FFPE or frozen) it returns six things:
1. **ER / PR / HER2 estimates:** a calibrated probability, the predicted *ESR1/PGR/ERBB2* expression with a percentile
   against the training cohort, and a conformal prediction set.
2. **A reliability status:** `CONFIDENT`, `UNCERTAIN → do IHC`, or `OUT-OF-DISTRIBUTION → do not use`.
3. **An IHC consistency check:** when the lab's IHC result is entered, `CONSISTENT` or `FLAG FOR RETEST`, with a
   discordance score and its rank in the cohort.
4. **Evidence:** a contribution heatmap, an attention map, a concept profile (14 bars vs the cohort mean, z-scores)
   and the 5 most similar reference patients with their labels.
5. **Scope statements:** "HER2-low (0 vs 1+/2+) not assessable from H&E"; "results within ±0.1 of the threshold are
   unreliable (section-to-section variability)".
6. **A research-use-only banner.**

### 7.2 Architecture

```
offline (Ada / laptop)                                             online (laptop)
slide ─► QC + tiling (TRIDENT) ─► encoder (X1 choice) ─► predictor (X2–X4 choice)
                                         │                    │
                                         ├─► OOD score         ├─► calibration + conformal (X6)
                                         ├─► maps (X5b)        ├─► IHC consistency (D)
                                         └─► retrieval (X5d)   └─► regression percentiles (X3)
                     ▼
     scripts/app_build_cases.py ─► app/cases/<patient>.json + app/cases/<patient>/*.png ─► app/receptorcheck.py (Streamlit)
```

- **The demo uses out-of-fold outputs only.** Each patient is predicted by a model that did not train on that
  patient's hospital, so the demo shows honest, held-out behaviour.
- **Live mode** (optional, after the demo works): `app/run_slide.sh <slide>` submits an Ada job that extracts features
  for a new slide and scores it with a model trained on all TCGA patients. The result is written as a new case bundle.

### 7.3 Modules (`app/receptorcheck_core/`)

| Module | Interface | Logic |
|---|---|---|
| `predict.py` | `predict(feats) -> {task: {p, p_cal, set, expr, expr_pct}}` | Chosen predictor + temperature + conformal quantile per task (from X6 artefacts) |
| `reliability.py` | `status(pred, ood) -> CONFIDENT / UNCERTAIN / OOD` | OOD: Mahalanobis distance of the patient mean (PCA-50 space of the training cohort) above the training 99th percentile. UNCERTAIN: conformal set = {+, −} or \|p_cal − 0.5\| < 0.1 |
| `consistency.py` | `check(pred, ihc) -> {state, score, rank}` | score = p_cal(opposite of IHC). `FLAG` if the conformal set is a singleton contradicting IHC **and** score ≥ τ. τ is chosen in experiment D on training folds for precision ≥ 0.5 vs RNA-defined discordance at the cohort level |
| `explain.py` | `assets(patient) -> paths` | X5b maps, concept z-scores, X5d neighbours |
| `bundle.py` | `write(patient)`, `read(patient)` | JSON schema §7.4 |

### 7.4 Case bundle schema (`app/cases/<patient>.json`)

```json
{"patient": "TCGA-XX-XXXX", "slide": "…-DX1", "slide_type": "DX", "hospital_seen_in_training": false,
 "model": {"encoder": "<X1 choice>", "head": "<X2/X4 choice>", "version": "<git sha>"},
 "ood": {"score": 2.1, "threshold": 3.4, "flag": false},
 "tasks": {"ER": {"p": 0.12, "p_cal": 0.18, "set": ["neg"], "status": "CONFIDENT",
                  "expr": {"gene": "ESR1", "pred_log2": 8.9, "pct_vs_ERpos": 4},
                  "ihc": "pos", "consistency": {"state": "FLAG", "score": 0.82, "rank_pct": 97}},
           "PR": {"…": "…"}, "HER2": {"…": "…", "her2_low": "not assessable"}},
 "concepts": {"tubules": -1.2, "pleomorphism": 2.0, "necrosis": 1.1, "…": 0},
 "neighbours": [{"patient": "…", "ER_IHC": 0, "ESR1_RNApos": 0, "PAM50": "Basal", "sim": 0.94}],
 "assets": {"thumb": "thumb.png", "contrib_ER": "contrib_er.png", "attention_ER": "att_er.png", "concept_maps": ["…"]},
 "truth_for_evaluation_only": {"ESR1_RNApos": 0, "PAM50": "Basal"}}
```

### 7.5 User interface (Streamlit, single page)

1. **Sidebar:**
   - patient selector, with filters: all / discordant / ER-low / frozen available;
   - IHC entry (pre-filled from TCGA, editable);
   - slide type (DX / TS).
2. **Row 1:** three cards (ER, PR, HER2) with p_cal, the prediction set, the status colour and the predicted expression
   with its percentile.
3. **Row 2:** the IHC consistency result, with its score and rank and a one-line reason (e.g. "high pleomorphism and
   necrosis; neighbours are Basal").
4. **Row 3:** tabs: H&E thumbnail · contribution map · attention map · concept maps · similar cases.
5. **Footer:** scope statements and the research-use banner. In "evaluation mode", the hidden molecular truth is shown
   for discussion.

Dependencies: `streamlit`, `pandas`, `numpy`, `pillow`, all pinned in `app/requirements.txt`. The app runs on the
laptop and reads only small files (bundles about 1–3 MB per patient with images; about 200 bundles for the demo set).

### 7.6 How the application is evaluated (retrospective, out-of-fold, site-grouped)

| Use case | Metric | Baseline |
|---|---|---|
| Triage | At abstention rates of 10 / 20 / 30%: accuracy, sensitivity and specificity on kept cases vs IHC and vs RNA; IHC tests "saved" = kept cases | No abstention |
| QA second reader | Precision@k and recall@k of `FLAG` for RNA-defined discordance (k = 1, 2, 5% of cohort); enrichment for ER-low | Random flagging (base rate 5.5%); flagging by IHC % stained alone |
| Frozen section | ER AUC and calibration on TS with the DX-trained model; conformal coverage on TS (G, X6) | DX performance |
| Robustness | Per-hospital AUC spread; OOD flag rate on held-out hospitals vs in-distribution | – |

**Acceptance for the demo:**
- every number in the UI traces to an OOF file;
- coverage within ±3% of the target on DX;
- `FLAG` precision above the base rate with a CI excluding it;
- the UI loads any bundle in < 2 s.

---

## 8. How the planned Week 4 work plugs in (unchanged content)

| Planned item | Role here |
|---|---|
| Extended statistics (TOST, Bayes, Wilson, McNemar, calibration) | `bench/stats.py` reuses `stats_extended.py`; calibration feeds X4f / §7 |
| **C** latent-class "truth" | Soft labels for X4d (fitted per training fold); an extra evaluation column for the leaderboard |
| **D** model as label auditor | Defines the consistency score and threshold τ in `consistency.py` |
| **G** DX ↔ frozen transfer | Frozen-slide mode of the app; conformal coverage on TS |

---

## 9. Compute and timeline

| Day | Work | Ada | Gate to pass |
|---|---|---|---|
| 1 | Harness (`scripts/bench/`), registry, inner folds + S16–S18; TRACE listing; `dx_manifest.csv`; `fetch_rna_panel.py`; ABMIL copy + reproduction check | CPU | S16–S18 pass; ABMIL reproduces (AUC ±0.003) |
| 1–3 | Start ResNet-50 DX extraction (3 array batches, background); stream UNI2-h, GigaPath, CTransPath, CHIEF, UNI 10x | 3 GPU + CPU jobs | Storage gate (§3.4); 5-slide equivalence check |
| 2–3 | X1 mean-pool + label efficiency on available encoders; X3 ridge; X5a, c, e | CPU | – |
| 3–5 | X2 methods on UNI (6 methods × ~3 h); X3 ABMIL regression; X4a, e, f | 3 GPU | – |
| 5–6 | X1 with ResNet-50 complete; X1 ABMIL on top 2; X4b, c, d (after C) | 3 GPU | – |
| 6–7 | X5b maps; X6 reliability | GPU + CPU | Contribution-map identity S19 |
| 7–9 | App: case bundles, Streamlit UI, application evaluation §7.6 | laptop | Demo acceptance §7.6 |
| 9–10 | Report tables and figures; plan / log / handoff updates | – | – |
| (stretch) | X4g LoRA; X5f SAE; DINOv2 natural | GPU | Only with ≥ 3 days left |

**Slurm budget:** at most 3 GPU jobs at once (the 40-CPU cap with `-c 6` → 12 CPUs each) and ≤ 20 submitted jobs
including array tasks. Peak memory for DX extraction is 17.5 GB (`--mem-per-cpu=3G` is fine).

---

## 10. Risks and mitigations

| Risk | Mitigation |
|---|---|
| TRACE folder names or dims differ from the registry | Day-1 listing; registry filled from the listing; S17 checks the dims |
| Home quota | Means-only for encoders not used in MIL; raw parquet only on `/scratch`; storage gate |
| UNI2-h / GigaPath weights gated, so no equivalence check | TRACE features used as provided (as for UNI v1); noted as a limitation |
| ResNet-50 extraction slower than estimated (IDC throughput) | Report ResNet-50 only on the full cohort, never a partial subset. If extraction would take > 4 days, extract ResNet-50 on a fixed random half of patients (all encoders then compared on that half) and say so |
| Too many comparisons | Secondary labelling, full tables, "≥ 2 encoders or both heads" rule (§4.5) |
| Leakage through inner tuning or calibration | All tuning and calibration inside training folds; latent-class soft labels per training fold; S16 |
| Interpretability over-read | Every map paired with a quantity (X5b map metrics, deletion curves) |
| App mistaken for a clinical tool | Research-only banner; OOF-only predictions; scope statements in the UI |

---

## 11. Open choices (defaults used unless changed)

| Choice | Default |
|---|---|
| UI framework | Streamlit (Gradio is equivalent) |
| DINOv2 natural-image encoder | Stretch only |
| Encoders with bags (MIL) | UNI, CONCH and the X1 winner |
| PAM50 gene panel regression | Secondary; *ESR1/PGR/ERBB2* primary |
| Conformal α | 0.10 |
| Demo patient set | All 54 ER-discordant + 10 ER-low + 150 random concordant patients |
