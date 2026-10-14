# Baseline report: Option 1 "Which truth does the slide see?"

*Status as of 14 Oct 2026 (Section 5d added: DeLong tests, combined 1b evidence, discordant-case description). All baseline runs are complete and have passed the sanity checks. Every number below can be regenerated
with the scripts listed in Section 7.*

## 1. Summary

**For ER, every model we tried predicts the RNA-defined label better than the clinical IHC label, and on discordant
patients it sides with RNA.** For HER2 the models track DNA amplification. For PR there is no reliable preference:
an apparent IHC preference disappears once the RNA cut-off is chosen fairly. HER2 0 vs 1+ is not visible in H&E at all.

| Question | Result | Robust? |
|---|---|---|
| 1a ER: is the RNA label more predictable than IHC? | **Yes, +0.03 to +0.06 AUC** in all 5 model variants | Paired 95% CI excludes 0 for all 4 patch-feature models and for TITAN; survives 3 RNA cut-off schemes |
| 1b ER: which side on IHC ≠ RNA patients (n = 54)? | **RNA, 34–39 of 54** in all 5 variants | Significant (p < 0.05) in 3 of 5; **patient-level consensus of all 5: 34 : 20, p = 0.076 (suggestive)**; 66% of calls side with RNA (CI 54–76%) |
| 1a ER under DeLong | RNA − IHC CI excludes 0 in all 5 variants (p ≤ 0.038) | DeLong and bootstrap CIs agree (5d) |
| Who are the ER-discordant patients? | Call follows the intrinsic subtype: 15/15 Basal IHC+/RNA− called ER−, 9/9 LumA IHC−/RNA+ called ER+; HER2-enriched IHC+/RNA− mostly called ER+ (8/11) | Descriptive (n = 54) |
| 1a/1b ER vs RPPA protein | RPPA sits between IHC and RNA; no side preference (30–35 of 63) | Consistent |
| 1a/1b PR | GMM cut-off: RNA *harder*, model sides with IHC. **Prevalence-matched cut-off: no difference, no side** | The GMM result is a cut-off artefact (all 5 variants) |
| 1a HER2: amplification vs clinical label | Amplification more predictable in all variants (+0.02 to +0.09) | CI excludes 0 for CONCH (mean, ABMIL); not for UNI |
| 1b HER2: side on discordant (n = 60) | **Amplification, 34–47 of 60** in all 5 variants | Significant in 3 of 5; consensus 42 : 18, p = 0.003 |
| 1d HER2 0 vs 1+ / 0 vs low | **Chance (0.40–0.55)**; ERBB2 mRNA alone 0.65 / 0.69 | All variants |
| Negative control (shuffled ER) | 0.46–0.51 | All variants |

## 2. Setup

| Item | Value |
|---|---|
| Cohort | TCGA-BRCA, FFPE diagnostic (DX) slides, **1,055 patients** with features |
| Encoders | **TITAN** (slide embedding, 768-d), **UNI v1** (patch, 1,024-d), **CONCH v1** (patch, 512-d), all frozen; TRACE precomputed features |
| Heads | **Logistic regression** on the TITAN embedding or on the mean of all patches (C = 0.01, class-balanced); **ABMIL** (gated attention, 256-d hidden, dropout 0.25, Adam 2e-4, 15 epochs, 2,048 random patches per step, up to 4,096 patches per patient at test, fp16) |
| Cross-validation | 5-fold, **grouped by tissue-source site** (40 hospitals), stratified; out-of-fold predictions |
| Repeats | TITAN 5 seeds; UNI and CONCH 3 seeds (mean ± SD across seeds) |
| CIs | Paired bootstrap (2,000 resamples) of AUC(molecular) − AUC(IHC) on seed 0 |
| 1b | Model trained on all concordant patients, scored on discordant ones; ABMIL probabilities averaged over 3 seeds; two-sided binomial test vs 50:50 |
| Labels | IHC from GDC clinical biotab; RNA (ESR1, PGR, ERBB2) from cBioPortal, positive = high component of a 2-component GMM; HER2 strict clinical = IHC 3+ or FISH+ vs IHC 0/1+ or 2+ FISH−; amplification = GISTIC 2; RPPA = GMM on protein level |
| Compute | Laptop (TITAN, robustness, sklearn 1.5.2) and Ada RTX 2080 Ti (UNI/CONCH mean + ABMIL, sklearn 1.9.1); ABMIL 4.2 h (UNI) / 4.0 h (CONCH) per full run |

## 3. Sanity checks

`scripts/sanity_check.py` → `results/sanity_check.txt`. Each check ran for both UNI and CONCH, mean-pool and ABMIL.

| # | Check | UNI | CONCH |
|---|---|---|---|
| S1 | Each out-of-fold block has exactly one prediction per patient; no NaN; probabilities in [0, 1] (48 blocks) | PASS | PASS |
| S2 | AUC recomputed from stored predictions equals the summary (max diff 1e-16) | PASS | PASS |
| S3 | 1a is truly paired: all labels in a group are predicted for the same patients, per seed | PASS | PASS |
| S4 | Labels stored with predictions equal the label table (0 mismatches) | PASS | PASS |
| S5 | Shuffled-label control at chance: mean 0.506 / ABMIL 0.487 (UNI), 0.490 / 0.463 (CONCH) | PASS | PASS |
| S6 | Mean-pool results identical between two runs | PASS (identical) | **Explained difference** (see below) |
| S7 | 1b side counts recomputed from stored probabilities equal the summary | PASS | PASS |
| S8 | Cohort sizes as expected: 1a n = 801 (ER), 986 (PR), 478 (HER2); 1b n = 54 / 63 / 169 / 60; PR prevalence-matched 1b n = 133 | PASS | PASS |
| — | All 4 Ada jobs COMPLETED, exit 0, no errors in logs (17536, 17537, 17660, 17661) | PASS | PASS |

**S6, CONCH.** The earlier CONCH mean-pool run was done on the laptop (scikit-learn 1.5.2), the new one on Ada
(scikit-learn 1.9.1). The 1b counts, which do not use cross-validation, are **identical** (35, 32, 61, 39), so features
and labels match. The cross-validated AUCs differ (by up to 0.08 on the small HER2-low task, ≤ 0.03 elsewhere) because
`StratifiedGroupKFold` shuffles differently in the two versions: the same seed gives different folds.

**Consequences:**
- Comparisons *within* a run are valid, because all labels share the same folds.
- Absolute AUCs on small tasks are fold-sensitive.
- **Action:** save fixed fold assignments to `data/splits/` and read them in every script, and pin scikit-learn.
- The CONCH mean-pool numbers on slide 4 came from the laptop run. The Ada numbers below differ slightly, but all
  conclusions are unchanged.

## 4. Results

### 4.1 Experiment 1a: which label is more predictable from H&E? (AUC, mean ± SD over seeds)

**ER (same 801 patients)**

| Model | IHC label | RNA label | Δ RNA − IHC | Seed-0 95% CI | RPPA label | Δ RPPA − IHC |
|---|---|---|---|---|---|---|
| TITAN linear (5 seeds) | 0.869 | **0.924** | **+0.055 ± 0.004** | +0.034, +0.088 | 0.892 | +0.023 ± 0.006 |
| UNI mean-pool | 0.873 ± 0.001 | **0.908** ± 0.003 | **+0.035 ± 0.002** | +0.011, +0.066 | 0.914 | +0.041 ± 0.002 |
| UNI ABMIL | 0.864 ± 0.006 | **0.908** ± 0.004 | **+0.044 ± 0.007** | +0.021, +0.079 | 0.896 | +0.032 ± 0.005 |
| CONCH mean-pool | 0.844 ± 0.001 | **0.876** ± 0.003 | **+0.031 ± 0.004** | +0.005, +0.061 | 0.876 | +0.032 ± 0.005 |
| CONCH ABMIL | 0.847 ± 0.006 | **0.896** ± 0.006 | **+0.049 ± 0.009** | +0.011, +0.071 | 0.869 | +0.022 ± 0.005 |

ER cut-off robustness (TITAN, 5 seeds): Δ RNA − IHC = +0.055 (GMM), +0.051 (prevalence-matched), +0.054 (confident-only).

**PR (same 986 patients)**

| Model | IHC label | RNA label, GMM cut-off | Δ (GMM) | RNA label, prevalence-matched | Δ (prevalence-matched) |
|---|---|---|---|---|---|
| TITAN linear | 0.800 | 0.793 | −0.007 ± 0.006 | 0.836 | +0.039 ± 0.006 |
| UNI mean-pool | 0.795 | 0.733 | −0.062 ± 0.001 | 0.789 | +0.003 ± 0.006 |
| UNI ABMIL | 0.798 | 0.740 | −0.058 ± 0.012 | 0.792 | **−0.006 ± 0.005** (CI −0.033, +0.036) |
| CONCH mean-pool | 0.776 | 0.750 | −0.026 ± 0.000 | 0.788 | +0.010 ± 0.005 |
| CONCH ABMIL | 0.784 | 0.755 | −0.029 ± 0.030 | 0.797 | **+0.011 ± 0.004** (CI −0.019, +0.048) |

The GMM threshold for PGR (log2 8.58, 614 positive) is much stricter than IHC (699 positive). With a threshold that
matches IHC prevalence (log2 6.95, 724 positive) the gap disappears. *(Prevalence-matched values for TITAN, UNI mean and
CONCH mean come from `robustness_mean.csv` (laptop run, 3 seeds); their IHC AUCs in that run were 0.797, 0.786 and 0.778,
so each Δ is computed against its own run's IHC AUC.)*

**HER2 (same 478 patients)**

| Model | Strict clinical | Amplification | Δ amp − clinical | Seed-0 95% CI | RPPA | Δ RPPA − clinical |
|---|---|---|---|---|---|---|
| TITAN linear | 0.678 | 0.719 | +0.041 ± 0.009 | −0.024, +0.108 | 0.696 | +0.018 ± 0.011 |
| UNI mean-pool | 0.627 | 0.679 | +0.052 ± 0.006 | −0.016, +0.114 | 0.665 | +0.038 ± 0.012 |
| UNI ABMIL | 0.671 | 0.692 | +0.021 ± 0.016 | −0.064, +0.076 | 0.648 | −0.023 ± 0.002 |
| CONCH mean-pool | 0.616 | 0.701 | **+0.085 ± 0.013** | +0.038, +0.157 | 0.667 | +0.051 ± 0.006 |
| CONCH ABMIL | 0.674 | 0.726 | **+0.053 ± 0.030** | +0.008, +0.143 | 0.730 | +0.056 ± 0.028 |

### 4.2 Experiment 1b: which side on discordant patients? (molecular : IHC)

| Discordant pair | n | TITAN | UNI mean | UNI ABMIL | CONCH mean | CONCH ABMIL |
|---|---|---|---|---|---|---|
| ER: IHC vs RNA | 54 | **39 : 15** (p = 0.0015) | 34 : 20 (0.076) | 34 : 20 (0.076) | **35 : 19** (0.040) | **36 : 18** (0.020) |
| ER: IHC vs RPPA | 63 | 33 : 30 (0.80) | 35 : 28 (0.45) | 34 : 29 (0.62) | 32 : 31 (1.00) | 30 : 33 (0.80) |
| PR: IHC vs RNA, GMM | 169 | 69 : 100 (0.021) | 55 : 114 (<0.001) | 51 : 118 (<0.001) | 61 : 108 (<0.001) | 70 : 99 (0.031) |
| PR: IHC vs RNA, prevalence-matched | 133 | 69 : 64 (0.73) | 61 : 72 (0.39) | 61 : 72 (0.39) | 65 : 68 (0.86) | 63 : 70 (0.60) |
| HER2: clinical vs amplification | 60 | 35 : 25 (0.25) | **44 : 16** (<0.001) | **47 : 13** (<0.001) | **39 : 21** (0.027) | 34 : 26 (0.37) |

ER robustness (TITAN): confident-only RNA 34 : 11 (p = 0.0008); prevalence-matched 42 : 19 (p = 0.004).

**Which ER cases drive it (mean-pool encoders):**

| Subgroup | TITAN | UNI mean | CONCH mean |
|---|---|---|---|
| IHC+/RNA− (n = 34) | 24 : 10 | 21 : 13 | 22 : 12 |
| IHC−/RNA+ (n = 20) | 15 : 5 | 13 : 7 | 13 : 7 |
| ER-low, < 10% cells stained (n = 11) | 9 : 2 | 9 : 2 | 9 : 2 |

### 4.3 Experiment 1d: HER2-low (same patients; H&E vs ERBB2 mRNA alone)

| Task | n (pos) | TITAN | UNI mean | UNI ABMIL | CONCH mean | CONCH ABMIL | ERBB2 mRNA |
|---|---|---|---|---|---|---|---|
| HER2 0 vs 1+ | 313 (257) | 0.464 | 0.466 | 0.503 | 0.506 | 0.465 | **0.649** |
| HER2 0 vs low | 447 (391) | 0.398 | 0.500 | 0.549 | 0.516 | 0.520 | **0.691** |

With random (not site-grouped) folds, TITAN scored 0.633 on 0 vs 1+: that apparent signal was hospital leakage.

### 4.4 Controls and preliminary checks

| Check | Result |
|---|---|
| Shuffled ER labels (AUC) | TITAN 0.462 · UNI mean 0.506 · UNI ABMIL 0.487 · CONCH mean 0.490 · CONCH ABMIL 0.463 |
| Label-only agreement, no images (AUC) | ESR1 → ER IHC 0.958 · PGR → PR IHC 0.921 · ERBB2 → strict HER2 0.873 |
| Hospital ID alone predicts label (AUC) | ER IHC 0.551 · ER RNA 0.529 · HER2 0.590 |
| Slide embedding identifies hospital (10 largest sites) | 81.6% vs 18.4% chance → site-grouped CV mandatory |

## 5. Interpretation

1. **ER: H&E morphology tracks ER gene expression more than the clinical stain.** The effect is consistent across three
   encoders, two classifier types and three RNA cut-offs, and it holds as a paired comparison on identical patients and
   folds. It is not an encoder or model artefact.
2. **The ER effect is not "protein amount".** RPPA protein sits between IHC and RNA and shows no side preference. The slide
   follows the transcriptional state of the tumour.
3. **ER-low is the likely driver.** Weakly stained IHC-positive tumours are called negative by the slide model (9 of 11),
   consistent with guidelines that treat ER-low as biologically ER-negative. Small n: a lead for the final report.
4. **PR: no preference.** The earlier "sides with IHC" was produced by a too-strict GMM threshold, and it vanishes in all
   5 variants with a prevalence-matched threshold. This is a methodological lesson worth reporting.
5. **HER2: the slide follows DNA amplification over the clinical call**, consistently in direction. Significance depends on
   the encoder. The clinical label set is small (n = 478, 60 discordant).
6. **HER2-low is invisible in H&E** under proper evaluation; random cross-validation would have suggested otherwise.
7. **ABMIL ≈ mean-pool.** Attention pooling does not change the conclusions, which shows the label effect is not tied to the classifier.

## 5b. Fixed-split rerun (13 Oct 2026): S6 resolved

All linear baselines were rerun on Ada in one environment (scikit-learn 1.9.1) using the new **fixed, library-independent,
site-grouped splits** `data/splits/folds.csv` (`scripts/make_splits.py`; every fold has ≥ 10 patients of each class for
every task, including HER2-low). Job 19737, 5 seeds. ABMIL results above used the earlier folds and were not rerun,
because the fixed-split linear results change no conclusion.

| Result | Before (old folds) | Fixed splits |
|---|---|---|
| TITAN 1a ER Δ RNA − IHC (n = 801) | +0.055 (CI +0.034, +0.088) | **+0.058** (CI +0.031, +0.084) |
| TITAN 1a ER Δ, prevalence-matched / confident-only cut-off | +0.051 / +0.054 | **+0.052 / +0.059** |
| TITAN paired ER (n = 989) | 0.882 → 0.930 | **0.884 → 0.937** (CI +0.031, +0.078) |
| UNI mean 1a ER Δ | +0.035 | **+0.037** |
| CONCH mean 1a ER Δ | +0.031 | **+0.030** |
| PR Δ, prevalence-matched (TITAN / UNI / CONCH mean) | +0.039 / +0.003 / +0.010 | **+0.033 / −0.008 / +0.010** |
| HER2 Δ amp − clinical (TITAN / UNI / CONCH mean) | +0.041 / +0.052 / +0.085 | **+0.044 / +0.071 / +0.073** |
| HER2 0 vs 1+ (TITAN / UNI / CONCH mean) | 0.46 / 0.47 / 0.51 | **0.45 / 0.50 / 0.51** |
| Shuffled control (TITAN / UNI / CONCH mean) | 0.46 / 0.51 / 0.49 | **0.46 / 0.51 / 0.50** |
| 1b counts (no CV involved) | | **identical** (ER 39/34/35 of 54; PR prevalence-matched 69/61/65 of 133; HER2 35/44/39 of 60) |

Files: `results/extended_titan.csv`, `results/robustness_mean.csv`, `results/baseline_results.csv`,
`results/paired_ihc_vs_rna.txt`, `results/mil_{uni,conch}_v1_20x_fixedsplits_summary.csv`.

## 5c. ABMIL attention maps for ER-discordant and ER-low patients (exploratory)

`scripts/attention_maps.py`, Ada job 19722 (1 h 42 min). ER ABMIL on UNI features was trained on ER-concordant patients
(3 seeds), holding out 10 concordant references (5 ER+, 5 ER−). For each of the 54 discordant patients and the 10 references,
the full set of patches (median ~13,000) was scored and the attention weights averaged over seeds.
Outputs: `results/attention/summary.csv`, one heatmap per patient (`results/attention/<patient>.png`) and
`results/attention/overview.png` (ER-low patients first). Colour = attention percentile within the patient.

| Group | n | Mean P(ER+) | Attention entropy (normalised) | Fraction of patches holding 50% of attention |
|---|---|---|---|---|
| Reference ER+ (concordant) | 5 | 0.83 | 0.82 | 6.9% |
| Reference ER− (concordant) | 5 | 0.22 | 0.85 | 7.4% |
| Discordant IHC+/RNA− | 34 | 0.40 | 0.81 | 7.1% |
| Discordant IHC−/RNA+ | 20 | 0.69 | 0.79 | 4.7% |

- **References:** 8 of 10 classified correctly (sanity check of the trained model).
- **Discordant:** sides with RNA in **35 of 54** (IHC+/RNA−: 21 of 34; IHC−/RNA+: 14 of 20), matching the 1b results.
- **ER-low (IHC+, < 10% stained):** sides with RNA in **8 of 10**. (The robustness table in 4.2 counts 11 because it also includes one IHC-negative
  patient whose record says < 10%; here ER-low requires IHC-positive.) The two
  exceptions (TCGA-3C-AALI, TCGA-AR-A0TX) are both PAM50 **HER2-enriched**; all six Basal-like ER-low tumours are called ER−.
- **P(ER+) tracks ESR1 expression** among discordant patients (r = 0.52): the model's output is graded with the tumour's
  ESR1 level, not just a binary flip.
- **Where attention goes:** visually it concentrates on dense, nodular central regions, with low attention on tissue edges
  and sparse areas. That is consistent with tumour regions, but it is **not yet verified against the H&E images** (the maps show
  patch positions only). Attention is moderately focused (about 5–7% of patches carry half of the attention) and similar across
  groups, so the model does not behave differently on discordant cases; it simply reads them as the molecular class.

## 5d. DeLong tests (1a), combined evidence across encoders (1b), who the discordant patients are (14 Oct 2026)

`scripts/titan_oof.py` + `scripts/stats_delong_1b.py`, Ada job 19989 (`jobs/stats_tier2.sbatch`); no new training except TITAN
out-of-fold predictions on the fixed splits. TITAN and mean-pool use the fixed splits; ABMIL predictions are from the earlier
folds (flagged "old"). Probabilities are averaged over seeds per patient (each seed's prediction is out-of-fold).

**DeLong.** The two AUCs in 1a use *different reference labels on the same patients*, so the covariance is computed from
per-patient DeLong influence terms; with identical labels this reproduces the classic DeLong test exactly (`--selftest`).
In every row below the DeLong 95% CI agrees with the paired bootstrap CI to within ±0.005.

| ΔAUC (molecular − IHC), seed-averaged | n | TITAN | UNI mean | CONCH mean | UNI ABMIL (old) | CONCH ABMIL (old) |
|---|---|---|---|---|---|---|
| **ER: RNA − IHC** | 801 | **+0.057** (+0.030, +0.084) p < 0.001 | **+0.036** (+0.007, +0.065) p = 0.014 | **+0.030** (+0.002, +0.058) p = 0.038 | **+0.034** (+0.007, +0.061) p = 0.014 | **+0.049** (+0.020, +0.077) p = 0.001 |
| ER: RPPA − IHC | 801 | +0.026, p = 0.12 | **+0.041**, p = 0.009 | **+0.033**, p = 0.023 | +0.023, p = 0.15 | +0.019, p = 0.29 |
| PR: RNA − IHC (GMM cut-off) | 986 | −0.006, p = 0.70 | **−0.048**, p = 0.007 | −0.027, p = 0.10 | **−0.064**, p < 0.001 | −0.032, p = 0.056 |
| PR: RNA − IHC (prevalence-matched) | 986 | **+0.033**, p = 0.032 | — | — | −0.006, p = 0.73 | +0.006, p = 0.72 |
| HER2: amplification − strict | 478 | +0.044, p = 0.17 | **+0.071**, p = 0.029 | **+0.074**, p = 0.015 | +0.024, p = 0.44 | +0.060, p = 0.067 |
| HER2: RPPA − strict | 478 | +0.010, p = 0.80 | +0.045, p = 0.26 | +0.044, p = 0.25 | −0.023, p = 0.53 | +0.064, p = 0.11 |

- **ER (primary claim): holds under DeLong in all 5 variants**, CIs excluding 0 (p ≤ 0.038; significant in 19 of 21 individual seeds).
- ER RPPA, PR (either cut-off) and HER2 are encoder-dependent: **not robust** under the decision rule (0.5). The single
  significant prevalence-matched PR result (TITAN, +0.033) is not replicated by UNI or CONCH ABMIL.

**Combined 1b.** The same discordant patients are scored by every variant, so the five counts are not independent and must
not be summed into one binomial test. Instead, one consensus call per patient (mean P over the 5 variants):

| Discordant pair | n | Per variant (TITAN / UNI m / CONCH m / UNI A / CONCH A) | Consensus molecular : IHC | Unanimous (mol / IHC) | Calls siding molecular (cluster-bootstrap CI) |
|---|---|---|---|---|---|
| **ER: IHC vs RNA** | 54 | 39 / 34 / 35 / 34 / 36 | **34 : 20** (p = 0.076) | 29 / 11 | **66%** (54–76%) |
| ER: IHC vs RPPA | 63 | 33 / 35 / 32 / 34 / 30 | 33 : 30 (p = 0.80) | 29 / 24 | 52% (41–64%) |
| PR: IHC vs RNA, GMM | 169 | 69 / 55 / 61 / 51 / 70 | 53 : 116 (p < 0.001) | 34 / 71 | 36% (30–42%) |
| PR: IHC vs RNA, prevalence-matched (3 variants) | 133 | 69 / – / – / 61 / 63 | 64 : 69 (p = 0.73) | 52 / 56 | 48% (41–56%) |
| **HER2: strict vs amplification** | 60 | 35 / 44 / 39 / 47 / 34 | **42 : 18** (p = 0.003) | 23 / 7 | **66%** (58–75%) |

- ER: every variant puts the majority on the RNA side, but pooled per patient the evidence is **suggestive, not significant**
  (34 : 20, p = 0.076). TITAN alone (39 : 15, p = 0.0015) is the most extreme variant; **report the consensus, not TITAN**.
  The fraction of calls siding with RNA (66%, CI excludes 50%) supports the same direction.
- HER2 combined: the model sides with amplification (42 : 18, p = 0.003), stronger than the 1a ΔAUC suggests.
- PR: no preference once the cut-off is prevalence-matched (confirms Section 4.2).

**Who the 54 ER-discordant patients are** (`scripts/describe_discordant.py` → `results/er_discordant_summary.txt`,
`results/er_discordant_cases.csv`; descriptive, n small):

| Direction | PAM50 | Consensus sides with RNA | sides with IHC |
|---|---|---|---|
| IHC+/RNA− (n = 34) | Basal-like | **15** | 0 |
| | HER2-enriched | 3 | **8** |
| | Normal-like / LumA / NA | 3 | 5 |
| IHC−/RNA+ (n = 20) | Luminal A | **9** | 0 |
| | LumB / Basal / HER2 / Normal / NA | 4 | 7 |

- The model's call follows the **intrinsic subtype**: every Basal-like IHC+/RNA− tumour is called ER− and every Luminal A
  IHC−/RNA+ tumour is called ER+. For reference, 72% of concordant ER− tumours are Basal-like and 68% of concordant ER+ are LumA.
- **HER2-enriched IHC+/RNA− tumours are the exception**: 8 of 11 are called ER+ (with IHC). The two ER-low exceptions in
  Section 5c are also HER2-enriched. A plausible reading is that HER2-enriched morphology is not "basal-looking" enough to
  be called ER−; to check with the heatmaps.
- Called-with-IHC discordant cases sit closer to the ESR1 cut-off (median ESR1 − cut-off −1.1 vs −2.6 log2 for IHC+/RNA−),
  i.e. the model disagrees with RNA mainly where RNA itself is borderline.
- ER-low (IHC+, < 10%): 10 patients (6 Basal, 3 HER2-enriched, 1 NA); consensus sides with RNA in 8.
- Histology and stage show no clear pattern (ductal 27 / lobular 4 / other 3 among IHC+/RNA−).

Interpretation for the report: on H&E, "ER status" is largely read through the tumour's intrinsic-subtype morphology
(basal vs luminal), which RNA reflects more faithfully than IHC in these discordant cases.

## 6. Caveats and next steps

| Caveat | Next step |
|---|---|
| ~~Fold assignment depends on the scikit-learn version (S6)~~ | **Done:** fixed splits in `data/splits/folds.csv`, scikit-learn pinned, linear baselines rerun (Section 5b) |
| Attention maps show patch positions, not tissue | Overlay on slide thumbnails (download a few DX slides) to confirm attention on tumour |
| Small discordant groups (ER 54, HER2 60, ER-low 11) | **Done (5d):** patient-level consensus across 5 variants; ER 34 : 20 (p = 0.076) is suggestive only, report it as such |
| HER2 significance encoder-dependent | **Done (5d):** DeLong per variant; HER2 1a not robust across encoders, but combined 1b 42 : 18 (p = 0.003) |
| CONCH mean-pool numbers on slide 4 are from the laptop run | Update slides with the Ada values (conclusions unchanged); change the PR bullet to "no reliable preference" |
| 1c and 1e not yet run | UNI access granted (Srivastava); manifest of 1,516 TS/BS slides built; extraction queued on Ada (Section 7) |
| Attention analysis is exploratory (n = 10 ER-low) | Report as a lead. **Done (5d):** discordant cases described by PAM50 / ER % / histology; check HER2-enriched exceptions in the heatmaps |

## 7. Files and how to reproduce

| What | File | Produced by |
|---|---|---|
| TITAN baselines (B0–B6) | `results/baseline_results.csv` | `scripts/baseline_titan.py`, `scripts/paired_ihc_vs_rna.py` |
| TITAN extended (E1–E5) | `results/extended_titan.csv` | `scripts/extended_titan.py` |
| Cut-off and ER-strata robustness | `results/robustness_mean.csv` | `scripts/robustness_mean.py` |
| UNI / CONCH mean + ABMIL | `results/mil_{uni,conch}_v1_20x_summary.csv`, `*_oof.csv.gz` | `jobs/train_mil.sbatch` → `scripts/train_mil.py` (Ada jobs 17536, 17537) |
| PR prevalence-matched ABMIL | `results/mil_{uni,conch}_v1_20x_prprev_summary.csv` | `train_mil.py --only pr --pr-cut prevmatch` (jobs 17660, 17661) |
| Sanity checks + paired CIs | `results/sanity_check.txt`, `results/paired_ci_patch.csv` | `scripts/sanity_check.py` |
| DeLong 1a, combined 1b | `results/stats_delong_1a.csv`, `results/stats_combined_1b.csv`, `results/stats_combined_1b_patients.csv` | `jobs/stats_tier2.sbatch` → `scripts/titan_oof.py`, `scripts/stats_delong_1b.py` (Ada job 19989) |
| TITAN out-of-fold predictions (fixed splits) | `results/titan_linear_{fixedsplits,prprev}_oof.csv.gz` | `scripts/titan_oof.py` |
| ER-discordant case description | `results/er_discordant_summary.txt`, `results/er_discordant_cases.csv` | `scripts/describe_discordant.py` |
| Frozen-slide manifest (1c, 1e) | `data/metadata/frozen_manifest.csv` | `scripts/make_frozen_manifest.py` |
| Job logs | `results/logs/17536.out` … `17661.out` | Ada |
