# Baseline results (preliminary)

**Setup.** TRACE TITAN slide embeddings (768-d) for 1,125 TCGA-BRCA DX slides → mean per patient →
standardised L2 logistic regression (C = 0.01, class-balanced). Out-of-fold predictions from 5-fold CV, either
random or grouped by tissue source site (40 sites). 95% CIs from 1,000–2,000 bootstrap resamples.
Laptop CPU, ~10 min. Reproduce: `scripts/baseline_titan.py`, `scripts/paired_ihc_vs_rna.py` (last run 12 Oct 2026).

These are **linear probes on slide-level embeddings**, not the final ABMIL models.

## B0–B4: data and sanity checks

| # | Check | Result | Verdict |
|---|---|---|---|
| B0 | Label audit (IHC + RNA + CNA + RPPA + PAM50) | 1,097 patients; see DATASET.md | pass |
| B1 | Label-only agreement | ESR1 mRNA → ER IHC 0.958; ERBB2 mRNA → strict HER2 0.873; ERBB2 mRNA, HER2 0 vs 1+ 0.652 | disagreement is real |
| B2 | Shuffled ER labels | 0.46 (site CV), 0.43 (random CV) | pass (below 0.5 = single-seed fold artefact) |
| B3 | Site-only model, no image | ER IHC 0.551 · ESR1 RNA 0.529 · HER2 strict 0.590 | pass, HER2 warning |
| B4 | Site identification from the embedding (10 largest sites, n = 772) | 81.6% vs 18.4% chance | strong site signal → site-grouped CV |

## Slide-level AUC, random vs site-grouped CV (`results/baseline_results.csv`)

| Label | n | pos | AUC random (95% CI) | AUC site-grouped (95% CI) |
|---|---|---|---|---|
| ER (IHC) | 1,004 | 777 | 0.893 (0.86–0.92) | 0.885 (0.86–0.91) |
| ER (ESR1 RNA) | 1,039 | 788 | 0.936 (0.92–0.95) | 0.921 (0.90–0.94) |
| PR (IHC) | 1,001 | 671 | 0.804 (0.77–0.83) | 0.792 (0.76–0.82) |
| PR (PGR RNA) | 1,039 | 587 | 0.797 (0.77–0.82) | 0.787 (0.76–0.82) |
| HER2 (strict clinical) | 590 | 139 | 0.733 (0.69–0.78) | 0.674 (0.62–0.73) |
| HER2 (ERBB2 amplified) | 1,027 | 119 | 0.757 (0.71–0.80) | 0.755 (0.71–0.80) |
| HER2 IHC 0 vs 1+ | 318 | 261 | 0.633 (0.55–0.71) | **0.433 (0.36–0.52)** |
| Shuffled ER (control) | 1,004 | 779 | 0.431 (0.39–0.47) | 0.460 (0.42–0.50) |

**HER2 0 vs 1+ collapses to chance under site-grouped CV: the apparent HER2-low signal is site leakage.**

## B5: preliminary 1a (same patients, same site-grouped folds, paired bootstrap)

| Label | n | AUC IHC-label model | AUC RNA/DNA-label model | Δ (95% CI) |
|---|---|---|---|---|
| ER | 989 | 0.882 | 0.930 | **+0.048 (+0.023, +0.074)** |
| PR | 986 | 0.789 | 0.799 | +0.009 (−0.022, +0.041) |
| HER2 (strict vs ERBB2 amp) | 568 | 0.689 | 0.736 | +0.047 (−0.015, +0.104) |

## B6: preliminary 1b (train on ER-concordant, score ER-discordant)

| Subset | n | Sides with RNA | Sides with IHC | Binomial p |
|---|---|---|---|---|
| All discordant | 54 | 39 | 15 | 0.0015 |
| ESR1 > 0.5 log2 from cut-off | 45 | 36 | 9 | 6.6e-5 |
| ESR1 > 1.0 log2 from cut-off | 44 | 35 | 9 | 1.1e-4 |

## Caveats (B0–B6)

- One encoder (TITAN), one RNA cut-off scheme (GMM), one CV seed.
- 57 ER-discordant patients exist; 54 have TITAN embeddings.

---

# Extended TITAN baselines (E1–E5): 5 seeds, site-grouped CV

`python scripts/extended_titan.py` (12 Oct 2026) → `results/extended_titan.csv`. RPPA labels = 2-component GMM on
RPPA protein (ER-α, HER2). Seed SDs are across 5 different fold assignments; CIs are paired bootstrap on seed 0.

## E1: per-label AUC (mean ± SD over 5 seeds)

| Label | n | AUC |
|---|---|---|
| ER (IHC) | 1,004 | 0.882 ± 0.004 |
| ER (ESR1 RNA) | 1,039 | **0.929 ± 0.006** |
| ER (RPPA protein) | 852 | 0.891 ± 0.004 |
| PR (IHC) | 1,001 | 0.804 ± 0.007 |
| PR (PGR RNA) | 1,039 | 0.792 ± 0.004 |
| HER2 (strict clinical) | 590 | 0.690 ± 0.015 |
| HER2 (ERBB2 amplified) | 1,027 | 0.749 ± 0.005 |
| HER2 (RPPA protein) | 852 | 0.719 ± 0.013 |
| HER2 IHC 0 vs 1+ | 318 | 0.447 ± 0.030 |
| HER2 IHC 0 vs low (1+ or 2+/FISH−) | 454 | 0.418 ± 0.024 |
| Shuffled ER (control) | 1,004 | 0.462 ± 0.011 |

## E2 / E4: 1a paired comparisons (same patients, same folds)

| Comparison | n | Ref AUC (IHC) | Δ mean ± SD (5 seeds) | Seed-0 95% CI |
|---|---|---|---|---|
| ER: RNA − IHC (patients with all 3 ER labels) | 801 | 0.869 | **+0.055 ± 0.004** | +0.034, +0.088 |
| ER: RPPA protein − IHC | 801 | 0.869 | +0.023 ± 0.006 | −0.005, +0.061 |
| ER: RNA prevalence-matched cut-off − IHC | 989 | 0.881 | **+0.051 ± 0.005** | +0.035, +0.085 |
| ER: RNA confident-only (GMM posterior > 0.9 / < 0.1) − IHC | 945 | 0.892 | **+0.054 ± 0.003** | +0.035, +0.085 |
| PR: RNA − IHC | 986 | 0.800 | −0.007 ± 0.006 | −0.044, +0.020 |
| HER2: ERBB2 amplification − strict clinical | 478 | 0.678 | +0.041 ± 0.009 | −0.024, +0.108 |
| HER2: RPPA protein − strict clinical | 478 | 0.678 | +0.018 ± 0.011 | −0.056, +0.087 |

## E3 / E4: 1b which side on discordant cases

| Pair | Discordant n | Sides with molecular | Sides with IHC | Binomial p |
|---|---|---|---|---|
| ER IHC vs ESR1 RNA (GMM) | 54 | **39** | 15 | 0.0015 |
| ER IHC vs ESR1 RNA (prevalence-matched) | 61 | **42** | 19 | 0.0044 |
| ER IHC vs ESR1 RNA (confident-only) | 45 | **34** | 11 | 0.0008 |
| ER IHC vs RPPA protein | 63 | 33 | 30 | 0.80 |
| PR IHC vs PGR RNA | 169 | 69 | **100** | 0.021 |
| HER2 strict vs ERBB2 amplification | 60 | 35 | 25 | 0.25 |

## E5: 1d HER2-low (same patients)

| Task | n (pos) | H&E (TITAN) AUC | ERBB2 mRNA alone AUC |
|---|---|---|---|
| HER2 0 vs 1+ | 313 (257) | 0.464 ± 0.038 | 0.649 |
| HER2 0 vs low | 447 (391) | 0.398 ± 0.046 | 0.691 |

## What the extended baselines say (preliminary, linear probe on TITAN)

- **ER: the slide tracks mRNA, robustly.** The RNA-label advantage (+0.05 AUC) and the RNA side in 1b survive all
  three cut-off schemes. RPPA protein sits in between and does not separate from IHC.
- **PR goes the other way:** no 1a difference, and on discordant cases the model sides with **IHC** (100 of 169,
  p = 0.02). 138 of the 172 PR-discordant cases are IHC+/RNA−, so the PGR cut-off needs the same sensitivity check.
- **HER2:** amplification looks more predictable than the clinical label, but CIs cross zero (n = 478).
- **HER2-low:** H&E is at or below chance under site-grouped CV; ERBB2 mRNA alone reaches 0.65–0.69.
  Below-chance values (0.40) suggest site-linked label prevalence that flips between folds.

---

# Mean-pool baselines on UNI v1 and CONCH v1 patch features (3 seeds, site-grouped CV)

Features: TRACE patch features for 1,055 patients (DX slides), patient feature = mean of all patches.
UNI: `results/mil_uni_v1_20x_summary_meanonly.csv` (Ada job 17536); CONCH: `results/mil_conch_v1_20x_summary_meanonly.csv`.

| | TITAN (5 seeds) | UNI v1 mean | CONCH v1 mean |
|---|---|---|---|
| 1a ER: IHC → RNA label (n = 801) | 0.869 → 0.924 | 0.873 → 0.908 | 0.843 → 0.866 |
| 1a PR: IHC → RNA label, GMM cut-off (n = 986) | 0.800 → 0.793 | 0.795 → 0.733 | 0.778 → 0.748 |
| 1a HER2: strict → amplification (n = 478) | 0.678 → 0.719 | 0.627 → 0.679 | 0.623 → 0.679 |
| 1b ER discordant (n = 54), RNA : IHC | 39 : 15 (p = 0.0015) | 34 : 20 (p = 0.08) | 35 : 19 (p = 0.04) |
| 1b PR discordant (n = 169), RNA : IHC | 69 : 100 (p = 0.02) | 55 : 114 (p < 0.001) | 61 : 108 (p < 0.001) |
| 1b HER2 discordant (n = 60), amp : clinical | 35 : 25 (p = 0.25) | 44 : 16 (p < 0.001) | 39 : 21 (p = 0.03) |
| Shuffled ER control | 0.46 | 0.51 | 0.46 |
| HER2 0 vs low (ERBB2 mRNA alone 0.69) | 0.40 | 0.50 | 0.43 |

# Robustness across encoders (`scripts/robustness_mean.py` → `results/robustness_mean.csv`)

## R1: the PR "sides with IHC" result is a cut-off artefact

PGR labels: GMM cut-off (614 positive) vs **prevalence-matched** cut-off (log2 ≥ 6.95, 724 positive, matching
PR IHC prevalence) vs confident-only (GMM posterior > 0.9 / < 0.1, n = 903).

| Encoder | PGR label | 1a Δ (RNA − IHC) | 1b n | RNA : IHC | p |
|---|---|---|---|---|---|
| TITAN | GMM | −0.004 | 169 | 69 : 100 | 0.021 |
| TITAN | prevalence-matched | +0.039 | 133 | 69 : 64 | 0.73 |
| TITAN | confident-only | +0.009 | 105 | 50 : 55 | 0.70 |
| UNI mean | GMM | −0.047 | 169 | 55 : 114 | < 0.001 |
| UNI mean | prevalence-matched | +0.003 | 133 | 61 : 72 | 0.39 |
| UNI mean | confident-only | −0.028 | 105 | 44 : 61 | 0.12 |
| CONCH mean | GMM | −0.030 | 169 | 61 : 108 | < 0.001 |
| CONCH mean | prevalence-matched | +0.010 | 133 | 65 : 68 | 0.86 |
| CONCH mean | confident-only | −0.000 | 105 | 43 : 62 | 0.08 |

**Conclusion:** with a threshold that matches IHC prevalence, the PR preference disappears in all three encoders.
The GMM threshold for PGR is too high (it calls many low-PGR IHC+ tumours RNA-negative). ER, by contrast, is
robust to all three cut-off schemes (E4 above). ABMIL check with the prevalence-matched PGR label: Ada jobs 17660/17661.

## R2: which ER-discordant cases drive "sides with RNA"

| Encoder | IHC+/RNA− (n = 34) | IHC−/RNA+ (n = 20) | ER-low < 10% (n = 11) | 10–39% (n = 9) | ≥ 40% (n = 4) | % unknown (n = 30) |
|---|---|---|---|---|---|---|
| TITAN | 24 : 10 (p = 0.02) | 15 : 5 (p = 0.04) | 9 : 2 | 7 : 2 | 3 : 1 | 20 : 10 |
| UNI mean | 21 : 13 | 13 : 7 | 9 : 2 | 6 : 3 | 1 : 3 | 18 : 12 |
| CONCH mean | 22 : 12 | 13 : 7 | 9 : 2 | 7 : 2 | 1 : 3 | 18 : 12 |

(RNA side : IHC side.) **ER-low tumours side with RNA in 9 of 11 cases in every encoder**, consistent with current
guidelines that treat ER-low as biologically ER-negative. Both discordance directions lean to RNA; subgroups are small.
