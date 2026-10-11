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

## Caveats

- One encoder (TITAN), one RNA cut-off scheme (GMM), one CV seed.
- 57 ER-discordant patients exist; 54 have TITAN embeddings.
