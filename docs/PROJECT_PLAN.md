# Project plan

## Research question

When IHC, RNA and PAM50 disagree, which one does an H&E foundation-model classifier follow, and does the answer
change between the FFPE diagnostic slide (DX) and the frozen slide (TS/BS) cut from the sequenced tissue?

Clinical relevance: before an H&E model is used as a second reader for receptor status (the Kamineni brief),
we must know whether it reproduces clinical IHC calls, molecular truth, or neither.

## Experiments (each standalone: its own model or test, none reuses another's output)

| ID | Design | Labels | Test | Patients | Possible findings (both informative) |
|---|---|---|---|---|---|
| **1a** | Same patients, same site-grouped folds; one model per label source | ER: IHC vs ESR1 RNA; PR: IHC vs PGR RNA; HER2: strict clinical vs ERBB2 amplification (+ RPPA protein) | Paired bootstrap / DeLong ΔAUC | ER 989 · PR 986 · HER2 568 | RNA > IHC: IHC is the noisier truth for morphology. IHC > RNA: morphology tracks protein |
| **1b** | Train only on IHC = RNA cases; score IHC ≠ RNA cases | ER, PR, HER2 | Binomial test on the side taken | Discordant: ER 57 · PR 172 · HER2 64 | Direct answer to "which truth" |
| **1c** | 2 × 2: {FFPE-DX, frozen-TS} × {IHC, RNA} | ER, PR, HER2 | Interaction (difference of ΔAUCs, bootstrap) | ~965 with DX + TS | Positive: tissue sampling location matters. Null: label mismatch is not a sampling effect |
| **1d** | HER2 IHC 0 vs 1+ (and 0 vs low = 1+ or 2+/FISH−) | HER2 IHC score | H&E model vs ERBB2-mRNA-only model, site-grouped only | 57 vs 261 (DX) | Is morphology or mRNA more separable on the hardest IHC call? |
| **1e** | Same trained model applied to paired TS and BS slides | ER, HER2 | ICC, Cohen's κ, % flipped calls | 374 pairs | First prediction-level within-patient reproducibility for breast receptors |

## Pipeline

```
Labels (no slides)  GDC clinical biotab (IHC) + cBioPortal (RNA, CNA, RPPA, PAM50)
                    -> scripts/build_labels.py -> data/metadata/labels_patient.csv

Features            DX slides:     TRACE precomputed UNI v1 (20x, 256 px) -> per-slide bags
                    Frozen TS/BS:  GDC download -> TRIDENT tissue seg + 256 px @ 20x -> UNI v1 (fp16, 2080 Ti)

Heads               bag -> ABMIL (gated attention) | mean-pool + logistic regression
                    patient-level, site-grouped 5-fold CV x 5 seeds

Statistics          paired bootstrap / DeLong, binomial, ICC / kappa, interaction via bootstrap

Robustness          CONCH v1 replication; RNA cut-off sensitivity (GMM / quartile / PAM50-matched);
                    ER-low vs ER-high strata; % tumour covariate for frozen slides
```

Details and step-by-step commands: [`PIPELINE_AND_MODELS.md`](PIPELINE_AND_MODELS.md).

## Timeline (proposed, adjust in comms/LOG.md)

| Week | Uniyal (DX track) | Srivastava (frozen track + stats) | Deliverable |
|---|---|---|---|
| 1 (to 12 Oct) | Dataset study, baselines on TITAN | Review handoff, request model access | Project update (done) |
| 2 (13–19 Oct) | TRACE UNI download + bag building; ABMIL code; 1a with ABMIL | GDC manifest of 1,517 tumour TS/BS slides; TRIDENT + UNI test on 1 slide, then batches | 1a (ABMIL), extraction running |
| 3 (20–26 Oct) | 1b for ER/PR/HER2; cut-off sensitivity; RPPA label | Finish frozen extraction; TRACE-equivalence check on ~50 DX slides; shared stats module | 1b, frozen features ready |
| 4 (27 Oct–2 Nov) | 1d; attention heatmaps on discordant cases | 1c and 1e | 1c, 1d, 1e first results |
| 5 (3–9 Nov) | CONCH replication | Repeated-CV + final statistics | Robustness tables |
| 6+ | Write-up, figures, final presentation | Write-up, figures, final presentation | Report + talk |

## Compute budget (Ada: 4× GTX 1080 Ti or 4× RTX 2080 Ti per node, 11 GB each, 1.8 TB scratch)

| Item | Estimate |
|---|---|
| TRACE UNI v1 20× download | 62 GB, 0 GPU-hours |
| Frozen TS/BS download | 430 GB in batches of ~300 slides (delete after encoding) |
| Frozen-slide extraction | ~6M patches; ~1 day wall time on one 4× 2080 Ti node |
| ABMIL grid (~10 labels × 2 slide types × 2 encoders × 5 seeds × 5 folds) | 40–80 GPU-hours |
| Peak scratch | ~550 GB |

## Risks

| Risk | Mitigation |
|---|---|
| RNA cut-offs arbitrary | Three cut-off schemes; 1b already holds after dropping borderline cases |
| Site leakage | Site-grouped CV everywhere; random CV only as a contrast |
| HER2 IHC score missing for 44% of patients; only 57 HER2-0 | Add "0 vs low"; frame 1d as exploratory |
| Frozen-slide artefacts confound 1c | % tumour covariate; extraction-equivalence check vs TRACE |
| Gated weights (UNI/CONCH need institutional e-mail) | Request today; TITAN/TRACE cover DX work meanwhile |
| Licences (CC-BY-NC-ND, TRACE unstated) | Academic use only; never commit features |

## Deferred

- **1f: LoRA fine-tuning of UNI** on IHC vs RNA labels (does adaptation change which truth is learned?).
  About 30–60 GPU-hours, uses frozen slides. Revisit after week 4.
