# Week 3: candidate foundation-model experiments that build on our findings

*Written 26 Oct 2026. Candidates only: pick which to run before starting (tick them in `docs/PROJECT_PLAN.md`, Tier 6b).*

**Starting point (what we know):** H&E foundation-model classifiers track the RNA label for ER more than the clinical IHC
label (P1), side with RNA on most discordant tumours (P2, frozen slides 39 : 13), and the answer does not change between
FFPE and frozen slides (P3). HER2 follows amplification. The model's call follows intrinsic subtype (Basal vs luminal),
with HER2-enriched tumours as the exception. Two frozen sections of one tumour disagree in 8-11% of calls, all near
the threshold. Embeddings carry hospital identity.

**Wider point these experiments would add:** in medical AI the "ground truth" is often one imperfect test among several.
Our data has three independent measurements of the same biology (IHC, RNA, protein) on the same patients, which is rare.
That lets us test ideas that most medical-AI papers cannot check, because they have only one reference standard.

**Budget:** about 1-2 weeks; Ada (3 × 11 GB GPUs at a time); UNI/CONCH/TITAN features and all out-of-fold predictions
already exist, so most of these are CPU work on saved outputs.

| # | Experiment | Question it answers | Uses | Effort | Novelty (literature search, 19 and 26 Oct) |
|---|---|---|---|---|---|
| A | **Label-free audit with zero-shot CONCH** (in plan) | If a model that never saw an IHC or RNA label also ranks tumours closer to RNA, then the RNA preference comes from morphology, not from the classifier fitting label noise | CONCH text + image towers on existing CONCH patch bags | 2 days, CPU (+15 min GPU check) | No zero-shot ER/HER2 audit against several reference standards found |
| B | **Concept bottleneck** (in plan, same job as A) | Which named morphology concepts (tubules, pleomorphism, mitoses, necrosis, TILs, apocrine, mucin…) drive the calls, and do they explain why HER2-enriched IHC+/RNA− tumours are called ER+? | CONCH concept scores → small interpretable model | with A | Concept bottlenecks exist; not used to explain a label-source preference |
| C | **Evaluation without a gold standard (latent class model)** | Treat IHC, RNA, protein and the H&E model as four imperfect "raters" of ER; estimate each one's sensitivity/specificity and each patient's probability of being truly ER+ (Dawid–Skene / Hui–Walter); re-score the models against that latent truth. Does the H&E model look better or worse than IHC once nobody is assumed perfect? | Saved predictions + labels | 1-2 days, CPU | Standard in diagnostic-test research; not found for pathology FMs |
| D | **The model as a label auditor, checked by an independent oracle** | Rank patients by how strongly the H&E model (and the 5-variant consensus) disagrees with IHC; how well does that find IHC≠RNA tumours (AUC, precision at top-k, ER-low enrichment)? Then drop or relabel the flagged cases and retrain: does IHC-label AUC rise? | Saved OOF predictions; mean-pool retraining | 1-2 days, CPU | Others re-stained flagged cases (Shamai 2024, ReceptorNet); RNA as an orthogonal oracle on a whole cohort not found |
| E | **Which label should you train on?** | Train on RNA labels, test against IHC (and vice versa), plus learning curves (50 → 800 training patients). If training on the cleaner label wins even when the hospital scores against IHC, that is a practical recommendation | Mean-pool features, fixed folds | 1 day, CPU | Closest: ESR1-expression training then clinical ER prediction (medRxiv 2024); paired same-fold comparison and label-efficiency curve not found |
| F | **Uncertainty that knows when the truth is ambiguous** | Uncertainty from patch resampling within a slide and from TS vs BS disagreement; risk–coverage curves scored against IHC and against RNA. Are the cases the model abstains on enriched for IHC/RNA discordance? | UNI bags (DX, frozen) on Ada | 2 days, CPU on Ada | Selective prediction is known; "abstention tracks reference disagreement" not found |
| G | **DX ↔ frozen transfer** | Train on FFPE, test on frozen (and back): how much does AUC drop, does calibration shift, and does the RNA preference survive the domain change? Relevant to intra-operative frozen-section AI | UNI means (DX + TS) on Ada | 0.5 day, CPU | Frozen ≈ FFPE accuracy shown (DeepHRD 2024); preference under transfer not found |
| H | **Does the reference standard change the leaderboard?** | Rank the 5 model variants under IHC, RNA and protein labels; Kendall τ and the bootstrap probability that the ranking flips. If it flips, benchmark results depend on which "truth" the benchmark used | Saved OOF predictions | 0.5 day, CPU | Not found |

**Recommendation.** A + B are already planned (running now). For new findings, add **C** (no-gold-standard
evaluation) and **D** (label auditor with an RNA oracle): both are cheap, use what we already have, and say something
general about evaluating medical AI when labels are imperfect. **H** and **E** are half-day add-ons with a clear
take-home message. Keep **F** and **G** for the time that is left after that. CPTAC stays on hold.

**Guardrails (same as the rest of the project):** fixed folds only; report n and a CI for every number; a claim needs
≥ 2 encoders and must survive the RNA cut-off schemes; subgroups with n < 20 are leads only; these are exploratory
(secondary) analyses and are labelled as such next to P1–P3.
