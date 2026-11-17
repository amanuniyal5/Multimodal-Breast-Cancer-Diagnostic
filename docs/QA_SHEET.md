# Q&A sheet for the presentation

Short, honest answers to the questions we expect, with the number to quote. Sources: `docs/FINAL_RESULTS.md` (final
tables, fixed splits; extensions in 1.1) and `docs/baseline_report.md` (details; extensions in 5k-5l).

## About the question and novelty

**Q: Isn't predicting ER/PR/HER2 from H&E already done?**
Yes, many times, always against one label (usually IHC). We do not claim a better predictor. Our contribution is the
evaluation: same patients, slides and folds against different reference standards, a discordant-case test, and FFPE vs
frozen slides from the same tumours. Nobody has asked which "truth" the model follows.

**Q: Why should a clinician care?**
Because IHC-trained H&E models are penalised when they agree with RNA, and IHC-only benchmarks understate how much ER
biology the slide carries. In low-expression tumours (ER-low, HER2-enriched) an H&E model should not be expected to
reproduce IHC.

**Q: Isn't the RNA label just easier because it is less noisy?**
That is one explanation and it is compatible with ours. Two checks argue it is not only noise: (1) on discordant
patients, a model trained only on concordant patients sides with RNA (34 : 20 FFPE, 39 : 13 frozen), and (2) the calls
follow the PAM50 subtype (Basal → ER−, Luminal A → ER+), i.e. a biological signal, not random label noise.

## About the methods

**Q: How do you avoid site leakage?**
All CV is grouped by tissue-source site (40 hospitals) using one fixed split file. Without it, HER2 0 vs 1+ looks
learnable (AUC 0.63); with it, it is at chance. Embeddings identify the hospital 82% of the time, so this matters.

**Q: Why frozen encoders and simple heads?**
The question is about what the pretrained representation supports, not about squeezing AUC; frozen encoders make the
comparison between labels clean and cheap. Fine-tuning (LoRA) was considered and dropped for time.

**Q: Your RNA "positive" threshold is arbitrary.**
Yes. We used a 2-component Gaussian mixture and repeated every RNA claim with prevalence-matched and confident-only
cut-offs. ER holds under all three; PR's apparent "IHC preference" did not, so we report PR as no preference.

**Q: Why is the DeLong test valid with two different labels?**
We use DeLong's structural components per patient; the covariance comes from per-patient influence terms, which reduces
exactly to the classic DeLong test when the labels coincide (verified in code) and matches the paired bootstrap within
0.005 in every comparison.

**Q: How do you know your frozen-slide features are comparable to the TRACE DX features?**
We re-extracted 46 DX slides with our pipeline (IDC copy → TRIDENT → UNI): cosine 1.0000 with TRACE (0.33 between
different patients), max per-feature difference ≤ 0.002, identical ER calls 46/46.

## About the results

**Q: One of your five models is not significant for ER.**
Correct: UNI ABMIL gives +0.018 (CI −0.008 to +0.044) on the fixed splits (it was +0.034 on older folds). The direction is
the same in all five; four of five exclude zero; our pre-specified primary model (UNI mean-pool) gives +0.036,
Holm-adjusted p = 0.042.

**Q: Is the 1b result significant?**
Not after correction on FFPE slides: the patient-level consensus of five models is 34 : 20 (p = 0.076, adjusted 0.15).
We call it suggestive. On frozen slides the same question gives 39 : 13 (p = 0.0004).

**Q: How many tests did you run? Couldn't some results be chance?**
About 80. At α = 0.05 we expect about 4 false positives, so only the three pre-specified ER claims are confirmatory
(Holm). Single borderline results (TITAN PR prevalence-matched p = 0.032, ER-RPPA interaction in 1c p = 0.045) are
explicitly labelled not robust.

**Q: What does 1c actually show?**
That the RNA preference is not explained by the FFPE block being a different piece of tumour from the sequenced one: the
frozen section next to the sequenced tissue shows exactly the same advantage (+0.032 vs +0.032; interaction 0.000,
CI ±0.024). A difference larger than 0.024 AUC is ruled out.

**Q: 1e — is ~10% flipped calls a lot?**
For two sections of the same tumour a few millimetres apart, yes: it is the floor on how reproducible a single-slide
binary call can be. Agreement of the probabilities is high (ICC ≈ 0.88); flips concentrate near the decision threshold.

**Q: What about HER2-low?**
At chance from H&E in every model (0.42-0.54), while ERBB2 mRNA alone gives 0.65-0.69 on the same patients. HER2-low is
not a visible morphology for these models.

**Q: Where does the attention model look?**
On tumour, in every slide we reviewed (fat and glass get the least). ER-low Basal tumours attract attention to
high-grade pleomorphic cells like the ER-negative reference; the HER2-enriched exceptions to cohesive or apocrine-like
cells. This is a qualitative, non-pathologist reading of 6 slides.

## About the extensions (Nov 2026; `docs/FINAL_RESULTS.md` 1.1)

**Q: Do you even need a pathology foundation model?**
Yes. An ImageNet ResNet-50 is 0.11-0.18 AUC below UNI on every receptor task and at chance for HER2 (436-patient
subset). But bigger or newer pathology models do not help for ER: none of 11 encoders beats UNI v1 (ceiling ≈ 0.87
with mean pooling; H-optimus-0 ABMIL 0.881 vs UNI ABMIL 0.891).

**Q: Why ABMIL and not a more modern MIL method?**
We tested CLAM-SB, TransMIL, DSMIL and ACMIL on UNI and repeated the best ones on CONCH: none beats ABMIL, and a linear
probe on mean + standard deviation of patch features matches all of them (UNI 0.893 vs 0.891). Caveat: paper defaults,
no tuning. The label preference is the same for every method.

**Q: Can you improve the model?**
The one thing that helped was better labels: training on soft labels from a latent-class model of IHC, RNA and protein
improved ER AUC by 0.015-0.028 in three encoders, even when scored against IHC. Fusion, multi-magnification, multi-task
training, DANN and noise-robust losses did not.

**Q: Is the RNA preference just hospital batch effects?**
No. Standardising each hospital's features with its own statistics removes hospital identity (from 82-96 % to chance)
and the ER preference is still there in all four encoders tested.

**Q: What is the "latent truth"?**
A Dawid–Skene model that treats IHC, RNA and protein as three imperfect tests and estimates each test's error rates and
each patient's probability of being truly ER+. RNA comes out most accurate (specificity 0.98 vs 0.89 for IHC); every H&E
model scores 0.03-0.05 AUC higher against the latent truth than against IHC, but is still less accurate than IHC itself.
Caveat: two of the three tests are molecular and their errors may be correlated.

**Q: Why does the leaderboard result matter?**
Because published benchmarks use one label. With 17 models on the same patients, the top-ranked model changes in 57-91 %
of bootstrap samples when IHC is replaced by RNA or protein. Benchmarks should state, and ideally vary, their reference
standard.

**Q: How do you know when not to trust a prediction?**
Split-conformal prediction gives sets with 90 % coverage overall (ER 0.900, PR 0.909), but coverage varies by hospital
(0.68-0.98). The abstained cases are the ambiguous ones: enriched for frozen-section call flips (odds ratio 4.5) and for
IHC/RNA discordance (odds ratio 2.5).

**Q: Does an FFPE model work on frozen sections?**
Partly: ER AUC falls by ~0.05 and calibration breaks (FFPE-calibrated sets cover only 80 % of frozen cases). The RNA
preference survives. A frozen-section tool would need its own calibration.

**Q: Are the heatmaps just pretty pictures?**
Each map is paired with a number. The linear contribution map is exact (it sums to the model's logit), the four
attention maps agree with it (signed Spearman 0.41-0.50, positive in 13 of 14 patients), 86 % of attention is on
tumour-concept patches, and deleting the top 5 % of patches moves the logit 700 times more than deleting random ones.

**Q: Is ReceptorCheck a clinical tool?**
No: a research-only demo, evaluated retrospectively on TCGA with out-of-hospital predictions. For ER its "flag for
retest" rule finds IHC/RNA-discordant tumours with precision 0.46 against a 5.5 % base rate; abstaining on the 20 % least
confident cases raises ER accuracy from 0.88 to 0.93. HER2 failed its acceptance checks and is not offered.

## Limitations to state before being asked

One cohort (external CPTAC check prepared, on hold); small discordant groups (54 ER, 10 ER-low); RNA cut-offs are not
clinically validated; frozen encoders only (no fine-tuning); some IDC slides (JPEG-2000) are unreadable by our pipeline;
attention review is qualitative; ResNet-50 compared on a 436-patient subset; MIL methods untuned; the latent-class model
assumes independent test errors; the application is retrospective and research-only.
