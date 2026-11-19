# Poster guide: Do We Need Pathology Foundation Models?

*19 Nov 2026 · Aman Srivastava, Aman Uniyal. Explanation of every panel and figure of the 20 Nov showcase poster
("Do We Need Pathology Foundation Models? Benchmarking 11 H&E encoders for breast-cancer receptor status"), with a
spoken script and likely questions. The poster itself is final and is not changed by this file. Numbers come from
`results/bench/` (details: [`baseline_report.md`](baseline_report.md) sections 5k-5l).*

## The poster in one minute

**The one sentence to remember:** pathology foundation models are needed for reading ER/PR/HER2 from H&E (an ImageNet
ResNet-50 loses 0.11-0.18 AUC), but a bigger or newer one does not help: none of the 11 beats UNI on ER.

Three supporting messages carry the rest of the poster:

1. **Accuracy is the same, the hospital signature is not.** All patch encoders reach about 0.83-0.87 ER AUC, but from
   their features a classifier can tell which hospital a slide came from 90-96 % of the time; TITAN, CONCH and CHIEF
   carry less (82-85 %).
2. **Slide-level TITAN learns most from few labels:** with 25 training patients it matches UNI trained on about 150.
3. **Cheap fixes work better than bigger models:** adding the patch standard deviation to the mean, training on ESR1
   mRNA instead of IHC, and recalibrating the scores.

**How the poster is laid out (read it left to right):**

| Column | Panels | Figures in it |
| --- | --- | --- |
| Left | Problem & context, Our approach | hospital-vs-ER scatter, accuracy vs hospital-ID plot, pipeline |
| Middle | Our solution, How it works, What we learned | ResNet-50 vs UNI vs TITAN bars, mean+SD and mRNA panels, CKA heat map, label-efficiency curves |
| Right | Impact and what's next, Acknowledgment, References | calibration bars |

**Start every conversation at the ResNet-50 bar chart** (top middle): it answers the title question in one look. Then
move to the CKA heat map and the label-efficiency curves.

**Where the numbers come from:** TCGA-BRCA, 1,055 patients with slides from 40 hospitals; every model is trained and
tested on the same patients with the same hospital-grouped folds.

## Panel by panel: what each text box says and why

### Title and subtitle

*"Do We Need Pathology Foundation Models? Benchmarking 11 H&E encoders for breast-cancer receptor status."* The question
comes from the professor's feedback: with so many large pathology models, does any of it matter compared with an
ordinary ImageNet network? The answer has two halves: **yes, a pathology model is needed; no, a bigger one does not
help.**

- **"11 encoders"** = 11 pathology feature sets: UNI, UNI2-h, Prov-GigaPath, H-optimus-0, Hibou-L, CONCH, UNI at 10x
  magnification, TITAN, GigaPath-slide, plus CTransPath and CHIEF. The last two are marked **†** because they were
  pretrained on data that includes TCGA, so they could have seen our slides: we show them only as a contrast and never
  call them best. The ImageNet ResNet-50 is the 12th, non-pathology baseline.
- **Patch vs slide encoders.** Most models describe 256-pixel tiles (patches); TITAN, GigaPath-slide and CHIEF produce
  one vector per slide. The colour code on the poster: blue = patch encoder, green = slide-level or vision-language
  model.

### The problem & context

- **ER/PR/HER2** decide breast-cancer treatment (hormone therapy, anti-HER2 drugs). Clinically they are measured by IHC
  staining.
- **IHC and RNA disagree in 5-17 % of tumours** in TCGA-BRCA. This is the motivation of our main project: when the
  tests disagree, which one does an H&E model follow?
- **H&E is cheap and always available**; IHC needs extra stains and RNA sequencing is rarely done. Reading receptor
  status from H&E would be a cheap screen or second opinion.
- **The gap:** there are now many large pathology foundation models, and nobody tells you which to use for this task,
  or whether you need one at all.

### Our approach

What makes the comparison fair:

- **Same patients, same folds, same labels for every model** (TCGA-BRCA, 1,055 patients). Any AUC difference is due to
  the encoder, not to the data split.
- **Hospital-grouped folds:** all patients from one hospital (tissue source site) are always in the same fold, so a
  model is always tested on hospitals it never saw. Without this the results are inflated, because the features
  recognise the hospital (see the hospital panel).
- **One simple classifier for all:** average the patch features of a slide (mean pooling) and fit a logistic
  regression. Simple heads make the encoders comparable; aggregation methods were tested separately.
- **Then the improvements:** mean + standard deviation pooling, training on mRNA expression, fusing encoders,
  calibration.
- **Then interpretation:** CKA (how similar two encoders' representations are), t-SNE maps and hospital
  identifiability.

### Our solution

*"An evaluation protocol plus a recipe."* The **protocol** is the fair comparison above; the **recipe** is what to do
in practice: use a pathology foundation model (ResNet-50 loses 0.11-0.18 AUC), do not chase the biggest one, add the SD
to mean pooling, and recalibrate. The bar chart under this box is the evidence for the headline.

### How it works

*"Patch encoders are near-identical (CKA 0.87-0.96), and TITAN differs."* The heat map under it shows that the main
patch encoders learn almost the same representation, which is why none of them beats the others. TITAN (a slide-level
model) learns something different (CKA ≤ 0.55 with every other model). The pipeline diagram on the left shows the steps
from slide to prediction.

### What we learned (four bullets)

| Bullet | What it means | Number behind it |
| --- | --- | --- |
| No model beats UNI on ER (IHC) | The ceiling is about 0.87 with mean pooling; differences are within the CIs | UNI 0.870, TITAN 0.873, H-optimus-0 0.869, CHIEF 0.868; lowest GigaPath 0.833 |
| TITAN with 25 patients matches UNI with about 150 | TITAN is far more label-efficient | ER AUC: TITAN 0.833 at n = 25; UNI 0.822 at 100 and 0.843 at 200 |
| Hospital identity is 90-96 % decodable from patch encoders vs 82-85 % for TITAN, CONCH and CHIEF | A linear classifier recognises the hospital from the features (10 largest hospitals, chance 18 %) | UNI 95.6 %, H-optimus-0 96.0 %, TITAN 81.6 %, CONCH 83.3 %, CHIEF 84.6 % |
| Mean + SD pooling and mRNA training give small gains | Two cheap changes, each about +0.01 to +0.035 AUC | Mean + SD: up to +0.033 (GigaPath), +0.029 (UNI 10x); ESR1-mRNA training: H-optimus-0 0.869 → 0.904 |

### Impact and what's next

- **Practical advice:** use a slide-level model (TITAN) when labels are scarce, and always recalibrate the scores
  before showing them as probabilities.
- **Next steps listed on the poster** (other aggregation methods, patch-level feature maps, the ReceptorCheck second
  reader) **have since been done**: see "Results since the poster" below.
- **Limitation: HER2-low.** Distinguishing HER2 0 from 1+ ("HER2-low") is at chance from H&E for almost every model,
  while ERBB2 mRNA alone reaches 0.65. Note: the printed sentence reads "goes here as a limitation"; say it as "one
  limitation is..." (see the script).

### Acknowledgment and references

- IIIT Hyderabad and Kamineni Hospitals (the base hospital project). Data are public (TCGA-BRCA); Kamineni data were not
  used.
- References: TITAN (Ding et al., arXiv 2024) and UNI (Chen et al., *Nature Medicine* 2024), the two models the poster
  leans on most. Full list in [`FINAL_REPORT.md`](FINAL_REPORT.md).

## Figure by figure: how to read each one

### 1. ResNet-50 vs UNI vs TITAN (top middle, under "Our solution")

**Read it:** one group of bars per task, AUC with 95 % CI; orange = ImageNet ResNet-50, blue = UNI, green = TITAN; the
dashed line at 0.5 is chance. **Point:** the orange bar is always lowest.

| AUC (436-patient subset) | ResNet-50 | UNI | TITAN |
| --- | --- | --- | --- |
| ER (IHC) | 0.736 | 0.867 | 0.912 |
| ER (RNA) | 0.714 | 0.885 | 0.956 |
| PR (IHC) | 0.663 | 0.777 | 0.827 |
| HER2 | 0.480 (chance) | 0.663 | 0.669 |

**Caveat to know:** ResNet-50 features had to be extracted by us from the slides, so this comparison uses a pre-set
random subset of 436 patients (n = 240-398 per task), not all 1,055. UNI and TITAN are scored on exactly the same subset.
Source: `results/bench/x1_subset/aucs.csv`.

### 2. "Less hospital, more biology" (top left scatter)

**Read it:** each dot is an encoder. For every patient we look at its 10 nearest neighbours in feature space. The
x-axis is the share of those neighbours from the **same hospital** (chance 7 %); the y-axis is the share with the
**same ER status** (RNA label; chance 63 %). **Good = top left.**

- TITAN is top left: 36 % same-hospital neighbours, 84 % same-ER neighbours.
- CHIEF (42 %, 79 %) and CONCH (39 %, 78 %) are next.
- The patch encoders sit bottom right: 63-73 % same-hospital neighbours, only 68-74 % same-ER.

Meaning: in the raw feature space, patch encoders group slides by hospital first and by biology second. A classifier
can still read ER, but the space itself is organised around the hospital. Source: `results/bench/x5/embedding_metrics.csv`.

### 3. "Accuracy is the same... hospital signature is not" (left, two panels)

**Left panel:** ER (IHC) AUC with 95 % CI per encoder. All 11 lie between 0.833 (GigaPath) and 0.873 (TITAN), with
overlapping CIs of about ±0.03. **Right panel:** % of slides whose hospital is correctly identified from the features
(10 largest hospitals, chance 18 %): patch encoders 90-96 %, TITAN 82 %, CONCH 83 %, CHIEF 85 %. **Point:** encoders
differ much more in how much hospital they carry than in how well they predict ER. Source: `results/bench/x1/aucs.csv`,
`results/bench/x1/site.csv`.

### 4. Pipeline (bottom left)

Five boxes, top to bottom: H&E slide (FFPE or frozen) → tile into 256-pixel patches → foundation model turns each patch
into a feature vector (UNI, TITAN, CONCH...) → aggregate the patches into one slide vector (mean + SD, or attention MIL)
→ calibrated ER / PR / HER2 probability. Use it to explain the method in 15 seconds.

### 5. "A cheap fix: add the SD" and "Train on expression, test on IHC" (bottom left, two small panels)

**Left:** each line is one encoder; ER (IHC) AUC with mean pooling (left end) and with mean + standard deviation of the
patch features (right end). Every line goes up: GigaPath 0.833 → 0.866, UNI 10x 0.866 → 0.895 (the best linear model
overall), H-optimus-0 0.869 → 0.881, CTransPath 0.833 → 0.843, UNI2-h 0.865 → 0.871. Why it helps: the mean says what
the average tissue looks like; the SD adds how varied it is (heterogeneity), which the mean throws away.

**Right:** the same encoder trained on the IHC label (left end) or trained to predict ESR1 mRNA level by regression
(right end), always **scored against IHC**. Training on the continuous mRNA signal is better: H-optimus-0 0.869 →
0.904 (p = 0.001), CHIEF 0.868 → 0.888, GigaPath 0.833 → 0.857, GigaPath-slide 0.844 → 0.866, UNI 0.870 → 0.882.
Why: mRNA is a cleaner, graded signal than a positive/negative IHC call. Source: `results/bench/x3/reg_to_clf.csv`.

### 6. CKA heat map (middle, "Patch encoders learn the same thing; TITAN differs")

**Read it:** CKA measures how similar two encoders' representations of the same patients are (1 = identical up to
rotation and scaling). Yellow = similar, purple = different. The core patch encoders (UNI, UNI 10x, GigaPath,
GigaPath-slide, CTransPath, H-optimus-0) are 0.87-0.96 with each other (UNI vs UNI 10x 0.95, GigaPath vs GigaPath-slide
0.96). Less similar: UNI2-h (0.46-0.83), CONCH (0.43-0.75), Hibou-L (0.33-0.65). TITAN is ≤ 0.55 with every other
model. So "0.87-0.96" on the poster refers to the core patch encoders, not every pair. **Point:** the big patch models
converge on almost the same representation; TITAN is genuinely different, which is why it behaves differently in the
other plots. Source: `results/bench/x5/cka.csv`.

### 7. Label efficiency: "TITAN with 25 patients ≈ UNI with ~150" (bottom middle)

**Read it:** x = number of training patients (25, 50, 100, 200, 400, all), y = ER (IHC) AUC on the full held-out folds,
averaged over 5 random subsets. TITAN is the top line; grey lines are the other encoders.

| ER (IHC) AUC | 25 | 50 | 100 | 200 | 400 | all |
| --- | --- | --- | --- | --- | --- | --- |
| TITAN | 0.833 | 0.841 | 0.859 | 0.874 | 0.873 | 0.885 |
| UNI | 0.721 | 0.756 | 0.822 | 0.843 | 0.860 | 0.878 |
| CONCH | 0.752 | 0.783 | 0.813 | 0.831 | 0.834 | 0.856 |

With 25 patients TITAN (0.833) is between UNI at 100 (0.822) and UNI at 200 (0.843): hence "about 150". With all
patients the lines converge. **Why it matters:** hospitals with few labelled cases should use a slide-level model.
Source: `results/bench/x1/curves.csv`.

### 8. "Raw scores are over-confident" (bottom right)

**Read it:** expected calibration error (ECE: average gap between predicted probability and observed frequency) of the
UNI classifier, grey = raw, blue = after Platt recalibration fitted inside the training folds. ER 0.090 → 0.044, PR
0.105 → 0.047, HER2 0.215 → 0.058, ERBB2 amplification 0.221 → 0.023. AUC is unchanged by recalibration; only the
probabilities become trustworthy. **Point:** never show raw scores to a clinician as probabilities, especially for HER2.
Source: `results/bench/x4/calibration.csv`.

## The script

Speak in plain language; the parts in *[brackets]* are where to point. Numbers are the ones printed on the poster.

### 60-second version (for someone passing by)

> *[Point at the title.]* Pathology now has many big foundation models. We asked a simple question: for reading
> breast-cancer receptor status from an ordinary H&E slide, do we need them, and which one should we use?
>
> *[Point at the bar chart.]* We tested 11 pathology models and a normal ImageNet network on the same 1,055 TCGA
> patients, always testing on hospitals the model never saw. The ImageNet network loses 0.11 to 0.18 AUC, and for
> HER2 it is at chance. So yes, a pathology model is needed.
>
> *[Point at the CKA heat map.]* But a bigger one does not help: none of the 11 beats UNI on ER, because the big patch
> models learn almost the same thing. The slide-level model TITAN is the different one *[point at the curves]*: with
> only 25 labelled patients it matches UNI trained on about 150.
>
> *[Point at the calibration bars.]* Our practical recipe: use a pathology model, use a slide-level one when labels
> are scarce, add the patch standard deviation to the average, and always recalibrate the scores before calling them
> probabilities.

### 3-minute walk-through (poster order, left to right)

**1. Problem (left top, ~25 s).**

> Doctors decide breast-cancer treatment on three markers: ER, PR and HER2. They are measured with special IHC stains,
> and in 5 to 17 percent of tumours IHC disagrees with RNA. H&E slides are cheap and exist for every patient, so
> predicting these markers from H&E is attractive. There are now many pathology foundation models, large networks
> pretrained on millions of tissue images, and nobody tells you which one to use, or whether you need one at all.

**2. Approach (left, ~30 s).** *[Point at the pipeline.]*

> We cut each slide into small tiles, turn each tile into a feature vector with a foundation model, average the
> vectors, and train a simple classifier. Every model gets exactly the same patients, the same labels and the same
> folds. The folds are grouped by hospital, so we always test on hospitals the model has never seen. That matters, as
> I will show in a moment.

**3. Main result (middle top, ~30 s).** *[Point at the bar chart.]*

> Here is the answer to the title. Orange is a standard ImageNet ResNet-50, blue is UNI and green is TITAN. ImageNet
> features lose 0.11 to 0.18 AUC on every task, and for HER2 they are at chance. So pathology-specific pretraining
> really matters.

**4. But bigger does not help (left middle and CKA, ~40 s).** *[Point at the accuracy / hospital panel, then the heat
map.]*

> On the left, all eleven models give almost the same ER accuracy, between 0.83 and 0.87, and the confidence intervals
> overlap. What differs is the hospital signature: from the features of most models a simple classifier can tell which
> hospital a slide came from 90 to 96 percent of the time. TITAN, CONCH and CHIEF carry less, 82 to 85 percent. The
> heat map explains why the accuracy is the same: the main patch encoders learn nearly identical representations, CKA
> 0.87 to 0.96. TITAN, which works at the slide level, is the different one. *[Point at the scatter.]* And it is the
> model whose nearest neighbours share ER status rather than hospital.

**5. Few labels (bottom middle, ~20 s).** *[Point at the curves.]*

> If a hospital has few labelled cases, the choice matters. With only 25 training patients TITAN reaches 0.83, which
> UNI needs about 150 patients to reach. With all patients they converge.

**6. Cheap improvements (bottom left, ~25 s).** *[Point at the two small panels.]*

> Two cheap fixes help more than a bigger model. First, add the standard deviation of the patch features to the
> average: every line goes up, by up to 0.03. Second, train the model to predict ESR1 mRNA instead of the IHC label:
> even when we score it against IHC it does better, for H-optimus-0 from 0.87 to 0.90.

**7. Recalibrate and conclude (right, ~30 s).** *[Point at the calibration bars.]*

> Raw scores are over-confident, especially for HER2; a simple recalibration inside the training folds cuts the
> calibration error by half or more, from 0.22 to 0.06 for HER2. So our recipe is: use a pathology foundation model, do
> not chase the biggest one, use a slide-level model when labels are scarce, add the SD, and recalibrate. One
> limitation: HER2-low, the difference between HER2 0 and 1+, is essentially invisible in H&E for almost every model,
> while ERBB2 mRNA alone reaches 0.65. Since the poster we also tested other aggregation methods, built patch-level maps
> and a second-reader demo; happy to show those.

### Delivery tips

- Lead with the bar chart, not the problem box: people decide in 5 seconds whether to stay.
- Say "tested on hospitals the model never saw" at least once; it is what makes the numbers credible.
- Keep the "which truth" story (IHC vs RNA) for people who ask what the bigger project is (see the questions section).

## Results since the poster (for "what's next?")

The three next steps printed on the poster are done. Mention them briefly; each has one number to quote.

| Poster said "next" | What we found | Number to quote |
| --- | --- | --- |
| Other aggregation methods | CLAM, TransMIL, DSMIL and ACMIL do not beat ABMIL; mean + SD pooling matches all of them, on UNI and on CONCH | ER (IHC): mean + SD 0.893 vs ABMIL 0.891 (UNI); 0.872 vs 0.868 (CONCH) |
| Patch-level feature maps | For 14 patients: exact per-patch contribution maps, 4 attention maps, concept and colour maps (`results/bench/x5/maps/`). Attention and evidence sit on tumour | 86 % of attention on tumour patches that cover 44 % of the tissue; the four attention maps agree (13 of 14 patients) |
| Second-reader application (ReceptorCheck, `app/`) | Research demo: calibrated ER/PR/HER2 with confidence sets, a reliability status and an "IHC may be wrong, retest" flag, all from out-of-hospital predictions | ER flag precision 0.46 vs 5.5 % base rate; HER2 failed its checks and is not offered |

**Two further results, if someone asks what else you found:**

- **Better labels beat bigger models.** Training on soft labels that combine IHC, RNA and protein (a latent-class
  model) raises ER AUC by 0.015-0.028 in all three encoders tested, even when scored against IHC.
- **The benchmark winner depends on the label.** Ranking 17 models against IHC, RNA or protein, the top model changes
  in 57-91 % of bootstrap samples. So a leaderboard of H&E models should say which "truth" it used.

**HER2-low update:** one model, H-optimus-0, does separate HER2 0 from 1+ (AUC 0.63 with mean pooling, 0.65 with
attention MIL, about the same as ERBB2 mRNA at 0.65), but no other encoder does. Report it as a lead, not a finding.

## Likely questions and honest answers

**Q: Why is ResNet-50 compared on fewer patients?**
We had to extract ResNet-50 features ourselves from the slide images, which is slow, so we used a pre-set random subset
of 436 patients. UNI and TITAN are scored on exactly the same subset, so the comparison is fair; the gap (0.11-0.18 AUC)
is far larger than the CIs.

**Q: Isn't UNI just the best model then?**
For ER it ties with the others: TITAN 0.873, UNI 0.870, H-optimus-0 0.869, all within each other's CIs. UNI is our
reference because it is the main encoder of the project. The honest answer is "any modern pathology model, then choose
by labels available and hospital signature".

**Q: What does "hospital identity is decodable" mean, and why care?**
A linear classifier given only the features names the hospital correctly 90-96 % of the time (chance 18 %). Slides
carry the hospital's staining and scanning style. If you test on the same hospitals you trained on, the model can cheat
by learning hospital, so all our tests use hospitals the model never saw.

**Q: What is CKA?**
Centred kernel alignment: a number from 0 to 1 saying how similar two models' representations of the same patients are,
ignoring rotation and scale. 1 means they encode the same information.

**Q: What are the daggers (†)?**
CTransPath and CHIEF were pretrained on data that includes TCGA, the same data we test on, so their numbers could be
optimistic. We show them only as a contrast.

**Q: Why does TITAN learn from so few labels?**
TITAN is a slide-level model: it already summarises the whole slide with a transformer trained on slides and pathology
reports, so the classifier has less to learn. Patch encoders give a raw average that the classifier must sort out with
more examples.

**Q: Why does adding the standard deviation help?**
The mean describes the average tissue; the SD describes how heterogeneous it is (for example, mixed tumour and stroma,
or varied nuclei). The mean alone throws that away.

**Q: Why is training on mRNA better even when scored against IHC?**
mRNA is a continuous, less noisy signal; IHC is a positive/negative call made by eye with known errors. Learning the
graded signal gives a better ranking of patients, including against IHC.

**Q: What is the bigger project behind this poster?**
Our main question is which "truth" an H&E model follows when IHC and RNA disagree. The RNA-defined ER label is more
predictable from H&E than the IHC label (+0.036 AUC, Holm-adjusted p = 0.042), and this holds for all 11 encoders.

**Q: Is this ready for clinical use?**
No. One public cohort (TCGA-BRCA), retrospective, no external hospital test yet. Our ReceptorCheck demo is
research-only.

### Caveats to say before being asked

- One cohort (TCGA-BRCA, 40 hospitals); no external validation yet.
- ResNet-50 compared on a 436-patient subset.
- Frozen encoders only (no fine-tuning).
- HER2-low is essentially invisible in H&E; the one H-optimus-0 signal is a lead from a single model.
- "11 encoders" includes UNI at a second magnification and two models marked † that may have seen TCGA in pretraining.
