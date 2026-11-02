# AI extensions: models, methods, interpretability and a target application

*Written 3 Nov 2026 in response to the course feedback. Proposal only: nothing here changes what is done or the
Week 4 plan (extended statistics; experiments C, D, G). Pick what to run before starting; a design / technical spec
follows once the selection is made.*

## 0. The feedback and where each point is answered

| Feedback | Answered in |
|---|---|
| Different models (UNI, CONCH): how are they performing, feature maps | §1 (current numbers), §2 (what each model is and why it behaves as it does), §4 (feature maps) |
| Regression, e.g. mRNA expression | §3, X3 |
| ResNet baseline; all foundation models | §3, X1 |
| Other methods apart from ABMIL; can we get an improvement | §3, X2 and X4 |
| Interpretability; t-SNE | §4 |
| Do we need pathology foundation models at all? | §3, X1 (ResNet / natural-image models, label-efficiency curves) |
| Focus more on AI than on the science | All of the above: the medical findings become *one evaluation axis* of an AI benchmark, not the headline |
| A proper application | §5 |

**Reframing.** So far the project answers a medical question with fixed models. The AI version asks: *how do the
encoder, the aggregation method, the training target (classification vs. regression, IHC vs. RNA vs. soft labels) and
the way we interpret and trust the output change what an H&E model can do, and what is the best system we can build
from it?* Our existing results (label preference, site leakage, frozen vs FFPE, section variability) stay. They become
the stress tests every new model is put through, and the building blocks of the application in §5.

---

## 1. Where the current models stand (as an AI system)

ER AUC on the same 801 patients, site-grouped CV (`docs/baseline_report.md` §4.1; ABMIL rows from the earlier folds,
fixed-split ΔAUC in `FINAL_RESULTS.md` T1):

| Encoder → head | Encoder size, pretraining | ER (IHC label) | ER (RNA label) | PR (IHC) | HER2 (strict) |
|---|---|---|---|---|---|
| TITAN slide encoder → linear | slide-level transformer on CONCH v1.5 patches | 0.869 | 0.924 | 0.800 | 0.678 |
| UNI → mean-pool + logistic regression | ViT-L, 100M patches, self-supervised (DINOv2) | 0.873 | 0.908 | 0.795 | 0.627 |
| UNI → ABMIL | same | 0.864 | 0.908 | 0.798 | 0.671 |
| CONCH → mean-pool + LR | ViT-B, 1.17M image–caption pairs (CLIP-style) | 0.844 | 0.876 | 0.776 | 0.616 |
| CONCH → ABMIL | same | 0.847 | 0.896 | 0.784 | 0.674 |
| CONCH zero-shot prompts (no training) | – | 0.48–0.52 | 0.46–0.53 | 0.50 | 0.53 |
| 14 CONCH concept scores → LR | 14 named features | 0.824 | 0.860 | – | – |

What this says about the *AI*, not the biology:

1. **The aggregator barely matters.** ABMIL does not beat mean-pooling (0.864 vs 0.873 for UNI). Either the ER signal
   is spread over the whole tumour (global morphology, not a few rare cells), or our ABMIL is under-tuned. Both are
   testable (X2).
2. **The encoder matters more.** UNI beats CONCH by ~0.03 on ER. CONCH is half the size and was trained to match
   captions, so it compresses away detail that captions never mention. Whether bigger or newer encoders keep improving
   is unknown for our tasks (X1).
3. **There is a ceiling.** ER ≈ 0.87–0.92, PR ≈ 0.80, HER2 ≈ 0.62–0.68, HER2-low = chance for every model. Part of the
   ER ceiling is label noise: the same model scores ~0.04 higher against RNA. That points to *training-target* methods
   (X3, X4), not only bigger encoders.
4. **The representations carry the hospital.** Slide embeddings identify the hospital in 81.6% of cases (chance 18.4%),
   and a random-CV HER2-low "signal" (0.63) was pure leakage. Any method must be checked for this (X4c).
5. **Language alignment does not give receptor status for free.** Zero-shot prompts are at chance, but 14 named
   concepts recover most of the signal (−0.05 AUC) with full transparency. That is an interpretability result
   (X5) and a design option for the application.

---

## 2. Interpreting the models: what each one is and what it should be good at

| Model | What it learned from | What its features probably encode | Expected behaviour here | Caveat |
|---|---|---|---|---|
| **ResNet-50 (ImageNet)** | 1.2M natural photos, supervised | Edges, textures, colour blobs; nothing pathology-specific | The classic baseline of pre-2023 papers. Expected clearly worse (~0.75–0.82 ER) but not useless: grade and texture are generic visual cues | Sensitive to stain colour → more hospital signature |
| **DINOv2 ViT-L (natural images)** *(optional)* | 142M natural images, self-supervised | General shape/texture; same recipe as UNI, different data | Separates "the recipe" from "the pathology data": if close to UNI, the data is not what helps | – |
| **UNI v1** (current main) | 100M patches from 100k slides (no TCGA), self-supervised DINOv2 | Cell and tissue morphology at 20x: nuclei, glands, stroma | Best single encoder so far | 256 px tiles: little context |
| **UNI2-h** | 200M+ patches, ViT-H | Same as UNI, more capacity and data | Tests "does scaling help this task?" | TRACE features exist (93 GB) |
| **Prov-GigaPath** | 1.3B patches, 171k slides (Providence), ViT-g | Same family, largest data | Second scaling point; different hospital system | TRACE features exist (93 GB) |
| **CONCH v1** (current replication) | 1.17M image–caption pairs, contrastive | What pathologists *write about*: grade, necrosis, patterns | Weaker classifier, but has a text side: zero-shot and concepts | Pre-projection features need `proj_contrast` |
| **TITAN** | Slide-level transformer over CONCH v1.5 patches, plus reports | Whole-slide context and layout | Strongest on ER-RNA (0.924) with only a linear head | One vector per slide: no patch maps |
| **CTransPath / CHIEF** (CHIEF builds on CTransPath) | Pretraining included TCGA | – | Use **only** as a "leakage" contrast: do they look suspiciously good? | Pretrained on our test slides |

**Interpreting the heads:**
- **Mean-pool + logistic regression** assumes "the slide is the average of its patches". Because it is linear, each
  patch's contribution is exactly `w · x_patch`. That gives an exact per-patch feature map at no extra cost (X5b).
- **ABMIL** assumes "a weighted subset of patches decides". Its attention says where it looked, not which way the
  patch pushed the prediction.
- **Slide encoders (TITAN)** model relations between regions (tumour next to stroma, etc.) but hide which region
  mattered.

---

## 3. Proposed AI experiments (all on the same patients, fixed site-grouped folds, same stats)

Every experiment reports: AUC for ER/PR/HER2 against IHC **and** the molecular label (so the label preference becomes
a robustness check for each new model), the hospital-identification accuracy of the representation, and DX → frozen
transfer where possible.

### X1. Encoder benchmark: do we need pathology foundation models? *(ResNet baseline, "all FMs")*

- **Encoders:**
  - ResNet-50 ImageNet: TRIDENT `resnet50` on DX slides from IDC, the same pipeline as the frozen slides.
  - Optional: DINOv2 ViT-L natural-image.
  - UNI, UNI2-h, Prov-GigaPath, CONCH, TITAN.
  - CTransPath / CHIEF, as the leakage contrast.
- **Head:** the same mean-pool LR and ABMIL for each, so the encoder is the only change.
- **Label-efficiency curves:** train on 25, 50, 100, 200, 400 and all training patients per fold. The real argument
  for foundation models is usually *fewer labels needed*, not the final AUC. If ResNet catches up at full data but
  loses badly at 50–100 patients, that is the answer to "do we need them". This matters for a hospital like Kamineni
  with a few hundred labelled cases.
- **Also measured per encoder:**
  - hospital-ID accuracy (does a stronger encoder carry more site signature?);
  - whether the RNA > IHC gap holds;
  - HER2-low (does any encoder see it?).
- **Cost:**
  - TRACE encoders: stream the parquet files on `/scratch` inside a job and keep only the mean and a ≤ 4,096-patch fp16
    bag. Home is 50 GB, so never store the raw 93 GB.
  - ResNet-50: ~1 day of GPU jobs (download dominates).
  - All heads: < 1 day.
- **Literature:** Neidlinger et al. 2025 (19 FMs) and Campanella et al. 2025 (FM vs ResNet) benchmark accuracy. Our
  additions are label efficiency under site-grouped CV, the effect of reference-standard choice on the leaderboard
  (experiment H) and site signature per encoder.

### X2. Aggregation methods beyond ABMIL

| Method | Idea | Why it could help here |
|---|---|---|
| Mean + std / top-k pooling + LR or small MLP | Add spread and extremes to the average | Cheap; tests whether "average" throws away information |
| **CLAM-SB** | ABMIL + instance-level clustering loss | Better attention on small datasets |
| **TransMIL** | Self-attention between patches (Nyström) | Models patch–patch context (tumour vs stroma layout) |
| **DSMIL** | Critical-instance + attention to it | Strong when a few patches decide |
| **ACMIL / DTFD-MIL** | Multiple attention branches / pseudo-bags | Reduces attention over-concentration and over-fitting |
| kNN / prototype classifier | Compare to class prototypes | Training-free reference; interpretable by example |

- Same encoder (UNI), same folds, 3 seeds. Hyper-parameters chosen on inner folds only.
- Expected outcome: small or no gain over mean-pooling. That is still a finding: ER is a global-morphology task. To be
  reported honestly either way.
- **Cost:** 11 GB GPUs are enough. ~2–4 h per method at the current ABMIL scale.

### X3. Regression on mRNA expression *(feedback item)*

- **Targets:**
  - predict continuous log2 *ESR1*, *PGR*, *ERBB2* (and optionally a panel, e.g. the 50 PAM50 genes) from H&E;
  - heads: ridge regression on means, and ABMIL with an MSE head.
- **Metrics:**
  - Pearson / Spearman r and R² per gene;
  - then **threshold the predicted expression** and score it as a classifier against IHC and RNA labels.
- **Questions:**
  - Does a model trained on the *continuous* signal classify better than one trained on a binary label? It sees more
    information per patient and needs no cut-off.
  - Does it remove the cut-off sensitivity we saw for PR?
- **Bonus:** the predicted expression is a continuous, interpretable score for the application ("predicted ESR1 is at
  the 12th percentile of ER+ tumours").

### X4. Ways to improve performance

| # | Method | Idea | Expected gain / question |
|---|---|---|---|
| a | **Encoder fusion and ensembling** | Concatenate UNI + CONCH + TITAN means, or average model probabilities | Usually +0.01–0.03 AUC; tests whether encoders carry complementary signal (check with CKA, X5) |
| b | **Multi-task learning** | One ABMIL with heads for ER, PR, HER2, PAM50 subtype and *ESR1/PGR/ERBB2* expression | Shared biology (ER and PR are correlated; subtype explains the calls) → better PR and HER2, which have fewer events |
| c | **Site-invariant features** | Per-hospital feature standardisation (ComBat-style) or domain-adversarial training (gradient reversal against hospital ID) | Hospital-ID accuracy should fall; does the AUC on unseen hospitals rise? Directly targets our leakage finding |
| d | **Learning with noisy labels** | Train on soft labels (the latent-class "probability truly ER+" from experiment C), on IHC∧RNA consensus, or with noise-robust losses (generalised cross-entropy, label smoothing, co-teaching) | Turns the label-disagreement finding into a training method; scored against IHC *and* RNA |
| e | **Multi-magnification** | Add 10x patches (TRACE has UNI 10x) to 20x | More context per patch; cheap with existing features |
| f | **Calibration** | Temperature scaling per model on inner folds | Needed for any threshold-based application; ABMIL probabilities are over-confident |
| g | *Parameter-efficient fine-tuning (LoRA on UNI)* | Adapt the encoder itself | Previously dropped (1f): heavy on 11 GB GPUs; only if time is left |

### X5. Interpretability and representation analysis *(feature maps, t-SNE)*

| # | Analysis | What it shows |
|---|---|---|
| a | **t-SNE / UMAP of slide embeddings**, per encoder, coloured by ER (IHC and RNA), PAM50 subtype and hospital; discordant patients highlighted | Where the disagreeing patients sit (inside the ER− cluster?), and how strongly each encoder clusters by hospital. Quantified, not just pictures: kNN label purity, kNN hospital purity, silhouette |
| b | **Patch-level feature maps** | (i) Exact linear-probe contribution map `w · x_patch` for mean-pool models. (ii) ABMIL attention (already have). (iii) CONCH concept maps (tubules, necrosis…) per patch. (iv) PCA-to-RGB map of patch embeddings (first 3 components as colours) to compare how UNI, CONCH and ResNet "see" the same slide |
| c | **Representation similarity (CKA)** between encoders | Are UNI, UNI2, GigaPath, CONCH and ResNet learning the same thing? Explains when fusion (X4a) can help |
| d | **Prototype / retrieval explanations** | For a slide, the most similar training patches and patients. "This looks like these ER− tumours" is an example-based explanation clinicians find natural |
| e | **Concept bottleneck (done) and concept probing per encoder** | How much of each named concept every encoder encodes (linear probe from embeddings to concept score) |
| f | *Sparse autoencoder on UNI patch features* (stretch) | Unsupervised dictionary of interpretable features; which ones drive ER calls and the discordant cases |

### X6. Reliability (extends experiment F, uses 1e)

- **Uncertainty:** deep ensembles (seeds × encoders) and MC-dropout in ABMIL.
- **Conformal prediction:** split conformal with site-grouped calibration folds, giving a prediction set {ER+}, {ER−}
  or {ER+, ER−} with guaranteed coverage.
- **Risk–coverage curves:** how much accuracy is gained by abstaining on 10–20%?
- **Checks against what we already measured:**
  - Do abstentions concentrate on the TS/BS flips from 1e?
  - Do they concentrate on IHC/RNA-discordant tumours?

---

## 4. Feature maps and t-SNE: what the figures will look like

- **One figure per encoder:** a 2-D map of all 1,055 slides.
  - Panel 1: coloured by ER RNA.
  - Panel 2: coloured by hospital.
  - The 54 discordant patients marked as crosses.
  - Expected pattern: ER and subtype as the main axis for UNI, a strong hospital pattern for ResNet.
- **One figure per example slide**, the same slide shown 5 ways:
  - H&E thumbnail;
  - UNI linear contribution map (red = pushes to ER+, blue = ER−);
  - ABMIL attention;
  - CONCH "tubules" and "necrosis" concept maps;
  - PCA-RGB map for UNI vs ResNet.
  - Uses the 10 slides already in `results/attention_overlay/`, including the ER-low Basal and HER2-enriched
    exceptions.
- **One CKA heatmap:** encoder × encoder similarity.

---

## 5. The application: an H&E receptor-status second reader with built-in trust checks

**Working name: *ReceptorCheck*.** A tool that takes a routine H&E slide (FFPE, or a frozen section during surgery)
and returns an ER/PR/HER2 estimate with an explanation and a reliability check. When the lab's IHC result is
available, it flags whether the result is consistent with what the tissue shows.

**Who uses it and why.**
1. **QA second reader in the pathology lab.** IHC is affected by fixation, antibodies and reader thresholds, and we
   showed the model's disagreements with IHC concentrate in cases where IHC itself is unreliable. The tool flags these
   for repeat staining or molecular testing.
2. **Triage in low-resource settings.** H&E is universal, IHC costs time and money. A confident H&E estimate lets a
   lab prioritise; uncertain cases go straight to IHC.
3. **Intra-operative early estimate.** We showed frozen sections carry the same signal as FFPE (1c), so an estimate
   can be available before the FFPE block is processed.

**Pipeline and which experiment justifies each part:**

| Component | What it does | Built from |
|---|---|---|
| 1. Slide QC + tiling | Tissue segmentation, 20x patches (FFPE or frozen) | Existing TRIDENT pipeline; DX ↔ frozen equivalence check |
| 2. Encoder | Patch features | Choice from **X1** (accuracy × label efficiency × site robustness) |
| 3. Aggregator / predictor | ER, PR, HER2 probabilities + predicted *ESR1/PGR/ERBB2* expression | **X2**, **X3**, **X4** (best of fusion, multi-task, noise-robust training) |
| 4. Which "truth" it is trained and evaluated on | Trained on the label that generalises best; reports agreement with IHC and with molecular truth | 1a, 1b, experiments **C** and **E** |
| 5. Reliability layer | Calibrated probability; conformal set; "abstain → do IHC"; warning when the slide looks unlike training hospitals (out-of-distribution score on the embedding) | **X6**, 1e (section variability), site-leakage finding, **X4c**, **G** (DX ↔ frozen transfer) |
| 6. IHC consistency check | If the IHC result disagrees with a confident H&E call, flag it for retest, with a ranked score | Experiment **D** (label auditor with RNA as the oracle), 1b, ER-low / HER2-enriched findings |
| 7. Explanation panel | Heatmaps, concept profile ("high pleomorphism, necrosis, few tubules"), similar reference cases | **X5**, concept bottleneck (Week 3), attention maps |
| 8. Scope limits, stated in the tool | No HER2-low calls (H&E cannot see it); results near the threshold marked uncertain | 1d, 1e |

**How the application itself is evaluated (retrospectively on TCGA, site-grouped):**
- **Triage:** at a fixed abstention rate (say 20%), the accuracy on the cases it keeps, and how many IHC tests would be
  saved.
- **QA:** if the top k% most-flagged IHC results are sent for retest, the fraction that are true IHC/RNA discordances
  (precision at k), compared with random.
- **Shift:** the drop in performance on held-out hospitals and on frozen slides, and whether the out-of-distribution
  warning fires on those cases.

**Demo:** a small local web app (Gradio or Streamlit) running on precomputed features for TCGA patients. Pick a
patient and see the prediction, the reliability check, the IHC consistency flag and the explanation panel. No new
slide processing is needed for the demo.

---

## 6. Suggested priority (time and compute)

| Priority | Items | Why | Rough cost |
|---|---|---|---|
| **Must** | X1 (TRACE encoders + ResNet-50, mean-pool LR + label-efficiency curves); X5a t-SNE/UMAP; X5b feature maps; X3 regression on *ESR1/PGR/ERBB2* | Directly answers the feedback; mostly existing features | 2–3 days (ResNet extraction runs in the background) |
| **Should** | X2 (3 aggregators: TransMIL, CLAM, DSMIL or ACMIL); X4a fusion; X4b multi-task; X4c site-invariance; X5c CKA; demo app | "Can we improve?" plus the application | 3–4 days |
| **Could** | X4d noisy-label training (after experiment C gives soft labels); X6 conformal / selective prediction; X5d retrieval | Strong links to the application's trust layer | 2–3 days |
| **Stretch** | X4g LoRA fine-tuning; X5f sparse autoencoder; DINOv2 natural-image encoder | Heavy or exploratory | – |

C, D, G (Week 4 plan) stay as planned. In this framing they become parts 4–6 of the application.

## 7. Risks and guardrails

- **Storage:** UNI2-h and GigaPath features are ~93 GB each and home is 50 GB. Stream them on node-local `/scratch`
  inside the job and keep only means and capped fp16 bags.
- **Licences and access:** UNI2-h and GigaPath weights are gated; TRACE's licence is not stated; CTransPath and CHIEF
  saw TCGA during pretraining (contrast only, never "best model").
- **Fair comparison:** every method uses the fixed folds, the same patients, inner-fold tuning only, and reports n and
  a CI. Model-selection results are labelled as engineering benchmarks (secondary), outside the Holm family P1–P3.
- **Over-interpretation of maps:** heatmaps and t-SNE are illustrations. Each is paired with a number (purity, CKA,
  per-patch contribution sums) so the claims rest on quantities, not pictures.
