> **Status (12 Oct 2026):** the team chose **Option 1** ("Which truth does the slide see?"). This file is kept as the record of why B1/B2 were dropped and which alternatives were considered. Current plan: `PROJECT_PLAN.md`.

# Novelty research: B1 / B2 audit and evaluation-style project options

*Compiled 2026-10-10 from two research agents. Sources: PubMed, arXiv, PMC, and the GDC, TCIA, cBioPortal, Xena and HF APIs.*

**How much to trust this:** the agents looked up live sources, but I have not re-opened every paper myself. Before pitching, open the papers marked ★ and confirm the numbers. Then run one 30-minute Google Scholar / arXiv check per chosen direction. "Searches empty" means the agents found nothing; it does not prove nobody has done it.

---

## Part 1: Are the existing proposals novel?

### B1 (FM + RNA fusion for ER/PR/HER2/PAM50, TCGA→CPTAC): not novel

| Claimed novelty | Already done by |
|---|---|
| Swap ResNet50 for a pathology FM | ★ Fernandez-Romero et al., *Med Biol Eng Comput* 2026 ([PMC13269319](https://pmc.ncbi.nlm.nih.gov/articles/PMC13269319/)): 13 FMs × 3 MIL methods, TCGA-BRCA → CPTAC-BRCA, PAM50 + ER/PR/HER2. With Virchow2, PAM50 macro-F1 drops from 0.54 to 0.36. ★ Neidlinger…Kather, *Nat Biomed Eng* 2025 ([doi](https://doi.org/10.1038/s41551-025-01516-3)): 19 FMs, TCGA-BRCA → CPTAC-BRCA for ER/PR/HER2. Campanella et al., *Nat Commun* 2025: FM vs ResNet50 for ER/PR/HER2 on about 2,000 Mount Sinai slides. |
| FM features + RNA, cross-attention, PAM50 | ★ Hezil et al., *Med Eng Phys* 2026 ([doi](https://doi.org/10.1088/1873-4030/ae449b)): CTransPath + RNA-seq; gated attention vs cross-attention vs concatenation; RNA-only, WSI-only and fusion baselines; 924 TCGA patients. |
| Public benchmark | Patho-Bench / THREADS (2025): TCGA-BRCA ER/PR/HER2 with BCNB as the external test set. |

**Scientific problem with B1: label leakage.** PAM50 is computed from RNA expression of 50 genes. ER/PR/HER2 status tracks ESR1/PGR/ERBB2 mRNA. An RNA branch therefore nearly hands the model the label, which is why in-distribution AUCs reach 0.96–0.99. Any reviewer who knows the field will raise this.

### B2 (tissue graph + pathway graph + modality dropout, survival): not novel as an architecture

| B2 element | Already done by |
|---|---|
| Spatial patch graph → GNN for survival | Patch-GCN (MICCAI 2021); Pathomic Fusion (TMI 2020) |
| Pathway-structured omics on TCGA-BRCA | ★ SurvPath (CVPR 2024): BRCA C-index 0.655 ± 0.089. ProtoPathway (arXiv 2605.21454, 2026): bipartite gene–pathway GNN fused with WSI, BRCA included. |
| **Tissue graph + gene/pathway graph + clinical, with masking** (B2's headline claim) | ★ **HGCN** (IEEE TMI 2023): WSI graph + genomic graph + clinical, plus an online masked autoencoder for missing modalities. ★ **SELECTOR** (2024): WSI patch graph + gene-set graph + clinical, plus a masked autoencoder. Both use 6 non-BRCA TCGA cohorts. Our pitch called SELECTOR "dropout-only"; that description is wrong, because SELECTOR is the closest prior work to B2. |
| Modality dropout | Cheerla & Gevaert, *Bioinformatics* 2019 (pan-TCGA); MultiSurv, *Sci Rep* 2021 (found that dropout did **not** help) |
| Missing-modality stress tests on BRCA | ★ HEALNet (NeurIPS 2024), DisPro (CVPR 2025), EMMS (2026), MUST (CVPR 2026), M2Surv, ProSurv, AdaMHF |
| TCGA → CPTAC survival test | **Impossible.** CPTAC-BRCA has no survival data: GDC shows vital_status "not reported" for all 134 cases, and cBioPortal has no OS/DFS fields. |
| "HGMSurvNet" citation | Not found anywhere. Likely a wrong citation; remove it. |

---

## Part 2: Has masking been done, and have all three modalities been used together?

**Masking: yes, extensively.** This covers modality dropout (Cheerla 2019, MultiSurv 2021), masked autoencoders that reconstruct a missing modality (HGCN 2023, SELECTOR 2024), and masked modelling across WSI + multi-omics including BRCA survival (MORPHEUS, arXiv 2508.00969, 2025).

**Reading A, three modalities = WSI + omics + clinical: already done.**
- With masking: HGCN and SELECTOR (non-BRCA), Cheerla, MultiSurv.
- Without masking, on TCGA-BRCA: MM-SurvNet, M3FusionNet.
- The only empty cell is "tri-modal + masking evaluated specifically on TCGA-BRCA". That is a change of cohort, not a contribution.

**Reading B, three modalities = MRI + WSI + genomics on TCGA-BRCA: not done, and not feasible.**

| Fact (verified via API) | Value |
|---|---|
| TCGA-BRCA patients with MRI on TCIA | 137 |
| Patients with MRI + diagnostic WSI + RNA-seq | **136** |
| Survival events in those patients | **OS: 3, PFI: 6**, too few for a Cox model or C-index |
| CPTAC-BRCA radiology | None |
| Closest work | MuTriM (*Eur J Cancer* 2026): MRI + WSI for recurrence, trained on 335 Chinese patients and tested on 126 TCGA patients; genomics used only post hoc |

**Verdict:** tri-modal + masking is not a viable novelty. Survival is impossible with n=136 and 3–6 events. Subtyping on 136 multi-vendor MRIs is weak, and MRI-based subtyping on this set has already been published by the TCGA radiogenomics group.

---

## Part 3: Facts that shape every option

| Item | Fact |
|---|---|
| TCGA-BRCA events (TCGA-CDR) | OS 151 / 1,091; PFI 145 / 1,091; DSS 83; DFI 84. That is about 29 events per test fold in 5-fold CV, and published fold SDs (±0.07–0.09) are larger than the gaps between methods. |
| TCGA-BRCA slides | 1,133 FFPE diagnostic slides (1,062 cases, 1.16 TB) plus 1,979 frozen slides (0.48 TB) |
| **Precomputed features** | **IUCompPath/TRACE** on Hugging Face, not gated (licence not stated; check it). BRCA sizes: UNI v1 20× = 62 GB, UNI2 = 93 GB, GigaPath = 93 GB, CTransPath = 47 GB, CONCH = 8 GB, plus TITAN slide embeddings. MahmoodLab/UNI2-h-features is gated and needs an **institutional email**. |
| FM weights | Open or auto-approved: Virchow, Prov-GigaPath, H-optimus-0, Hibou, Kaiko (Apache/MIT). Non-commercial licence: UNI/UNI2/CONCH/Virchow2 (Mahmood-lab models may need an institutional email). |
| Ada GPUs | The **GTX 1080 Ti has very slow fp16, so run fp32 there.** Use the **RTX 2080 Ti for fp16**. Neither supports bf16. ViT-L inference fits easily in 11 GB. |
| Own extraction cost | All 1.16 TB of diagnostic slides would take about 1.5–3 days on one 4×1080 Ti node and fill about 65% of scratch. **Use TRACE instead.** Frozen slides (0.48 TB) and BCNB are feasible to extract yourselves. |

---

## Part 4: Novelty options (evaluation / finding type, no new architecture)

Ranked by (confidence that it is unexplored) × (sure-shot: informative whichever way it comes out) × (fit with compute).

### Option 1: "Which truth does the slide see?" (label source × slide type)

**Confidence it is unexplored: HIGH.** Kamineni fit: pathology + omics, diagnosis. Natural owner: B1.

TCGA-BRCA has three competing ground truths for the same patient: clinical IHC (ER/PR/HER2 plus HER2 IHC 0/1+/2+/3+ and ISH), RNA (ESR1/PGR/ERBB2 expression, ERBB2 copy number), and PAM50. Up to ~20% of cases disagree. TCGA also has two kinds of slide: FFPE diagnostic slides (DX), and frozen slides (TS/BS) cut from the same tissue piece that was sequenced.

- **Closest prior work:** Shamai 2024 (H&E model flags IHC errors, no RNA comparison); ESR1-expression-from-H&E models (RNA targets only); Path2Omics 2025 (frozen vs FFPE for gene expression, but not receptors); Dawood 2026 (labels taken as given).
- **Searches that came up empty:** "IHC vs mRNA discordant H&E deep learning", "RNA-derived vs IHC labels histopathology", "frozen FFPE TCGA receptor transfer".
- **Compute:** DX features come from TRACE. Frozen slides are a 0.48 TB download plus about 20–40 GPU-hours of extraction on 2080 Ti. The heads (ABMIL / logistic regression) train in minutes.

**Standalone inferences** (each one is its own trained model or test, none derived from another):

| ID | Experiment | Why either outcome is a result |
|---|---|---|
| 1a | Train ER (and HER2) models on IHC labels vs RNA labels, using the same patients and folds. Which label is more predictable from H&E? | RNA > IHC points to IHC noise. IHC > RNA means morphology tracks protein. |
| 1b | Train only on concordant cases, then test on IHC≠RNA cases. Does the model side with IHC or RNA? (binomial test) | Either side gives a direct answer to "which truth". |
| 1c | 2×2 design: FFPE-DX vs frozen-TS × IHC vs RNA labels (interaction test) | A positive result means tissue sampling location matters. A null means the label mismatch is not explained by sampling location. |
| 1d | HER2 0 vs 1+ (HER2-low): is H&E or ERBB2 mRNA more separable? | The clinically hot HER2-low question, compared against an independent molecular reference. Pathologists reach consensus on only 26% of these calls. |
| 1e | Within-patient reproducibility: TS vs BS frozen-slide predictions (ICC, κ, share of calls that flip) | No prior prediction-level study. |

**Risks:** RNA cutpoints are arbitrary (run sensitivity analyses); there may be only about 50–100 discordant cases; HER2 IHC fields are sparse in TCGA (verify in week 1); frozen-slide quality is poor.

### Option 2: Strict re-evaluation of multimodal survival on TCGA-BRCA (rescues B2)

**Confidence it is unexplored: MEDIUM (about 55–65%).** Kamineni fit: risk / outcome. Natural owner: B2.

Use existing open models (HEALNet, SurvPath, ABMIL-Cox, plus a linear Cox on FM embeddings) on TRACE features.

- **Closest prior work:** Sucre 2025 (imaging adds nothing beyond clinical + omics for OS, but no CIs); Song 2025 (ML4H); Ghawami 2026 (calibration audit); Gustafsson 2026 (CIs overlap across FMs); SurvBoard (omics only).
- **Compute:** low (features are precomputed and the models are small).

| ID | Experiment | Standalone? |
|---|---|---|
| 2a | Power simulation: what is the minimum detectable ΔC-index for OS, PFI and DSS with 145–151 events? Tabulate published BRCA gains against that threshold. | Yes. Pure statistics plus meta-research, no model training. |
| 2b | Random vs **site-preserved** folds for multimodal survival models: how large is the inflation? | Yes |
| 2c | **Realistic missingness** (whole sites missing a modality, or missingness correlated with stage) vs random masking: do robustness rankings flip? | Yes. Every BRCA missing-modality paper masks at random. |
| 2d | Incremental value of FM-WSI and of RNA over a **full** clinical Cox (age, stage, T, N, PAM50/receptors) for PFI, using a likelihood-ratio test and compareC | Yes |
| 2e | 0–90% missing-rate curves with vs without modality dropout, using modern FM features | Weakest; partly covered by HEALNet / MUST / MultiSurv |

**Risks:** reviewers may say "known in general". The novelty rests on combining power, site and realistic missingness specifically for BRCA, so frame it as an audit.

### Option 3: Resection → core-biopsy realism

**Confidence it is unexplored: MEDIUM-HIGH.** Kamineni fit: strong, since hospitals see biopsies first.

- **Data:** TCGA (surgical resections, TRACE features) and **BCNB** (1,058 core biopsies with ER/PR/HER2, HER2 level, Ki67, grade and lymph-node status; free for non-commercial use).
- **Searches that came up empty:** "trained on TCGA tested on BCNB", "resection to core-needle biopsy receptor transfer", "biomarker performance vs number of patches".
- **Compute:** BCNB extraction about 5–15 GPU-hours.

| ID | Experiment |
|---|---|
| 3a | Zero-shot TCGA → BCNB AUC vs BCNB-internal AUC for ER/PR/HER2 |
| 3b | Tissue dose–response: TCGA AUC vs "virtual core" area (1, 2, 4, … cores' worth of contiguous tissue) |
| 3c | Within BCNB: is model correctness associated with the amount of biopsy tissue? |
| 3d | Does biopsy H&E predict lymph-node metastasis beyond tumour size, age and subtype? (likelihood-ratio test) |

**Risk:** BCNB is single-centre, so scanner and specimen type are confounded. Report this openly.

### Option 4: Does pretreatment MRI add anything beyond receptor subtype for pCR?

**Confidence it is unexplored: MEDIUM-HIGH.** Fit: radiology. This direction belongs to the MNJ project (pCR = pathologic complete response after neoadjuvant chemotherapy), not Kamineni.

- **Prior results:** the MAMA-MIA challenge's best *external* pCR balanced accuracy was 0.558 (chance is 0.50), and the challenge reported **no clinical-only baseline**. BreastDCEDL claims AUC 0.94 in HR+/HER2− without a clinical baseline.
- **Data:** MAMA-MIA (1,506 cases across 4 TCIA cohorts, with masks; Synapse registration needed) and BreastDCEDL (2,070 cases).
- **Compute:** mostly CPU (PyRadiomics on the provided masks, then logistic regression / gradient boosting).
- **Inferences:**
  - Subtype-only clinical model under leave-one-cohort-out validation.
  - Likelihood-ratio test / ΔAUC for MRI over subtype.
  - Separately trained within-subtype MRI models.
  - Can radiomics predict the source cohort, and does random CV inflate AUC compared with leave-one-cohort-out?
  - Re-test the 0.94 claim under cohort-held-out splits.

### Option 5: FM-era site-leakage and pretraining-overlap audit for breast slide-level labels

**Confidence it is unexplored: MEDIUM / LOW-MEDIUM.** Cheapest option: TRACE features only.

- **Prior work:** Howard 2021 (pre-FM); PathoROB 2026 (patch-level only, and it names slide-level analysis as future work); Neidlinger flags pretraining overlap only as a caveat.
- **Inferences:**
  - Random vs site-preserved AUC for ER/PR/HER2/PAM50 across 5+ encoders.
  - A site-only baseline with no image input.
  - Site-classification accuracy per encoder.
  - FMs pretrained on TCGA (Phikon, Kaiko, CTransPath) vs FMs not trained on TCGA: does the ranking from TCGA-internal CV hold on CPTAC?
- **Good fallback** if frozen-slide or BCNB extraction stalls.

### Lower-ranked (partly done, or weak)

- **"Is it just grade?"**: largely covered by Dawood et al., *Nat Biomed Eng* 2026. What is left is testing modern FMs, a likelihood-ratio framing, and tumour purity.
- **Calibration / fairness by race and age**: low statistical power (about 180 Black patients); calibration is partly covered by Ghawami 2026.
- **Attention faithfulness**: concept already covered by HIPPO and xMIL.
- **MRI vs H&E on the same 136 TCGA patients**: underpowered.
- **Perturbation robustness (stain / JPEG / blur)**: saturated.

---

## Part 5: Suggested team shape (2 people, Kamineni base)

**Theme:** *Stricter tests for a multimodal breast assistant: against which truth, on which specimen, beyond which baseline?*

| Person | Owns | Inferences (all standalone) |
|---|---|---|
| B1 owner (pathology / omics) | Option 1, plus Option 3 if time allows | 1a, 1b, 1c, 1d, 1e; 3a, 3b |
| B2 owner (survival / statistics) | Option 2 | 2a, 2b, 2c, 2d |

Both share a single feature pipeline: TRACE UNI features for diagnostic slides, self-extracted features for frozen slides and BCNB, the same ABMIL / logistic / Cox heads, and patient-level, site-aware splits. **Fallback:** if the frozen-slide extraction fails, swap 1c and 1e for Option 5.

**Week-1 checks before committing:**
1. Exact TCGA HER2 IHC score / ISH fields and their counts, and the number of IHC-vs-RNA discordant cases.
2. The TRACE licence, and whether it covers frozen TS/BS slides.
3. One Google Scholar pass per chosen direction.
4. Open the ★ papers and confirm the numbers quoted above.
