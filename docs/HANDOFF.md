# Handoff: full context up to 12 Oct 2026

Written so either of us (and our Claude Code sessions) has the full context from the planning sessions.
There are no fixed owners: whoever is working picks the next open item in `docs/PROJECT_PLAN.md`.
Update the "Current state" section whenever you finish something.

---

## 1. How we got here

1. **Original pitches (individual, same base project: Kamineni).**
   - **B1 (Uniyal):** WSI + RNA-seq cross-attention fusion with a pathology foundation model for ER/PR/HER2/PAM50,
     TCGA-BRCA → CPTAC-BRCA. [`archive/pitch_B1_fusion_subtyping_uniyal.pdf`](archive/pitch_B1_fusion_subtyping_uniyal.pdf)
   - **B2 (Srivastava):** tissue graph + pathway graph GNN survival model with modality dropout and
     missing-modality stress tests. [`archive/pitch_B2_graph_survival_srivastava.pdf`](archive/pitch_B2_graph_survival_srivastava.pdf)
2. **Novelty audit (11 Oct).** Both were found **not novel**. Details and citations: [`NOVELTY_RESEARCH.md`](NOVELTY_RESEARCH.md).
   - B1 already done: Fernandez-Romero 2026 (13 FMs, TCGA→CPTAC), Neidlinger 2025 (19 FMs, TCGA→CPTAC),
     Hezil 2026 (CTransPath + RNA cross-attention PAM50). B1 also has **label leakage**: PAM50 is computed from RNA.
   - B2 already done: HGCN (TMI 2023) and SELECTOR (2024) = WSI graph + gene graph + clinical + masked autoencoder.
     Modality dropout: Cheerla 2019, MultiSurv 2021. **CPTAC-BRCA has no survival data**, so the TCGA→CPTAC
     survival test is impossible. "HGMSurvNet" (cited in B2) does not appear to exist.
   - **Three-modality idea** (MRI + WSI + genomics on TCGA-BRCA): only 136 patients have all three, with
     3 OS / 6 PFI events. Not viable.
3. **Course constraint:** novelty must be a *new question, evaluation or finding*, not a new architecture.
4. **Options considered** (ranked in NOVELTY_RESEARCH.md Part 4): (1) which truth does the slide see,
   (2) strict re-evaluation of TCGA-BRCA survival claims, (3) resection → core biopsy, (4) MRI beyond subtype
   for pCR, (5) site-leakage audit.
5. **Decision (11 Oct): Option 1 only**, all five experiments (1a–1e). Fine-tuning experiment "1f"
   (LoRA on UNI) was discussed and **deferred**.

## 2. Key decisions and why

| Decision | Why |
|---|---|
| Option 1 "Which truth does the slide see?" | Highest confidence of being unexplored; runs on public data; informative either way |
| Main encoder **UNI v1**, robustness encoder **CONCH v1**, quick baselines **TITAN** | UNI has precomputed TCGA-BRCA features (TRACE) *and* downloadable weights, so frozen slides can be encoded with the same model; not pretrained on TCGA; fits 11 GB GPUs |
| Avoid Phikon / CTransPath / Kaiko as main encoders | Their pretraining included TCGA slides (leakage) |
| **Site-grouped CV everywhere** (group = tissue source site, barcode chars 6–7) | Embeddings identify the site 82% of the time; HER2-low signal is site leakage |
| Strict HER2 label (IHC 3+ or FISH+ = pos; 0/1+ or 2+/FISH− = neg) | TCGA `her2_status_by_ihc = Positive` includes 11 IHC 1+ and 27 IHC 2+ cases |
| RNA labels = 2-component Gaussian mixture on log2(RSEM+1) | Data-driven cut-off (ESR1 10.17, PGR 8.58, ERBB2 15.14); sensitivity analysis still to do |
| Encoder stays **frozen**; we train heads (ABMIL, logistic regression) | Standard for this kind of study; compute fits Ada |
| Use TRACE precomputed features for DX slides; extract frozen slides ourselves | TRACE covers DX only (verified); avoids downloading 1.16 TB of DX slides |
| Do **not** re-extract all DX slides for 1c | Instead check our extraction matches TRACE on ~50 DX slides |

## 3. Current state (12 Oct 2026)

**Done**
- [x] Metadata pipeline: `scripts/fetch_metadata.py`, `scripts/build_labels.py` (reproduced from a clean run)
- [x] Dataset study and feasibility check for all five experiments → [`DATASET.md`](DATASET.md)
- [x] Baselines B0–B6 on TITAN slide embeddings → [`BASELINE_RESULTS.md`](BASELINE_RESULTS.md)
- [x] Project-update deck (5 Beamer slides) → [`../slides/project_update_slides.tex`](../slides/project_update_slides.tex)

**Not started**
- [ ] UNI / CONCH weight access on Hugging Face (both of us, with IIIT e-mails; Gmail is refused)
- [ ] Download TRACE UNI v1 20× features (62 GB) to Ada and regroup into per-slide bags
- [ ] ABMIL training code
- [ ] Frozen-slide (TS/BS) download + TRIDENT/UNI extraction on Ada (1c, 1e)

## 4. Next actions (see PROJECT_PLAN.md for the full week plan)

No fixed owners. The authoritative to-do list, with checkboxes and the order to work in, is **`docs/PROJECT_PLAN.md`**
(tiers 0–7). Before starting an item, mark it `[~] in progress (Name, date)` and push; see section 0.6 of the plan.

## 5. Gotchas learned the hard way

- **GDC downloads can truncate silently.** The first clinical-table download had 771 of 1,097 patients.
  `fetch_metadata.py` now checks Content-Length; `build_labels.py` asserts more than 1,000 rows.
- Clinical biotab has **3 header rows**: read with `skiprows=[1, 2]`.
- Use primary-tumour samples only: barcode sample code `01`. Sample `11` = normal tissue (many TS slides are `11A-…`).
- TS and BS slides from the same `vial + portion` (e.g. `01A-01-TS1` and `01A-01-BS1`) bracket the sequenced tissue.
- TRACE parquet rows are **chunks of up to 256 patches**, not slides: group by `slide_id` to build bags.
  TRACE UNI v1 = 256 px patches at 20× (and 10×, 5×); CONCH v1 = 512 px patches.
- GTX 1080 Ti (Pascal) has very slow fp16: use fp32 there. RTX 2080 Ti: fp16 autocast. No bf16 on either.
- The shuffled-label control gives AUC slightly below 0.5 (0.43–0.46) with single-seed CV: a known
  fold-prevalence artefact. Repeated CV should bring it to 0.50.
- On Windows, set `LOKY_MAX_CPU_COUNT=4` to silence a harmless joblib warning.

## 6. Where else context lives

- A Claude Docs page with the same content (Uniyal's account, may not be shared):
  https://claude.ai/code/artifact/5564474f-a2e2-429d-bd38-979832982984. **This repo is the source of truth.**
- Course submission requirements for the update: final proposal (clinical problem, input dataset, rough pipeline),
  1–2 baseline experiments with analysis, and the plan for the next week.
