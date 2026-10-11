# Dataset study (TCGA-BRCA)

All numbers were produced by `scripts/fetch_metadata.py` + `scripts/build_labels.py` on 11–12 Oct 2026 from
under 3 MB of metadata, without opening a slide.

## Sources

| Source | What it gives | Size | Access |
|---|---|---|---|
| GDC files API (`data_type = Slide Image`) | Every slide: file name, size, DX vs TS/BS, sample/portion codes | 3,112 rows | open |
| GDC BCR Biotab `nationwidechildrens.org_clinical_patient_brca.txt` (file id `8162d394-8b64-4da2-9f5b-d164c54b9608`) | IHC ER/PR/HER2, ER/PR % positive, HER2 IHC score, FISH, histology, race, age, stage | 1,097 patients, 1.7 MB | open |
| GDC Biotab `biospecimen_slide_brca.txt` (`c10090b8-fa42-4e32-86a0-9dedceb5b605`) | Per frozen slide: % tumour, stroma, necrosis, lymphocytes | 1,579 slides | open |
| cBioPortal `brca_tcga_pan_can_atlas_2018` | PAM50; ESR1/PGR/ERBB2 RSEM; ERBB2 GISTIC; RPPA protein | 1,082 tumours | open |
| HF `IUCompPath/TRACE` | Precomputed features for **DX slides only** (UNI v1, UNI2, CONCH, GigaPath, CTransPath, CHIEF, TITAN) | UNI v1 20× = 62 GB | not gated; licence not stated |

**Slide name anatomy:** `TCGA-B6-A0RG-01A-01-TS1` = patient `B6-A0RG` · sample `01` (primary tumour) · vial `A` ·
portion `01` · top frozen section. `-01Z-00-DX1` = FFPE diagnostic slide. Sample `11` = normal tissue (exclude).

## Cohort

1,097 patients, 40 tissue source sites (median 7 patients per site; largest `BH` = 150, `A2` = 102, `E2` = 90).
Median age 58. Histology: 784 ductal, 203 lobular, 46 other, 30 mixed. Race: 757 White, 183 Black, 61 Asian, 92 not available.

## Labels

| Label | Source / definition | Positive | Negative |
|---|---|---|---|
| ER | IHC | 808 | 238 |
| ER | ESR1 RNA (GMM cut-off log2 10.17) | 816 | 265 |
| PR | IHC | 699 | 344 |
| PR | PGR RNA (cut-off log2 8.58) | 614 | 467 |
| HER2 | Strict clinical (IHC 3+ or FISH+ / 0, 1+ or 2+ FISH−) | 145 | 468 |
| HER2 | ERBB2 GISTIC amplified | 123 | 946 |
| HER2 IHC score | 0 / 1+ / 2+ / 3+ | 61 / 272 / 199 / 90 | |
| PAM50 (n = 981) | LumA 499 · LumB 197 · Basal 171 · Her2 78 · Normal 36 | | |
| RPPA protein | ER-α and HER2, continuous | 875 patients | |

Label-only agreement (no images): ESR1 mRNA → ER IHC AUC 0.958; ERBB2 mRNA → strict HER2 AUC 0.873;
ERBB2 mRNA, HER2 IHC 0 vs 1+ AUC 0.652.

## How often the truths disagree

| Pair | Both known | Discordant | Direction |
|---|---|---|---|
| ER IHC vs ESR1 RNA | 1,031 | 57 (5.5%) | 37 IHC+/RNA−, 20 IHC−/RNA+ |
| PR IHC vs PGR RNA | 1,028 | 172 (16.7%) | 138 IHC+/RNA−, 34 IHC−/RNA+ |
| HER2 strict vs ERBB2 amplification | 591 | 64 (10.8%) | 55 clinical+/not amplified, 9 clinical−/amplified |
| ER IHC+ but PAM50 Basal | | 20 | |
| HER2 strict− but PAM50 Her2-enriched | | 13 | |

- ER IHC+/RNA− cases are mostly **ER-low**: of 26 with a known %, 22 have < 40% stained cells (10 of the 17
  patients with < 10% are RNA-negative).
- TCGA `her2_status_by_ihc = Positive` includes 11 IHC 1+ and 27 IHC 2+ cases: never use that field alone.

## Slides (primary tumour only)

| Type | Slides | Patients | Raw size | Features |
|---|---|---|---|---|
| DX (FFPE diagnostic) | 1,131 | 1,060 | 1.16 TB | TRACE (precomputed) |
| TS (frozen, top) | 1,045 | 1,043 | 287 GB | extract ourselves |
| BS (frozen, bottom) | 472 | 427 | 143 GB | extract ourselves |
| TS + BS from the same portion | 374 pairs | 374 | | for 1e |

Frozen slides are tumour-rich (biospecimen table): median 80% tumour cells on TS, 75% on BS.
TRACE TITAN covers 1,125 DX slides from 1,055 patients.

## Feasibility per experiment

| Exp | Analysis set | Smallest class | Still needed | Verdict |
|---|---|---|---|---|
| 1a | ER 989 · PR 986 · HER2 568 paired | HER2+ 139 | nothing | feasible now |
| 1b | Discordant ER 57 · PR 172 · HER2 64 | ER 57 | nothing | feasible now |
| 1c | DX + TS with labels: ER-IHC 965 · ER-RNA 997 · HER2 559 | HER2+ ~130 | TS features | feasible, week 2–3 |
| 1d | HER2 0 vs 1+ with DX: 57 vs 261 | 57 | nothing | feasible but small; add "0 vs low" |
| 1e | TS + BS pairs 374 (ER-IHC 368, HER2 247) | | TS + BS features | feasible, week 3–4 |

## Example patients

| Patient | IHC | RNA | PAM50 | Slides | Why interesting |
|---|---|---|---|---|---|
| TCGA-AR-A1AH | ER+ (10–19%) | ESR1 log2 2.8 | Basal | DX1, TS1 | IHC ER+, molecularly triple-negative-like |
| TCGA-A2-A3XV | ER+ (< 10%) | ESR1 log2 2.3 | Her2 | DX1, TSB | ER-low: guidelines now treat as ER-negative |
| TCGA-B6-A0RG | ER− | ESR1 log2 14.0 | LumA | DX1, TS1, BS1 | Possible IHC false negative; has both frozen slides |
| TCGA-A2-A0EQ | HER2 3+, FISH− | ERBB2 log2 12.6 | Her2 | | Protein and DNA disagree |
