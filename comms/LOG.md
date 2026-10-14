# Team log (append only, newest at the bottom)

## 2026-10-12 · Uniyal (via Claude Code)

**Did:**
- Ran a novelty audit of our two pitches (B1 fusion subtyping, B2 graph survival). Both are already published;
  details in `docs/NOVELTY_RESEARCH.md`.
- We picked **Option 1: "Which truth does the slide see?"** with five experiments 1a–1e (`docs/PROJECT_PLAN.md`).
- Built the metadata pipeline, the dataset study (`docs/DATASET.md`) and baselines on TITAN slide embeddings
  (`docs/BASELINE_RESULTS.md`).
- Wrote the 5-slide project update in Beamer (`slides/project_update_slides.tex`).
- Set up this repo, with `CLAUDE.md` rules and this log.

**Results (preliminary, site-grouped CV):**
- 1a ER: RNA label AUC 0.930 vs IHC label 0.882, Δ +0.048 (95% CI +0.023, +0.074), n = 989.
- 1b ER: model sides with RNA in 39/54 discordant cases (p = 0.0015).
- HER2 0 vs 1+ falls from 0.63 (random CV) to 0.43 (site-grouped): site leakage.

**Blocked on:** UNI / CONCH access on Hugging Face (request with IIIT e-mail).

**Next for me:** download TRACE UNI v1 features to Ada, build per-slide bags, write ABMIL, run 1a with ABMIL.

**Ask for Srivastava:**
1. Read `docs/HANDOFF.md`, then confirm or change the proposed split (you: frozen-slide pipeline + 1c + 1e +
   shared stats code) in a reply entry here.
2. Request UNI and CONCH access on Hugging Face with your IIIT e-mail.
3. Try one frozen TS slide through TRIDENT + UNI on an Ada 2080 Ti node and log the seconds per slide.

## 2026-10-12 · Uniyal (via Claude Code)

**Did:** connected Claude Code to Ada (key + Windows ssh-agent), cloned this repo to `~/project` on Ada,
ran the smoke test `jobs/check_env.sbatch` (job 17358, gnode043). Wrote `docs/ADA.md`.

**Results:** Ada changed since the official guide: SSH via `ada-gw1/ada-gw2.iiit.ac.in`; GPU partition is `u22`
(Ubuntu 22.04), RTX 2080 Ti via `--constraint=2080ti`; modules `u22/cuda/12.x`, `u22/python/3.12.x`;
**compute nodes have internet** (GDC, Hugging Face, PyPI); `/share1` is gateway-only.

**Blocked on:** my Ada GPU quota (4 GPUs / 40 CPUs per user) is mostly used by my other projects right now.

**Ask for Srivastava:** your Ada account has its own 40 CPU / 4 GPU quota, so the frozen-slide extraction
(1c/1e) is probably best run under your account. `docs/ADA.md` has the job template and the key-setup steps.

## 2026-10-13 · Uniyal (via Claude Code)

**Did:** all Ada baseline jobs finished (17536 UNI mean+ABMIL, 17537 CONCH mean+ABMIL, 17660/17661 PR
prevalence-matched ABMIL). Ran `scripts/sanity_check.py` and wrote **`docs/baseline_report.md`**.
Cancelled my unrelated DFT jobs and stopped their submit loops on ada-gw1 (so our project has the full quota).

**Results:** ER: RNA label more predictable than IHC in all 5 model variants (+0.03 to +0.06, CI excludes 0) and the
model sides with RNA on discordant cases (34–39 of 54). HER2: sides with amplification. PR: no preference once the RNA
cut-off is prevalence-matched (the earlier "sides with IHC" was a GMM-threshold artefact). HER2-low: chance from H&E.
Sanity checks all pass except one explained item: scikit-learn 1.5.2 (laptop) vs 1.9.1 (Ada) give different
StratifiedGroupKFold splits for the same seed.

**Next:** fixed split files in `data/splits/` + pinned scikit-learn; update slide 4 (PR bullet, Ada CONCH numbers);
frozen-slide pipeline for 1c/1e (needs UNI weights).

**Ask for Srivastava:** read `docs/baseline_report.md`; please take the frozen-slide extraction (1c/1e) on your Ada account.

## 2026-10-13 (evening) · Uniyal (via Claude Code)

**Did:**
- Fixed, library-independent site-grouped splits `data/splits/folds.csv` (`scripts/make_splits.py`); all scripts
  read them via `common.site_folds()`; scikit-learn pinned to 1.9.1.
- Reran every linear baseline on Ada with the fixed splits (job 19737).
- ABMIL attention maps for 54 ER-discordant + 10 reference patients (job 19722, `results/attention/`).
- Rewrote `docs/PROJECT_PLAN.md` as a tiered plan with checklists and guardrails. **Please read it; it is now the plan.**

**Results:** fixed splits reproduce every conclusion (TITAN ER Δ +0.058, CI +0.031 to +0.084). Attention model sides
with RNA in 35/54 discordant and 8/10 ER-low (the 2 exceptions are HER2-enriched); P(ER+) tracks ESR1 (r = 0.52).
Details: `docs/baseline_report.md` Sections 5b and 5c.

**Next:** Tier 2 (DeLong, combined 1b evidence), Tier 3 (thumbnail overlay), Tier 4 frozen slides.

**Ask for Srivastava:** take Tier 4 (frozen slides, 1c/1e) and the DeLong / combined-evidence items in Tier 2.

## 2026-10-13 (night) · Uniyal (via Claude Code)

**Did:** removed all task assignments from the docs. **There are no fixed owners**: either of us does whatever is
needed next. Earlier "Ask for Srivastava" items in this log are withdrawn as assignments.
How it works now (`docs/PROJECT_PLAN.md` section 0.6): `git pull` → pick the first open item in the lowest unfinished
tier → mark it `[~] in progress (Name, date)` and push → tick it `[x]` and log the result when done.

**Next open items:** Tier 0 (Ada + Hugging Face access for both), Tier 2 (DeLong test, combined 1b evidence),
Tier 3 (thumbnail overlay, discordant-case description), Tier 4 (frozen slides for 1c/1e).

## 2026-10-14 (night) · Srivastava (via Claude Code)

**Did (in progress, results tomorrow):**
- Ada access under `aman.srivastava` works (smoke test job 19986, RTX 2080 Ti). Own venv in `~/venvs/whichtruth`
  (a fresh install took ~1 h on the NFS home; `jobs/setup_env.sbatch` now has 3 h, pip retries, skips torch if present).
- Tier 2: `scripts/stats_delong_1b.py` (DeLong for 1a with different reference labels on the same patients, checked
  against the classic DeLong; combined 1b evidence at patient level across variants) and `scripts/titan_oof.py`
  (TITAN out-of-fold predictions on the fixed splits). Ada job 19989. TITAN reproduces the 1b counts exactly (ER 39:15).
- Tier 3: `scripts/describe_discordant.py` (PAM50 / ER % / histology of the 54 ER-discordant patients vs model calls).
  First look (4 variants): all 15 Basal IHC+/RNA− called ER−, all 9 LumA IHC−/RNA+ called ER+, Her2-enriched mostly follow IHC.
- Tier 4 prep: `data/metadata/frozen_manifest.csv` (`scripts/make_frozen_manifest.py`: 1,516 TS/BS slides, 429 GB,
  373 TS+BS pairs, 30 batches of ~15 GB); `scripts/extract_frozen.py` + `jobs/extract_frozen.sbatch` (GDC download →
  TRIDENT → per-slide mean + capped bag, raw slide deleted). TRIDENT install job 19991; timing test with a random-weight
  ViT-L (same architecture as UNI) job 20061; DX UNI features rebuilt under my account (job 19993).

**Blocked on:** UNI access on Hugging Face for my account (requested by me, not yet approved / logged in on Ada).

**Note:** the attention-map `.npz` files (coords + weights) are only on Uniyal's Ada account, so the thumbnail overlay
(Tier 3) is easiest from there.

## 2026-10-14 (late night) · Srivastava (via Claude Code)

**Correction to my previous entry:** UNI v1 access for my Hugging Face account was granted tonight (it previously had
only UNI2-h); CONCH access is still pending. Token is on Ada (`~/.cache/huggingface/token`).

**Results (Ada job 19989, fixed splits for TITAN/mean-pool, old folds for ABMIL):** files `results/stats_delong_1a.csv`,
`results/stats_combined_1b.csv`, `results/er_discordant_summary.txt`.
- 1a DeLong, ER RNA − IHC (n = 801): positive in all 5 variants, Δ +0.030 to +0.056, DeLong p ≤ 0.04; DeLong and bootstrap CIs agree.
- 1b combined over 5 variants (one consensus call per patient, mean P over variants): ER 34:20 (p = 0.076; 66% of
  individual calls side with RNA, cluster-bootstrap CI 54–76%); HER2 amplification 42:18 (p = 0.003); PR GMM IHC 116:53;
  PR prevalence-matched 64:69 (p = 0.73). TITAN alone (39:15) is the strongest variant, so the pooled ER evidence is weaker.
- ER-discordant patients: all 15 Basal IHC+/RNA− called ER−; all 9 LumA IHC−/RNA+ called ER+; Her2-enriched IHC+/RNA−
  mostly called ER+ (8/11). The slide follows the intrinsic subtype.

**Ada:** PyPI downloads from compute nodes stall for long periods (setup jobs 19985 timed out, 20011 failed, 19991 stalled).
Overnight chain: `jobs/requeue_frozen.sh` → TRIDENT setup 20078 → UNI test on 3 slides 20079 → frozen extraction array
20080 (6 tasks × 5 batches, 2 at a time). Old chain 20061/20076/20077 is on hold (superseded). DX UNI features rebuild: 19993.

**Next open items:** check 20078–20080 and the seconds per slide; extraction-equivalence check vs TRACE on ~50 DX slides;
write the Tier 2/3 numbers into `docs/baseline_report.md`; Tier 3 thumbnail overlay (needs Uniyal's attention `.npz`).

## 2026-10-14 (afternoon) · Srivastava (via Claude Code) · Tier 4 PAUSED, resume steps below

**Found:** downloading frozen slides from **GDC is not feasible from Ada**: 0.02 MB/s per connection on a compute node,
0.29 MB/s on the gateway, 0.3–1.6 MB/s with 16–32 parallel connections (429 GB → weeks). The **NCI Imaging Data
Commons (IDC)** has every TCGA-BRCA slide as DICOM-WSI (readable by OpenSlide 4.0) in public cloud buckets; from a
compute node **Google Cloud Storage gives 1.8 MB/s per connection, 10 MB/s with 16** (gateway: 15 MB/s single).
All 1,516 manifest slides are in IDC (458 GB). AWS S3 copy of IDC is as slow as GDC: use `storage.googleapis.com`.
Hugging Face `chadvanderbilt/TCGA_foundation_model_features` (UNI BRCA, 53.7 GB) looks DX-only (first 82 entries all DX).

**Also fixed today:** TRIDENT's install pulled torch 2.14+cu130 (too new for the Ada driver) and CUDA-13 nvidia
libraries that overwrote cuDNN 9.1. torch is back to 2.5.1+cu121; cuDNN/CUDA come from Ada modules
(`u22/cuda/12.1 u22/cudnn/9.1.0-cuda-12`; `extract_frozen.py` re-execs with them on `LD_LIBRARY_PATH`). Verified conv
fp32/fp16 + UNI load on a 2080 Ti. TRIDENT DataLoaders capped to the job's CPUs.

**State of Ada jobs when paused (Srivastava's account):**
- 20507 `uni_test` RUNNING but useless (old single-connection GDC code, stuck on slide 1). Cancel approved, but the
  `scancel` did not go through (VPN dropped). **First thing on resume: `scancel 20507 20080`.**
- 20080 `extract_frozen` array (6 × 5 batches) PENDING on 20507; would run the GDC path from the current Ada checkout. Cancel it too.
- 20710 `list_hf_archive` RUNNING (lists the HF archive into `results/hf_archive_listing.txt`; harmless; check TS/BS counts at the end of `logs/20710.out`).
- 20061, 20076, 20077 held (superseded); can be cancelled.
- Done: 19993 DX UNI features in `~/project/data/features/uni_v1_20x_*` (my account).

**Code state:** `scripts/extract_frozen.py` now has `--source idc` (default) and parallel ranged downloads with prefetch;
**the IDC path is not yet tested**. `scripts/make_idc_objects.py` (needs `pip install idc-index`, run on the laptop)
writes `data/metadata/idc_slides.csv`; it was interrupted, so that CSV does not exist yet.

**Resume steps (Tier 4):**
1. `ssh ada "scancel 20507 20080"` (and optionally `scancel 20061 20076 20077`).
2. Laptop: `pip install idc-index` (in a short path venv; Windows long-path limit), `python scripts/make_idc_objects.py`,
   commit `data/metadata/idc_slides.csv`, push; `ssh ada "cd ~/project && git pull"`.
3. Test 3 slides from IDC: `sbatch -J uni_test --export=ALL,ENCODER=uni_v1,LIMIT=3,BATCH=0 jobs/extract_frozen.sbatch`;
   check `logs/uni_test_<id>_*.out` for "N patches … download/segment/patch/encode s" and that TRIDENT reads the
   DICOM mpp (column `mpp` in `data/features/frozen_uni_v1/timing_<id>.csv`).
4. If OK: `sbatch --array=0-5%2 --dependency=afterok:<test id> jobs/extract_frozen.sbatch` (QoS: max 20 submitted jobs).
5. Extraction-equivalence check: ~50 DX slides from IDC through the same pipeline vs TRACE features (cosine similarity,
   same predictions); IDC rows with `slide_type == DX` are in `idc_slides.csv`. Run with `--slides <DX barcodes>`.
6. Then 1c (DX vs TS × IHC vs RNA, fixed splits) and 1e (373 TS/BS pairs: ICC, κ, % flipped).

## 2026-10-14 (evening) · Srivastava (via Claude Code) · Tier 4 resumed

**Did:**
- Hugging Face archive `chadvanderbilt/TCGA_foundation_model_features` (UNI, BRCA) fully listed (job 20710):
  1,117 DX slides, **0 frozen** → no shortcut for TS/BS features; we extract them ourselves.
- Built `data/metadata/idc_slides.csv` (`scripts/make_idc_objects.py`): all 1,516 manifest TS/BS slides + 1,133 DX
  slides are in IDC (13,003 DICOM files; TS 306 GB, BS 152 GB, DX 1,205 GB).
- `extract_frozen.py` now downloads from IDC via Google Cloud (parallel ranges, prefetch). Cancelled the stuck GDC test 20507.
- **UNI test on 3 TS slides from IDC passed (job 21299, gnode054):** mpp 0.2525 read from DICOM (40× scans, patches
  at 20×); 1,258–4,260 patches per slide; per slide download 21–39 s, segmentation 7–23 s, encoding 10–24 s.
  Feature files OK (1,024-d mean, ≤ 4,096-patch fp16 bag, no NaN).

**Running on Ada:** array 20080 (all 30 batches, 2 tasks at a time; est. 8–12 h); extraction-equivalence check 21300
(50 single-DX-slide patients from IDC through our pipeline vs TRACE: cosine of slide means vs different-patient baseline,
ER P(+) agreement; output `results/equivalence_check.csv`).

**Note:** frozen-slide mean features have a larger norm than TRACE DX means (~24–27 vs ~20); the equivalence check will
show whether that is the pipeline or the tissue (frozen vs FFPE).

**Next open items:** check 20080 / 21300; then 1c (DX vs TS × IHC vs RNA, fixed splits, % tumour cells as covariate)
and 1e (373 TS/BS pairs); Tier 3 thumbnail overlay + heatmap review; ER-low write-up.

## 2026-10-14 (19:15) · Srivastava (via Claude Code) · Tier 4 progress note

- Extraction array 20080: tasks 0 and 1 at ~1 slide/min each, 71 slides done, 0 failures; all 1,516 expected ~06:00–08:00 on 15 Oct.
- IDC codec check (idc-index `sm_instance_index`): **all TS and BS slides are JPEG baseline** (fine). 48 of 1,132 DX slides
  are JPEG 2000 (YBR_FULL_422), which OpenSlide's DICOM reader cannot open: one of the 50 equivalence-check slides
  (TCGA-AO-A1KO-01Z-00-DX1) failed for that reason, so the check will use 49. Only matters if DX ever has to be
  re-extracted ourselves (then take those 48 from GDC).

## 2026-10-14 (20:30) · Srivastava (via Claude Code) · Extraction-equivalence check passed

**Result (Ada job 21300, `scripts/equivalence_check.py`, `results/equivalence_check.csv`):** 50 single-DX-slide patients
re-extracted from IDC with our TRIDENT + UNI v1 pipeline; 46 succeeded (4 JPEG-2000 slides, all site AC, unreadable).
Our slide means are **numerically the same as TRACE**: cosine 1.0000 (min 1.0000) vs median 0.33 for different patients;
max absolute difference per feature 0.00004–0.0018 (fp16 rounding); identical patch counts in 39/46 (others ±1);
ER P(+) from a TRACE-trained probe r = 1.000, same call 46/46. Independently computed (bag mean vs stored mean 0.9997).

**Meaning:** TRACE DX features and our frozen-slide features come from the same pipeline, so 1c (DX vs TS) compares
tissue, not pipelines. The higher feature norm of frozen slides (~24–27 vs ~20) is a tissue difference.

**Running:** extraction array 20080 (frozen slides), expected done ~06:00–08:00 on 15 Oct. Next: 1c and 1e scripts.

## 2026-10-14 (21:00) · Srivastava (via Claude Code) · 1c / 1e scripts ready

- `scripts/frozen_1c.py`: same patients and fixed folds; linear probe on DX (TRACE) vs primary TS (ours) features ×
  IHC vs molecular label (ER RNA, ER RPPA, PR prevalence-matched, HER2 amplification); interaction
  (TS Δ − DX Δ) with paired bootstrap; split by % tumour cells; 1b on TS vs DX for the same discordant patients.
- `scripts/frozen_1e.py`: ER probe trained on primary TS of out-of-fold patients, applied to both slides of each TS+BS
  pair; ICC(2,1) of logit P, r, κ, % flipped, AUC TS vs BS, |ΔP| vs |Δ tumour %|, feature cosine.
- `jobs/frozen_analysis.sbatch` runs both (CPU, < 1 min). Smoke test on the 258 slides extracted so far (job 21482,
  output in `results/partial_test/` on Ada only): runs end to end. **Numbers not meaningful yet** (partial, non-random set).

**To do when extraction 20080 is complete:** `sbatch jobs/frozen_analysis.sbatch`, copy `results/frozen_1c*.csv*`
and `results/frozen_1e*.csv` to the repo, write up 1c/1e.
