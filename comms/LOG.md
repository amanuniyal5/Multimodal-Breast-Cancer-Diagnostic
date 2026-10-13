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
