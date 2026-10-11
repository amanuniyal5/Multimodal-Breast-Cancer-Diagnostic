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
