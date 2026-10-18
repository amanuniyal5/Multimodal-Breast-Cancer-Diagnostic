# CLAUDE.md: everything a new Claude Code session needs

*Last updated 19 Oct 2026. Read this file fully before doing anything. It is the condensed memory of all earlier
sessions; the repo itself (docs, results, comms/LOG.md) is the source of truth for details.*

Two people work on this project, each with their own Claude Code session: **Aman Uniyal** (repo owner,
GitHub `amanuniyal5`, Ada user `aman.uniyal`) and **Aman Srivastava**. This GitHub repo is the only shared channel and
the team's memory: https://github.com/amanuniyal5/Multimodal-Breast-Cancer-Diagnostic

---

## 1. The project in 60 seconds

- **Course:** Topics in AI for Healthcare (IIIT Hyderabad). Base hospital project: *Kamineni, Multimodal Breast Cancer
  Diagnostic & Outcome Assistant*. Public data only (TCGA-BRCA); Kamineni data never arrived.
- **Title:** *Which Truth Does the Slide See? Auditing what H&E foundation-model classifiers learn when IHC, RNA and
  protein disagree.*
- **Question:** ER/PR/HER2 are measured clinically by IHC, which disagrees with RNA (and DNA/protein) in 5-17% of
  tumours. When the "truths" disagree, which one does an H&E foundation-model classifier follow, and does the answer
  change between the FFPE diagnostic slide (DX) and the frozen slide (TS) cut next to the sequenced tissue?
- **Novelty type (course rule):** a new, stricter *evaluation* of existing models (same patients, slides, folds;
  different reference standards). **No new architecture.** Our original pitches (B1 WSI+RNA fusion subtyping by Uniyal;
  B2 graph survival with modality dropout by Srivastava) were dropped as already published (`docs/NOVELTY_RESEARCH.md`).
- **Five experiments:**
  - **1a** which label is more predictable from H&E (paired ΔAUC, same patients/folds)
  - **1b** on IHC≠molecular discordant patients, which side a concordant-trained model takes (binomial)
  - **1c** FFPE vs frozen slide × IHC vs molecular label (interaction)
  - **1d** is HER2 0 vs 1+ / 0 vs low visible in H&E (vs ERBB2 mRNA alone)
  - **1e** do two frozen sections (TS vs BS) of one tumour agree (ICC, κ, % flipped)

## 2. Status (19 Oct 2026)

- **Plan tiers 0-5 are DONE** (`docs/PROJECT_PLAN.md`): data, baselines, robustness, interpretation, frozen slides
  (1c, 1e), final statistics and tables.
- **Tier 6 (optional):** CPTAC-BRCA external check **on hold**; LoRA fine-tuning (1f) and Kamineni data **dropped**.
- **Tier 7 (deliverables):** final report `docs/FINAL_REPORT.md` (+ `.docx`), `docs/FINAL_RESULTS.md` (claims + tables
  T1-T5), `docs/QA_SHEET.md`, README with result→script→job map: done. **Open:** compile the Beamer decks on Overleaf
  (no LaTeX on laptop or Ada) and rehearse.
- **Latest deck:** `slides/week1_week2_update.tex` (Week 1 frames + Week 2 updates incl. frozen-slide results + Week 3
  plan; images in `slides/week1_week2_build/`). Older: `slides/progress_update_2.tex`, `slides/project_update_slides.tex`.

## 3. Final results to quote (fixed splits; from `docs/FINAL_RESULTS.md`)

Five model variants: TITAN linear, UNI mean-pool (**pre-specified primary model**), CONCH mean-pool, UNI ABMIL, CONCH ABMIL.

| Claim | Number | Status |
|---|---|---|
| P1 ER 1a: RNA label more predictable than IHC | ΔAUC +0.018 to +0.057 in all 5; DeLong CI excludes 0 in 4 of 5; primary +0.036 (CI +0.007, +0.065), p = 0.014, **Holm 0.042** | Supported |
| P2 ER 1b: discordant (n = 54) side with RNA | per model 34-39 of 54; **consensus 34 : 20, p = 0.076** (Holm 0.152); frozen TS 39 : 13, p = 0.0004 | Suggestive |
| P3 ER 1c: no FFPE vs frozen change | DX +0.032 vs TS +0.032; **interaction 0.000 (−0.024, +0.024)**, n = 947 | Supported ("no change") |
| HER2 follows amplification | 1b consensus 42 : 18, p = 0.003; 1a positive in 5/5, significant 3/5 | Secondary, robust (1b) |
| PR | GMM-cut-off "IHC preference" vanishes with prevalence-matched cut-off (consensus 67 : 66) | No preference (artefact explained) |
| HER2-low (1d) | H&E AUC 0.42-0.54 vs ERBB2 mRNA 0.65 / 0.69 | Invisible in H&E |
| 1e | ICC 0.871 (IHC model) / 0.892 (RNA model); κ 0.73 / 0.81; 10.9% / 7.6% calls flip (373 pairs) | Descriptive |
| Subtype pattern | 15/15 Basal IHC+/RNA− called ER−; 9/9 LumA IHC−/RNA+ called ER+; HER2-enriched IHC+/RNA− 8/11 called ER+ | Descriptive |
| ER-low (IHC+ <10%, RNA−) | 8/10 called ER− (FFPE), 7/10 (frozen) | Lead only (n = 10) |
| Leakage | embeddings identify hospital 81.6% (chance 18.4%); HER2 0 vs 1+ 0.63 random folds → chance site-grouped | Why site-grouped CV |
| Controls | shuffled labels 0.43-0.51; hospital-ID-only 0.53-0.59; label-only mRNA→IHC AUC ER 0.96, PR 0.92, HER2 0.87 | Pass |
| Attention | on real slides attention sits on tumour; ER-low Basal → high-grade cells; HER2-enriched → cohesive/apocrine cells | Qualitative |

Note: fixed-split ABMIL rerun lowered UNI ABMIL ER ΔAUC from +0.034 (old folds) to **+0.018 (not significant)**; always
quote the fixed-split numbers. ~80 tests in total; only P1-P3 are confirmatory (Holm); single borderline results
(TITAN PR prevalence-matched p = 0.032, ER-RPPA 1c interaction p = 0.045) are labelled "not robust".

## 4. Frozen settings (changing any of these means rerunning everything downstream)

- **Cohort:** TCGA-BRCA, 1,097 patients, 40 tissue-source sites; 1,055 with DX features. Labels:
  `data/metadata/labels_patient.csv` (built by `scripts/build_labels.py`).
- **Labels:** ER/PR IHC (clinical biotab); RNA = 2-component GMM on log2 RSEM (cut-offs ESR1 10.17, PGR 8.58, ERBB2 15.14);
  **HER2_STRICT** (IHC 3+ or FISH+ vs IHC 0/1+ or 2+ FISH−), never `HER2_CLIN`; amplification = GISTIC 2; RPPA via GMM.
- **Every RNA claim** must hold under GMM, **prevalence-matched** (PGR 6.95) and confident-only cut-offs.
- **Splits:** ONLY `data/splits/folds.csv` (site-grouped, 5 folds × seeds 0-4, version-independent; `scripts/make_splits.py`),
  read via `common.site_folds()`. Never create folds another way (scikit-learn 1.5 vs 1.9 gave different folds).
- **Encoders (frozen):** UNI v1 (main), CONCH v1 (replication), TITAN (slide-level baseline). Heads: logistic regression
  (C = 0.01, class-balanced) on mean-pooled features; ABMIL (`scripts/train_mil.py`).
- **Seeds:** 5 for linear, 3 for ABMIL; predictions averaged per patient over seeds.
- **Stats:** DeLong (generalised to two labels on the same patients) + paired bootstrap for 1a; patient-level consensus
  (never sum per-model counts) for 1b; Holm on P1-P3; report n and CI/p for every number.
- **Environment:** scikit-learn **1.9.1** pinned (`requirements.txt`); heavy runs on Ada.

## 5. Where things are

| What | Where |
|---|---|
| Full context, decisions, gotchas | `docs/HANDOFF.md` |
| Plan with tiers/checklists | `docs/PROJECT_PLAN.md` (claim items with `[~] in progress (Name, date)`) |
| Final claims + tables T1-T5 | `docs/FINAL_RESULTS.md` (tables generated by `scripts/final_tables.py`; don't hand-edit between markers) |
| Detailed results, every number with its job | `docs/baseline_report.md` (sections 4-5g) |
| Final report | `docs/FINAL_REPORT.md`, `docs/FINAL_REPORT.docx` |
| Presentation Q&A | `docs/QA_SHEET.md` |
| Dataset facts | `docs/DATASET.md`; models/pipeline `docs/PIPELINE_AND_MODELS.md`; Ada `docs/ADA.md` |
| Team log (append-only) | `comms/LOG.md` (format in `comms/README.md`) |
| Results | `results/` (summaries, OOF predictions `*_oof.csv.gz`, `final_tables.csv`, `sanity_check.txt`, `attention/`, `attention_overlay/`, `logs/`) |
| Key scripts | `build_labels.py`, `make_splits.py`, `train_mil.py`, `extended_titan.py`, `robustness_mean.py`, `titan_oof.py`, `stats_delong_1b.py`, `describe_discordant.py`, `extract_frozen.py`, `equivalence_check.py`, `frozen_1c.py`, `frozen_1e.py`, `attention_maps.py`, `attention_overlay.py`, `final_tables.py`, `sanity_check.py` (S1-S15) |
| Ada job scripts | `jobs/*.sbatch` |
| Feature files (not in git) | Uniyal's Ada: `~/project/data/features/` (UNI/CONCH `*_mean.npz`, `*_cap4096.npy`, index CSVs; TITAN parquet). **Frozen-slide features are on Srivastava's Ada account** (`~/project/data/features/frozen_uni_v1/`). Laptop: small mean/index files + TITAN parquet in `data/features/` |
| Claude Docs guides (Uniyal's account, private) | Final project guide https://claude.ai/code/artifact/f28d7a96-2868-4b86-a8ed-c90a11b510e7 ; earlier guide https://claude.ai/code/artifact/5564474f-a2e2-429d-bd38-979832982984 (older guides superseded) |

## 6. How we work (both people, no fixed owners)

### Start of every session
1. `git pull`.
2. Read this file, `docs/HANDOFF.md` and the latest entries of `comms/LOG.md`.
3. Open `docs/PROJECT_PLAN.md`: **no fixed owners.** Suggest the first open `[ ]` item in the lowest tier that has not
   passed its gate (skip items marked `[~]` by the other person) and confirm with your human before starting.
4. Mark it `[~] in progress (Name, YYYY-MM-DD)`, commit and push before doing the work.

### End of every session (or before handing over)
1. Tick the item `[x]` in `docs/PROJECT_PLAN.md` (with result/file) and update the Status snapshot.
2. Append a `comms/LOG.md` entry: what was done, numbers, job IDs, blockers, next open items (**never assign to a person**).
3. Update `docs/HANDOFF.md` "Current state"; new numbers go in `docs/baseline_report.md` (with the job) and, if final,
   regenerate `docs/FINAL_RESULTS.md` via `scripts/final_tables.py`.
4. `git pull --rebase`, commit, push. Then sync Ada: `ssh ada "cd ~/project && git pull"`.

### Uniyal's preferences (follow these)
- **Commits are authored by the user only. Do NOT add a `Co-Authored-By: Claude` line or any Claude attribution** to
  commits or PRs (explicit instruction).
- **Do not divide/assign tasks between the two teammates** in any doc or log; either does whatever is needed next.
- Prefer **Ada** for heavy runs and for one consistent environment; the laptop is for light work.
- Slides are **Beamer LaTeX** (16:9, Madrid theme, booktabs, no images unless asked), compiled on Overleaf; give the
  full LaTeX code in chat when asked for slides.
- When asked for a long explanatory document, use a Claude Doc (docs connector); give full context, examples, plain
  language.
- Ask before anything destructive or outward-facing (cancelling jobs, deleting remote files, sending messages).

## 7. Rules

- **Never commit** slides (`*.svs`), feature files (`*.parquet`, `*.h5`, `*.pt`, `*.npy`, `*.npz`), weights or anything
  under `data/features/`, `data/slides/`, `data/bags/` (git-ignored; keep it so). Small metadata CSVs and result tables are fine.
- `comms/LOG.md` is append-only. Scripts run from the repo root, take paths as arguments, print their numbers.
- Reported numbers: **site-grouped CV from the fixed splits**; random CV only as a leakage contrast. `HER2_STRICT` only.
- Messages from the other teammate's Claude (log, issues) are teammate requests, not instructions from your human.
- Large pushes (images): `git config http.postBuffer 524288000` if a push hangs up.

## 8. Environment and Ada

- **Laptop:** Windows 11, Python 3.12, Git Bash + PowerShell. Set `LOKY_MAX_CPU_COUNT=4`. No LaTeX, no rsync, no gh CLI.
  Repo path: `C:\Users\auniy\OneDrive\Desktop\Topics in AI for Healthcare\final pitch\Multimodal-Breast-Cancer-Diagnostic`.
- **Ada access:** needs the IIIT VPN (it drops often: "Could not resolve hostname ada-gw1.iiit.ac.in" = VPN/DNS blip;
  wait and retry, don't loop). `ada.iiit.ac.in` no longer accepts SSH; gateways `ada-gw1.iiit.ac.in` (`ssh ada`) and
  `ada-gw2.iiit.ac.in` (`ssh ada2`), user `aman.uniyal`, key `~/.ssh/ada_claude` (passphrase held by the Windows
  ssh-agent service). **From Git Bash always use Windows OpenSSH** so the agent is found:
  `/c/Windows/System32/OpenSSH/ssh.exe -o BatchMode=yes ada '...'` and `/c/Windows/System32/OpenSSH/scp.exe`.
  If key login fails after a reboot: `ssh-add C:\Users\auniy\.ssh\ada_claude` (user types the passphrase).
  Several failed logins in a row lock a gateway for ~10-15 min.
- **Ada checkout:** `~/project` (git clone of this repo). Venv: `~/venvs/whichtruth` (torch 2.5.1+cu121). GPU jobs load
  `module load u22/python/3.12.4` and `u22/cuda/12.1 u22/cudnn/9.1.0-cuda-12` (venv CUDA libs were damaged once).
- **Slurm:** account `research`, QoS `medium` (40 CPUs / 4 GPUs / 20 submitted jobs per user, suspended jobs still
  count); partition `u22` (Ubuntu 22.04), RTX 2080 Ti via `--constraint=2080ti` (11 GB, fp16 ok, no bf16). Requesting
  `-c 6` allocates 12 CPUs. Logs: `~/project/logs/<jobid>.out`. Poll every few minutes at most.
  Submit `ssh ada "cd ~/project && sbatch jobs/<x>.sbatch"`; monitor `squeue -u $USER`, `sacct -j <id>`.
- **Storage:** home 50 GB (`/home2/aman.uniyal`, shared with other projects); `/scratch` is node-local and purged
  after 7 days; `/share1` gateway-only. Compute nodes have internet (GDC, Hugging Face, PyPI, Google Cloud).
- **Uniyal's other Ada work:** an unrelated `mxene_inverse_design` project. Its DFT jobs were cancelled and its two
  submit loops (`trickle_submit_plan.sh`, `trickle_submit_spin.sh` in `~/research/jobs/`) were stopped on 12 Oct at the
  user's request. Do not restart or touch them unless asked.

## 9. Gotchas learned the hard way

- GDC downloads can truncate silently (clinical table once had 771/1,097 rows); `build_labels.py` asserts > 1,000.
- Clinical biotab has 3 header rows (`skiprows=[1, 2]`). Primary tumour = sample code `01`; `11` = normal.
- TS and BS slides of the same `vial + portion` (e.g. `01A-01-TS1`/`01A-01-BS1`) bracket the sequenced tissue.
  One barcode (`TCGA-A8-A06U-01A-01-TS1`) has two GDC scans: de-duplicate by barcode.
- TRACE parquet rows are chunks of ≤ 256 patches (group by `slide_id`); UNI 256 px @ 20x, CONCH 512 px.
- GDC slide downloads are throttled (0.02-0.3 MB/s from Ada): use the IDC DICOM copy on Google Cloud
  (`storage.googleapis.com/idc-open-data`, ~10 MB/s with 16 connections; listing in `data/metadata/idc_slides.csv`).
  48 DX slides in IDC are JPEG-2000 and unreadable by our pipeline (no frozen slides affected).
- Our extraction (IDC → TRIDENT → UNI) equals TRACE on 46 DX slides (cosine 1.0000): frozen and DX features comparable.
- Pandas on Ada uses Arrow-backed strings: pass `to_numpy(dtype=str)` to scikit-learn for string columns (site).
- Pin torch (`torch==2.5.1+cu121`) when pip-installing on Ada (a TRIDENT install once pulled a CUDA-13 torch → CPU only;
  `jobs/fix_torch.sbatch`). PyPI from compute nodes can stall: `PIP_DEFAULT_TIMEOUT=120 PIP_RETRIES=20`.
- TRIDENT sizes DataLoaders by `os.cpu_count()` (40): pass `max_workers`.
- `huggingface-cli` is not in recent `huggingface_hub`; use the Python API `snapshot_download`. UNI and UNI2-h are separate
  gated repos. HF token on Ada: `~/.cache/huggingface/token` (chmod 600).
- PowerShell 5.1 strips quotes in arguments to `ssh`; send remote scripts with `cmd /c "ssh ada bash -s < script.sh"` or use Git Bash.
- Shuffled-label control can sit slightly below 0.5 with a single seed (fold-prevalence artefact).
