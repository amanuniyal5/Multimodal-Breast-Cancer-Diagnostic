# Foundation models and pipeline

## Model choice

| Model | Role | Architecture / dim | TCGA-BRCA features (TRACE) | Weights | Pretrained on TCGA? |
|---|---|---|---|---|---|
| **UNI v1** | Main patch encoder | ViT-L/16, 1024-d, 256 px | 20×: 62 GB · 10×: 16 GB · 5× | `MahmoodLab/UNI`, gated (IIIT e-mail), CC-BY-NC-ND | No (Mass-100K) |
| **CONCH v1** | Robustness encoder | ViT-B/16, 512-d, 512 px | 20×: 8 GB | `MahmoodLab/CONCH`, gated, CC-BY-NC-ND | No |
| **TITAN** | Laptop baselines | slide encoder on CONCH v1.5, 768-d per slide | 6.8 MB (1,125 DX slides) | gated, CC-BY-NC-ND | No (verify in paper) |
| UNI2-h / Virchow2 / GigaPath | Not used | ViT-H / ViT-g | 93 GB each | gated / open | No; 2–4× slower to extract |
| Phikon / CTransPath / Kaiko | Avoid | ViT-B / Swin-T | partly | open | **Yes** (leakage) |

The encoder is **frozen**. We train only heads: ABMIL (gated attention, ~0.5–1M parameters) and a mean-pool
logistic-regression baseline.

## Hardware notes (Ada)

- Nodes: 4× GTX 1080 Ti or 4× RTX 2080 Ti (11 GB each), 40 cores, 128 GB RAM, 1.8 TB local scratch, SLURM, Ubuntu 18.04.
- 1080 Ti: fp32 only (fp16 is very slow on Pascal). 2080 Ti: `torch.autocast("cuda", dtype=torch.float16)`. No bf16.
- UNI ViT-L inference fits easily at batch 64–128.

## Step by step

1. **Access.** Request `MahmoodLab/UNI` and `MahmoodLab/CONCH` on Hugging Face with an IIIT e-mail as the
   account's primary address. Then `huggingface-cli login` on Ada.
2. **DX features (no extraction).**
   ```bash
   huggingface-cli download IUCompPath/TRACE --repo-type dataset \
       --include "data/TCGA-BRCA/uni_v1/20x/*" --local-dir /scratch/$USER/trace
   ```
   63 parquet files of ~1 GB. Columns: `slide_id, case_submitter_id, magnification, patch_size_px, chunk_id,
   patch_start_idx, num_patches, embedding_dim, coords, embeddings, svs_filename, gdc_file_id`.
3. **Regroup into per-slide bags.** Each row is a chunk of up to 256 patches. Stream one file at a time with
   `pyarrow`, concatenate chunks per `slide_id` (ordered by `chunk_id`), save `bags/<slide_id>.pt` holding
   `{"feats": N×1024, "coords": N×2}` (fp16 halves the 62 GB).
4. **Frozen slides (1c, 1e).**
   - Build a GDC manifest of tumour TS/BS slides from `data/metadata/gdc_slides.csv`
     (sample code `01`, `TS`/`BS` in the name): 1,517 slides, 430 GB.
   - Download in batches of ~300 with `gdc-client download -m manifest_batchK.txt`.
   - Run [TRIDENT](https://github.com/mahmoodlab/TRIDENT): tissue segmentation → 256 px patches at 20× → UNI,
     fp16 on 2080 Ti. Save bags in the same format as step 3. Delete raw slides after each batch.
   - Run the same TRIDENT settings on ~50 DX slides and compare with TRACE (cosine similarity per patch,
     agreement of predictions) so FFPE vs frozen differs only in tissue, not preprocessing.
5. **Heads.** ABMIL and mean-pool LR, one model per label source; patient-level, site-grouped 5-fold CV × 5 seeds.
   Fixed split files shared in `data/splits/` (to be created once, used by both of us).
6. **CONCH replication** with TRACE `conch_v1/20x` for the headline experiments.

## Quick laptop path (already working)

`scripts/fetch_titan.py` → `scripts/baseline_titan.py` / `scripts/paired_ihc_vs_rna.py`: TITAN slide embeddings +
logistic regression. Used for the preliminary results only.
