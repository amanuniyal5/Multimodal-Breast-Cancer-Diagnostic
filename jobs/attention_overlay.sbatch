#!/bin/bash
# Tier 3: attention heatmaps overlaid on the real DX slides (scripts/attention_overlay.py). CPU only; downloads the
# example slides from the IDC copy to node-local /scratch and deletes them afterwards.
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -p u22
#SBATCH -J attn_overlay
#SBATCH -c 4
#SBATCH --mem-per-cpu=3G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%j.out
set -e
module load u22/python/3.12.4
source $HOME/venvs/whichtruth/bin/activate
cd $HOME/project
P=${PATIENTS:-TCGA-A2-A4S1,TCGA-AR-A0TP,TCGA-AR-A0TX,TCGA-E2-A152,TCGA-A1-A0SH,TCGA-B6-A0RG,TCGA-AN-A0FS,TCGA-E2-A1LL,TCGA-LL-A5YP,TCGA-BH-A0WA}
W=/scratch/$USER/overlay_${SLURM_JOB_ID}
python -u scripts/attention_overlay.py --att-dir ${ATT:-results/attention_v2} --patients $P --work $W --out results/attention_overlay
rm -rf $W
echo OVERLAY_OK
