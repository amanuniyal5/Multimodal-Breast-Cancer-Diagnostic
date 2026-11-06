#!/bin/bash
# X1: stream one TRACE encoder folder (47-121 GB of parquet) through scripts/trace_means.py on node-local /scratch and
# keep only patient means (+SD) in ~/project/data/features (a few MB). CPU only, download-bound (~10 MB/s per job).
#   sbatch --export=ALL,P=uni_v2/20x,NAME=uni_v2_20x jobs/trace_means.sbatch
#   sbatch --export=ALL,P=chief/10x,NAME=chief_slide,EXTRA_ARGS=--slide-level jobs/trace_means.sbatch
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -p u22
#SBATCH -J trace_means
#SBATCH -c 2
#SBATCH --mem-per-cpu=3G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.out
set -e
module load u22/python/3.12.4
source $HOME/venvs/whichtruth/bin/activate
cd $HOME/project
W=/scratch/$USER/tm_${SLURM_JOB_ID}
echo "node $(hostname); $P -> $NAME; scratch free $(df -h /scratch | tail -1 | awk '{print $4}')"
export HF_HUB_DISABLE_PROGRESS_BARS=1
python -u scripts/trace_means.py --path data/TCGA-BRCA/$P --work $W --out data/features --name $NAME ${EXTRA_ARGS}
rm -rf $W
