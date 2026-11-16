#!/bin/bash
# X5b patch-level maps (scripts/maps.py): needs data/features/maps_uni_v1_20x/ and maps_conch_v1_20x/ from
# jobs/trace_keep.sbatch, the UNI / CONCH cap4096 bags and the CONCH package (concept maps). Slide thumbnails from IDC.
#   sbatch --dependency=afterok:<keep_uni>:<keep_conch> jobs/maps.sbatch
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -p u22
#SBATCH --constraint=2080ti
#SBATCH -J maps
#SBATCH -c 6
#SBATCH --gres=gpu:1
#SBATCH --mem-per-cpu=3G
#SBATCH --time=8:00:00
#SBATCH --output=logs/%x_%j.out
set -e
module load u22/python/3.12.4
module load u22/cuda/12.1 u22/cudnn/9.1.0-cuda-12
source $HOME/venvs/whichtruth/bin/activate
cd $HOME/project
W=/scratch/$USER/maps_${SLURM_JOB_ID}
mkdir -p $W
trap "rm -rf $W" EXIT
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
python -u scripts/maps.py --work $W ${EXTRA_ARGS}
echo MAPS_JOB_OK
