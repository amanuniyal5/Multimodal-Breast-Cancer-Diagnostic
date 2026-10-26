#!/bin/bash
# 1c and 1e with the ABMIL head (scripts/frozen_abmil.py). Array of 5 training tasks, then the CPU analysis:
#   sbatch --array=0-4%3 jobs/frozen_abmil.sbatch
#   sbatch --dependency=afterok:<array id> jobs/week3_stats.sbatch
# Tasks: 0 DX er, 1 TS er, 2 DX prher2, 3 TS prher2, 4 TS 1e. Needs data/features/uni_v1_20x_* and frozen_uni_v1/.
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -p u22
#SBATCH --constraint=2080ti
#SBATCH -J frozen_abmil
#SBATCH -c 6
#SBATCH --gres=gpu:1
#SBATCH --mem-per-cpu=3G
#SBATCH --time=10:00:00
#SBATCH --output=logs/%x_%A_%a.out
set -e
module load u22/python/3.12.4
module load u22/cuda/12.1 u22/cudnn/9.1.0-cuda-12  # venv CUDA libs were damaged (docs/HANDOFF.md)
source $HOME/venvs/whichtruth/bin/activate
cd $HOME/project
export LOKY_MAX_CPU_COUNT=4
SLIDES=(DX TS DX TS TS)
PARTS=(er er prher2 prher2 1e)
i=${SLURM_ARRAY_TASK_ID:-0}
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
python -u scripts/frozen_abmil.py --slide ${SLIDES[$i]} --part ${PARTS[$i]} --seeds 3
echo FROZEN_ABMIL_OK
