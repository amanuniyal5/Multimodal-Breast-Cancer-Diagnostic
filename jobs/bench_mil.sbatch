#!/bin/bash
# X2 aggregation methods (scripts/bench_mil.py). One method per array task, e.g.
#   sbatch --array=0-4%3 jobs/bench_mil.sbatch                       # abmil clam transmil dsmil acmil on UNI
#   sbatch --array=0 --export=ALL,METHODS="clam",EXTRA_ARGS=--smoke jobs/bench_mil.sbatch
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -p u22
#SBATCH --constraint=2080ti
#SBATCH -J bench_mil
#SBATCH -c 6
#SBATCH --gres=gpu:1
#SBATCH --mem-per-cpu=3G
#SBATCH --time=16:00:00
#SBATCH --output=logs/%x_%A_%a.out
set -e
METHODS=(${METHODS:-abmil clam transmil dsmil acmil})
M=${METHODS[${SLURM_ARRAY_TASK_ID:-0}]}
NAME=${NAME:-uni_v1_20x}
module load u22/python/3.12.4
module load u22/cuda/12.1 u22/cudnn/9.1.0-cuda-12  # venv CUDA libs were damaged (docs/HANDOFF.md)
source $HOME/venvs/whichtruth/bin/activate
cd $HOME/project
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
python -u scripts/bench_mil.py --method $M --name $NAME --out ${OUT:-results/bench/x2} ${EXTRA_ARGS}
echo BENCH_MIL_JOB_OK
