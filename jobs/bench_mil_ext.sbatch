#!/bin/bash
# X3 ABMIL regression, X4b multi-task (+ single-task reference), X4c DANN, X4d noisy labels, X6 MC-dropout / half-bags
# (scripts/bench_mil_ext.py) on UNI bags. One mode per array task:
#   sbatch --array=0-5%2 jobs/bench_mil_ext.sbatch
#   sbatch --array=0 --export=ALL,MODES="noisy",EXTRA_ARGS=--smoke jobs/bench_mil_ext.sbatch
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -p u22
#SBATCH --constraint=2080ti
#SBATCH -J mil_ext
#SBATCH -c 6
#SBATCH --gres=gpu:1
#SBATCH --mem-per-cpu=3G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%A_%a.out
set -e
MODES=(${MODES:-regress multitask single dann noisy reliab})
M=${MODES[${SLURM_ARRAY_TASK_ID:-0}]}
case $M in regress) DEF_OUT=results/bench/x3 ;; reliab) DEF_OUT=results/bench/x6 ;; *) DEF_OUT=results/bench/x4 ;; esac
module load u22/python/3.12.4
module load u22/cuda/12.1 u22/cudnn/9.1.0-cuda-12
source $HOME/venvs/whichtruth/bin/activate
cd $HOME/project
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
python -u scripts/bench_mil_ext.py --mode $M --name ${NAME:-uni_v1_20x} --out ${OUT:-$DEF_OUT} ${EXTRA_ARGS}
echo MIL_EXT_JOB_OK
