#!/bin/bash
# Tier 6c: CPU benchmarks on saved patient means (X1 encoders + label efficiency, X1 ResNet-50 subset, X3 regression,
# X4a/f fusion + calibration, X5a/c/e embeddings; x2pool = X2 pooling baselines, x2an = X2 tables). Pick steps with STEPS (default all), e.g.
#   sbatch --export=ALL,STEPS="x1 x3" jobs/bench_cpu.sbatch
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -p u22
#SBATCH -J bench_cpu
#SBATCH -c 4
#SBATCH --mem-per-cpu=3G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.out
set -e
module load u22/python/3.12.4
source $HOME/venvs/whichtruth/bin/activate
cd $HOME/project
export LOKY_MAX_CPU_COUNT=8 SLURM_CPUS_PER_TASK=8
STEPS=${STEPS:-"x1 x1sub x3 x4 x5"}
for s in $STEPS; do
  echo "===== $s ====="
  case $s in
    x1)    python -u scripts/bench_encoders.py --std --out results/bench/x1 ${X1_ARGS} ;;
    x1sub) python -u scripts/bench_encoders.py --cohort subset --out results/bench/x1_subset ${X1SUB_ARGS} ;;
    x3)    python -u scripts/bench_regression.py --x1 results/bench/x1 --out results/bench/x3 ;;
    x4)    python -u scripts/bench_fusion.py --x1 results/bench/x1 --out results/bench/x4 ;;
    x5)    python -u scripts/embeddings.py --out results/bench/x5 ;;
    x2pool) python -u scripts/bench_mil.py --method pool --name ${NAME:-uni_v1_20x} --out results/bench/x2 ;;
    x2an)  python -u scripts/analyse_x2.py --dir results/bench/x2 --name ${NAME:-uni_v1_20x} ;;
  esac
done
echo BENCH_CPU_OK
