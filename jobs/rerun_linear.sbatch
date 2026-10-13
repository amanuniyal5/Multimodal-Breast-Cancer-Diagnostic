#!/bin/bash
# Rerun every linear (logistic-regression) baseline on the fixed splits in data/splits/folds.csv, in one environment:
# TITAN (B2-B6, E1-E5), robustness (PR cut-off, ER strata), UNI + CONCH mean-pool (all experiments). CPU only.
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -p u22
#SBATCH -J rerun_linear
#SBATCH -c 8
#SBATCH --mem-per-cpu=3G
#SBATCH --time=03:00:00
#SBATCH --output=logs/%j.out
set -e
module load u22/python/3.12.4
source $HOME/venvs/whichtruth/bin/activate
cd $HOME/project
pip install -q -r requirements.txt
python -c "import sklearn; print('scikit-learn', sklearn.__version__)"
export LOKY_MAX_CPU_COUNT=8
[ -f data/features/titan_brca_20x.parquet ] || python scripts/fetch_titan.py
echo "=== baseline_titan (B2-B4, 1b)"; python -u scripts/baseline_titan.py
echo "=== paired_ihc_vs_rna (B5, 1b robustness)"; python -u scripts/paired_ihc_vs_rna.py | tee results/paired_ihc_vs_rna.txt
echo "=== extended_titan (E1-E5, 5 seeds)"; python -u scripts/extended_titan.py --seeds 5
echo "=== robustness_mean (R1, R2, 5 seeds)"; python -u scripts/robustness_mean.py --seeds 5
for N in uni_v1_20x conch_v1_20x; do
  echo "=== mean-pool $N (5 seeds)"
  python -u scripts/train_mil.py --name $N --models mean --seeds 5 --tag _fixedsplits
done
echo RERUN_OK
