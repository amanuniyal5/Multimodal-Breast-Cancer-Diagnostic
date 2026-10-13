#!/bin/bash
# Tier 2 statistics: TITAN out-of-fold predictions on the fixed splits (GMM and prevalence-matched PR cut-off),
# then DeLong tests (1a) and combined evidence across the five model variants (1b). CPU only, ~10 min.
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -p u22
#SBATCH -J stats_tier2
#SBATCH -c 4
#SBATCH --mem-per-cpu=3G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%j.out
set -e
module load u22/python/3.12.4
source $HOME/venvs/whichtruth/bin/activate
cd $HOME/project
python -c "import sklearn; print('scikit-learn', sklearn.__version__)"
export LOKY_MAX_CPU_COUNT=4
[ -f data/features/titan_brca_20x.parquet ] || python scripts/fetch_titan.py
echo "=== titan_oof (GMM)"; python -u scripts/titan_oof.py --seeds 5
echo "=== titan_oof (PR prevalence-matched)"; python -u scripts/titan_oof.py --seeds 5 --pr-cut prevmatch --tag _prprev
echo "=== stats_delong_1b"; python -u scripts/stats_delong_1b.py --selftest
python -u scripts/stats_delong_1b.py --boot 2000
echo STATS_OK
