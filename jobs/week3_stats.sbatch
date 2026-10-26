#!/bin/bash
# Week 3 analysis on saved outputs (CPU, minutes): frozen-slide ABMIL summaries (scripts/frozen_abmil.py --analyse,
# after jobs/frozen_abmil.sbatch) and the extended statistics (scripts/stats_extended.py).
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -p u22
#SBATCH -J week3_stats
#SBATCH -c 4
#SBATCH --mem-per-cpu=3G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%j.out
set -e
module load u22/python/3.12.4
source $HOME/venvs/whichtruth/bin/activate
cd $HOME/project
export LOKY_MAX_CPU_COUNT=4
python -c "import sklearn; print('scikit-learn', sklearn.__version__)"
if ls results/frozen_abmil_raw_*.csv.gz >/dev/null 2>&1; then python -u scripts/frozen_abmil.py --analyse; fi
python -u scripts/stats_extended.py
echo WEEK3_STATS_OK
