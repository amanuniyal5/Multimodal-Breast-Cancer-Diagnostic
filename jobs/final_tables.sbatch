#!/bin/bash
# Tier 5: final tables from fixed-split runs only (DeLong 1a + combined 1b with fixed-split ABMIL, then
# scripts/final_tables.py into docs/FINAL_RESULTS.md) and the full sanity check (S1-S15). CPU only.
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -p u22
#SBATCH -J final_tables
#SBATCH -c 4
#SBATCH --mem-per-cpu=3G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%j.out
set -e
module load u22/python/3.12.4
source $HOME/venvs/whichtruth/bin/activate
cd $HOME/project
export LOKY_MAX_CPU_COUNT=4
python -u scripts/stats_delong_1b.py --abmil fixed --tag _final --boot 2000
python -u scripts/final_tables.py
python -u scripts/sanity_check.py
echo FINAL_OK
