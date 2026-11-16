#!/bin/bash
# X5b: stream one TRACE encoder folder and keep ALL DX patch features + coords of the map patients
# (data/metadata/maps_patients.csv) in data/features/maps_<NAME>/ (~0.3-1 GB). Means go to check_<NAME>_mean.npz
# (never overwrites the stored means). CPU only.
#   sbatch --export=ALL,P=uni_v1/20x,NAME=uni_v1_20x jobs/trace_keep.sbatch
#   sbatch --export=ALL,P=conch_v1/20x,NAME=conch_v1_20x jobs/trace_keep.sbatch
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -p u22
#SBATCH -J trace_keep
#SBATCH -c 2
#SBATCH --mem-per-cpu=3G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.out
set -e
module load u22/python/3.12.4
source $HOME/venvs/whichtruth/bin/activate
cd $HOME/project
W=/scratch/$USER/tk_${SLURM_JOB_ID}
trap "rm -rf $W" EXIT
echo "node $(hostname); $P -> maps_$NAME; home used $(du -sh $HOME 2>/dev/null | cut -f1)"
export HF_HUB_DISABLE_PROGRESS_BARS=1
python -u scripts/trace_means.py --path data/TCGA-BRCA/$P --work $W --out data/features --name check_$NAME \
  --keep data/metadata/maps_patients.csv --keep-dir data/features/maps_$NAME
echo TRACE_KEEP_JOB_OK
