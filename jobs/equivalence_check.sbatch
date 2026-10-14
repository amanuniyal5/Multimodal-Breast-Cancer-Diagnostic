#!/bin/bash
# Tier 4 extraction-equivalence check: 50 DX slides from IDC through our TRIDENT + UNI pipeline, compared with TRACE.
# Needs data/features/uni_v1_20x_* (jobs/uni_features.sbatch) and data/metadata/idc_slides.csv.
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -p u22
#SBATCH --constraint=2080ti
#SBATCH -J equiv_check
#SBATCH -c 6
#SBATCH --gres=gpu:1
#SBATCH --mem-per-cpu=3G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%j.out
set -e
module load u22/python/3.12.4
module load u22/cuda/12.1 u22/cudnn/9.1.0-cuda-12
source $HOME/venvs/whichtruth/bin/activate
cd $HOME/project
N=${N:-50}
SLIDES=$(python scripts/equivalence_check.py --select $N)
echo "selected $N DX slides: $SLIDES"
W=/scratch/$USER/equiv_${SLURM_JOB_ID}
python -u scripts/extract_frozen.py --slides "$SLIDES" --work $W --out data/features/dx_check_uni_v1 --encoder uni_v1
rm -rf $W
python -u scripts/equivalence_check.py --compare --ours data/features/dx_check_uni_v1
echo EQUIV_OK
