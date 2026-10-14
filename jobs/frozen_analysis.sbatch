#!/bin/bash
# Experiments 1c (DX vs frozen TS x IHC vs RNA) and 1e (TS vs BS agreement) on the extracted frozen-slide features.
# CPU only, minutes. Run after jobs/extract_frozen.sbatch has finished (it also works on a partial set, for testing).
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -p u22
#SBATCH -J frozen_analysis
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
echo "frozen slides extracted: $(ls data/features/frozen_uni_v1/*.npz | wc -l)"
echo "=== 1c"; python -u scripts/frozen_1c.py ${ARGS_1C}
echo "=== 1e"; python -u scripts/frozen_1e.py ${ARGS_1E}
echo ANALYSIS_OK
