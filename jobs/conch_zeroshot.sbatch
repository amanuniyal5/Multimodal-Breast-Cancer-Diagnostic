#!/bin/bash
# Week 3 label-free audit with zero-shot CONCH (scripts/conch_zeroshot.py):
#   1. install the CONCH package if missing (torch/torchvision/numpy/timm pinned to the venv's versions);
#   2. day-1 check: re-extract N DX slides with TRIDENT conch_v1 (512 px at 20x, as TRACE) and find out whether the
#      TRACE CONCH features are before or after the contrastive projection (results/conch_feature_check.json);
#   3. zero-shot scores, 1a/1b without training labels, concept scores and the concept-only model.
# Needs data/features/conch_v1_20x_* and the Hugging Face token with CONCH access (~/.cache/huggingface/token).
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -p u22
#SBATCH --constraint=2080ti
#SBATCH -J conch_zeroshot
#SBATCH -c 6
#SBATCH --gres=gpu:1
#SBATCH --mem-per-cpu=3G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%j.out
set -e
module load u22/python/3.12.4
module load u22/cuda/12.1 u22/cudnn/9.1.0-cuda-12  # venv CUDA libs were damaged (docs/HANDOFF.md)
source $HOME/venvs/whichtruth/bin/activate
cd $HOME/project
export LOKY_MAX_CPU_COUNT=4 PIP_DEFAULT_TIMEOUT=120 PIP_RETRIES=20
N=${N:-5}
if ! python -c "import conch.open_clip_custom" 2>/dev/null; then
    C=$HOME/pip_constraints_$SLURM_JOB_ID.txt
    pip freeze | grep -iE '^(torch|torchvision|numpy|timm)==' > $C
    cat $C
    pip install -c $C "https://github.com/Mahmoodlab/CONCH/archive/refs/heads/main.zip"
    rm -f $C
fi
python -c "import torch, conch.open_clip_custom; print('torch', torch.__version__, 'cuda', torch.cuda.is_available()); \
assert torch.__version__.startswith('2.5.1'), 'torch was changed by the CONCH install'"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

if [ ! -f results/conch_feature_check.json ]; then
    SLIDES=$(python scripts/equivalence_check.py --select $N)
    echo "CONCH space check on $N DX slides: $SLIDES"
    W=/scratch/$USER/conch_check_${SLURM_JOB_ID}
    python -u scripts/extract_frozen.py --slides "$SLIDES" --work $W --out data/features/dx_check_conch_v1 \
        --encoder conch_v1 --patch-size 512
    rm -rf $W
    python -u scripts/equivalence_check.py --compare --name conch_v1_20x --ours data/features/dx_check_conch_v1 \
        --out results/equivalence_check_conch.csv
    python -u scripts/conch_zeroshot.py --check --ours data/features/dx_check_conch_v1
fi
python -u scripts/conch_zeroshot.py
echo CONCH_ZEROSHOT_OK
