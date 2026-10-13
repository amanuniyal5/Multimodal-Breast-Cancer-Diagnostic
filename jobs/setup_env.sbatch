#!/bin/bash
# One-time: Python venv in home with PyTorch (CUDA 12.1 wheels) + project requirements.
# A fresh install on the NFS home took more than 1 h (job 19985 timed out), hence 3 h.
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -p u22
#SBATCH --constraint=2080ti
#SBATCH -J setup_env
#SBATCH -c 4
#SBATCH --gres=gpu:1
#SBATCH --mem-per-cpu=3G
#SBATCH --time=03:00:00
#SBATCH --output=logs/%j.out
set -e
module load u22/python/3.12.4
VENV=$HOME/venvs/whichtruth
[ -d "$VENV" ] || python3 -m venv "$VENV"
source "$VENV/bin/activate"
# PyPI reads from compute nodes sometimes time out (job 20011): longer timeout, more retries; rerunnable
export PIP_DEFAULT_TIMEOUT=120 PIP_RETRIES=10
python -c "import torch" 2>/dev/null || python -m pip install --no-cache-dir -q torch==2.5.1 --index-url https://download.pytorch.org/whl/cu121
python -m pip install --no-cache-dir -q -r requirements.txt huggingface_hub
python -c "import torch; print('torch', torch.__version__, 'cuda', torch.cuda.is_available(), torch.cuda.get_device_name(0))"
du -sh "$VENV"
echo SETUP_OK
