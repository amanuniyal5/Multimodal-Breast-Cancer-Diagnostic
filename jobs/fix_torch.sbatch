#!/bin/bash
# Restore the pinned PyTorch build after a dependency upgraded it. TRIDENT's `pip install -e .` pulled torch 2.14.1+cu130
# (job 20078), which the Ada driver (570, CUDA <= 12.8) cannot use, so everything silently ran on CPU (job 20079).
# Fails unless the GPU is usable, so dependent jobs (afterok) only start with a working GPU stack.
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -p u22
#SBATCH --constraint=2080ti
#SBATCH -J fix_torch
#SBATCH -c 4
#SBATCH --gres=gpu:1
#SBATCH --mem-per-cpu=3G
#SBATCH --time=03:00:00
#SBATCH --output=logs/%j.out
set -e
module load u22/python/3.12.4
source $HOME/venvs/whichtruth/bin/activate
cd $HOME/project
export PIP_DEFAULT_TIMEOUT=120 PIP_RETRIES=20
python -m pip install --no-cache-dir -q --no-deps --force-reinstall torch==2.5.1+cu121 torchvision==0.20.1+cu121 \
    --index-url https://download.pytorch.org/whl/cu121
python -m pip install --no-cache-dir -q --no-deps triton==3.1.0
python - <<'PY'
import torch, torchvision, timm, trident
print("torch", torch.__version__, "| torchvision", torchvision.__version__, "| timm", timm.__version__)
assert torch.cuda.is_available(), "CUDA not available"
x = torch.randn(1024, 1024, device="cuda", dtype=torch.float16)
print("GPU", torch.cuda.get_device_name(0), "matmul ok", float((x @ x).float().abs().mean()) > 0)
from trident.patch_encoder_models import encoder_factory
print("UNI forward:", encoder_factory("uni_v1").cuda().eval()(torch.randn(2, 3, 224, 224, device="cuda")).shape)
PY
echo FIX_TORCH_OK
