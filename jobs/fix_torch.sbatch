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
# torch 2.14 also brought CUDA-13 nvidia-* packages; they share the nvidia/ folders with the CUDA-12 ones and replaced
# e.g. libcudnn.so.9 (cuDNN 9.24 for CUDA 13 -> CUDNN_STATUS_NOT_INITIALIZED, job 20506). Remove them, then reinstall
# torch 2.5.1 with all of its own CUDA-12 dependencies (force) from the PyTorch index.
python -m pip uninstall -y -q nvidia-cublas nvidia-cuda-cupti nvidia-cuda-nvrtc nvidia-cuda-runtime nvidia-cudnn-cu13     nvidia-cufft nvidia-cufile nvidia-curand nvidia-cusolver nvidia-cusparse nvidia-cusparselt-cu13 nvidia-nccl-cu13     nvidia-nvjitlink nvidia-nvshmem-cu13 nvidia-nvtx || true
python -m pip install --no-cache-dir -q --force-reinstall torch==2.5.1+cu121 torchvision==0.20.1+cu121     --index-url https://download.pytorch.org/whl/cu121 --extra-index-url https://pypi.org/simple
python -m pip list 2>/dev/null | grep -iE "^(torch|torchvision|triton|nvidia-cudnn|numpy|pillow) "
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
