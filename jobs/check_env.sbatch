#!/bin/bash
# Smoke test: GPU visibility, CUDA modules, Python, internet access from a compute node.
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -p long
#SBATCH -J check_env
#SBATCH -c 10
#SBATCH --gres=gpu:1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:15:00
#SBATCH --output=logs/%j.out

echo "node: $(hostname)  date: $(date)"
echo "CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv
module load u18/cuda/11.6 u18/cudnn/8.4.0-cuda-11.6 2>&1
module list 2>&1
which python3; python3 --version
df -h /scratch /ssd_scratch 2>/dev/null
echo "--- internet from compute node:"
curl -s -o /dev/null -w "GDC %{http_code}\n" --max-time 15 https://api.gdc.cancer.gov/status || echo "GDC unreachable"
curl -s -o /dev/null -w "HuggingFace %{http_code}\n" --max-time 15 https://huggingface.co || echo "HF unreachable"
