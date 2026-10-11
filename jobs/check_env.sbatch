#!/bin/bash
# Smoke test: GPU visibility, OS, modules, Python, internet access and disk from a compute node.
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -p u22
#SBATCH --constraint=2080ti
#SBATCH -J check_env
#SBATCH -c 10
#SBATCH --gres=gpu:1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:15:00
#SBATCH --output=logs/%j.out

echo "node: $(hostname)  date: $(date)"
grep PRETTY_NAME /etc/os-release
echo "CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv
echo "--- modules (cuda/python/conda):"
module avail 2>&1 | tr -s ' ' '\n' | grep -i -E "cuda|cudnn|python|conda|anaconda|miniforge" | sort -u | head -40
echo "--- python:"
which python3; python3 --version
echo "--- disk:"
df -h /scratch /ssd_scratch /share1 $HOME 2>/dev/null
echo "--- internet from compute node:"
curl -s -o /dev/null -w "GDC %{http_code}\n" --max-time 15 https://api.gdc.cancer.gov/status || echo "GDC unreachable"
curl -s -o /dev/null -w "HuggingFace %{http_code}\n" --max-time 15 https://huggingface.co || echo "HF unreachable"
curl -s -o /dev/null -w "PyPI %{http_code}\n" --max-time 15 https://pypi.org/simple/ || echo "PyPI unreachable"
echo DONE
