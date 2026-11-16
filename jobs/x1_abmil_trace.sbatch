#!/bin/bash
# X1 ABMIL on a second encoder without using home quota: stream the TRACE folder once into node-local /scratch
# (patient means + a 2,048-patch reservoir bag per patient + coords), then run bench_mil.py ABMIL on those bags in the
# same job, and delete /scratch at the end. Default: H-optimus-0 (joint best mean-pool ER IHC AUC with UNI, X1).
#   sbatch jobs/x1_abmil_trace.sbatch
#   sbatch --export=ALL,P=uni_v2/20x,NAME=uni_v2_20x jobs/x1_abmil_trace.sbatch
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -p u22
#SBATCH --constraint=2080ti
#SBATCH -J x1_abmil
#SBATCH -c 6
#SBATCH --gres=gpu:1
#SBATCH --mem-per-cpu=3G
#SBATCH --time=20:00:00
#SBATCH --output=logs/%x_%j.out
set -e
P=${P:-hoptimus0/20x}
NAME=${NAME:-hoptimus0_20x}
K=${K:-2048}
module load u22/python/3.12.4
module load u22/cuda/12.1 u22/cudnn/9.1.0-cuda-12
source $HOME/venvs/whichtruth/bin/activate
cd $HOME/project
W=/scratch/$USER/x1a_${SLURM_JOB_ID}
mkdir -p $W/feats
trap "rm -rf $W" EXIT
echo "node $(hostname); $P -> $NAME cap $K; scratch free $(df -h /scratch | tail -1 | awk '{print $4}')"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
python -c "import torch; assert torch.cuda.is_available(), 'no usable GPU'"
export HF_HUB_DISABLE_PROGRESS_BARS=1
python -u scripts/trace_means.py --path data/TCGA-BRCA/$P --work $W/dl --out $W/feats --name $NAME --bags $K
python - <<EOF
import numpy as np
a = np.load("$W/feats/${NAME}_mean.npz"); b = np.load("data/features/${NAME}_mean.npz")
ia = {p: i for i, p in enumerate(a["patients"])}
d = max(abs(a["X"][ia[p]] - b["X"][i]).max() for i, p in enumerate(b["patients"]) if p in ia)
print(f"means vs stored ${NAME}_mean.npz: {len(set(a['patients']) & set(b['patients']))} common patients, max abs diff {d:.2e}")
EOF
python -u scripts/bench_mil.py --method abmil --name $NAME --feats $W/feats --cap $K --out results/bench/x1_abmil
echo X1_ABMIL_JOB_OK
