#!/bin/bash
# Option-1 experiments (1a, 1b, 1d, shuffled control) with mean-pool LR and ABMIL on patch features.
# Override with: sbatch --export=ALL,NAME=conch_v1_20x,MODELS=mean jobs/train_mil.sbatch [extra train_mil.py args]
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -p u22
#SBATCH --constraint=2080ti
#SBATCH -J train_mil
#SBATCH -c 6
#SBATCH --gres=gpu:1
#SBATCH --mem-per-cpu=3G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%j.out
set -e
NAME=${NAME:-uni_v1_20x}
MODELS=${MODELS:-mean,abmil}
SEEDS=${SEEDS:-3}
module load u22/python/3.12.4
module load u22/cuda/12.1 u22/cudnn/9.1.0-cuda-12  # venv CUDA libs were damaged (docs/HANDOFF.md)
source $HOME/venvs/whichtruth/bin/activate
cd $HOME/project
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
python -u scripts/train_mil.py --feats data/features --name $NAME --cap 4096 \
    --labels data/metadata/labels_patient.csv --out results --models $MODELS --seeds $SEEDS "$@"
echo TRAIN_OK
