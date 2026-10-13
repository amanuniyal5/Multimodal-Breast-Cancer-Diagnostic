#!/bin/bash
# Tier 4: download + encode one manifest batch (~50 slides, ~15 GB) of frozen TS/BS slides with TRIDENT + UNI v1.
# Raw slides live on node-local /scratch only while being encoded; per-slide features go to data/features/.
# Timing test (random-weight ViT-L, 3 slides):
#   sbatch --export=ALL,ENCODER=vitl_random,LIMIT=3,BATCH=0 jobs/extract_frozen.sbatch
# All 30 batches (needs UNI access) as 6 array tasks of 5 batches each (the QoS allows 20 submitted jobs per user):
#   sbatch --array=0-5%2 jobs/extract_frozen.sbatch
# Slides already in the output folder are skipped, so any task can simply be resubmitted.
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -p u22
#SBATCH --constraint=2080ti
#SBATCH -J extract_frozen
#SBATCH -c 6
#SBATCH --gres=gpu:1
#SBATCH --mem-per-cpu=3G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/%x_%A_%a.out
set -e
module load u22/python/3.12.4
source $HOME/venvs/whichtruth/bin/activate
cd $HOME/project
ENCODER=${ENCODER:-uni_v1}
PER_TASK=${PER_TASK:-5}
if [ -z "$BATCH" ] && [ -n "$SLURM_ARRAY_TASK_ID" ]; then
  BATCH=$(( SLURM_ARRAY_TASK_ID * PER_TASK ))-$(( SLURM_ARRAY_TASK_ID * PER_TASK + PER_TASK - 1 ))
fi
BATCH=${BATCH:-0}
LIMIT=${LIMIT:-0}
MANIFEST=${MANIFEST:-data/metadata/frozen_manifest.csv}
OUT=${OUT:-data/features/frozen_${ENCODER}}
W=/scratch/$USER/frozen_${SLURM_JOB_ID}
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
echo "batch $BATCH encoder $ENCODER limit $LIMIT -> $OUT; scratch free $(df -h /scratch | tail -1 | awk '{print $4}')"
python -u scripts/extract_frozen.py --manifest $MANIFEST --batch $BATCH --limit $LIMIT --work $W --out $OUT \
    --encoder $ENCODER ${EXTRA_ARGS}
rm -rf $W
echo EXTRACT_OK
