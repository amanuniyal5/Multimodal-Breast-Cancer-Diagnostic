#!/bin/bash
# Download TRACE UNI v1 20x patch features for TCGA-BRCA (62 GB) to node-local /scratch, build patient-level
# features with scripts/build_bags.py, keep only the outputs (~9 GB) in ~/project/data/features.
# Override the encoder with: sbatch --export=ALL,MODEL=conch_v1 jobs/uni_features.sbatch
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -p u22
#SBATCH -J uni_features
#SBATCH -c 4
#SBATCH --mem-per-cpu=3G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%j.out
set -e
MODEL=${MODEL:-uni_v1}
MAG=${MAG:-20x}
CAP=${CAP:-4096}
module load u22/python/3.12.4
source $HOME/venvs/whichtruth/bin/activate
cd $HOME/project
W=/scratch/$USER/trace_$SLURM_JOB_ID
mkdir -p $W
echo "node $(hostname); scratch free: $(df -h /scratch | tail -1)"
t0=$(date +%s)
huggingface-cli download IUCompPath/TRACE --repo-type dataset \
    --include "data/TCGA-BRCA/${MODEL}/${MAG}/*" --local-dir $W --max-workers 8 > /dev/null
echo "download done in $(( $(date +%s) - t0 ))s: $(du -sh $W | cut -f1)"
python scripts/build_bags.py --parquet-dir $W/data/TCGA-BRCA/${MODEL}/${MAG} --work $W/bags \
    --out $HOME/project/data/features --name ${MODEL}_${MAG} --cap $CAP
rm -rf $W
ls -la $HOME/project/data/features
echo FEATURES_OK
