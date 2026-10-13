#!/bin/bash
# ABMIL attention heatmaps for ER-discordant (incl. ER-low) and reference patients (scripts/attention_maps.py).
# Downloads the full TRACE UNI v1 20x patch features (with coordinates) to node-local /scratch, then deletes them.
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -p u22
#SBATCH --constraint=2080ti
#SBATCH -J attention
#SBATCH -c 6
#SBATCH --gres=gpu:1
#SBATCH --mem-per-cpu=3G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%j.out
set -e
module load u22/python/3.12.4
source $HOME/venvs/whichtruth/bin/activate
cd $HOME/project
pip install -q -r requirements.txt
W=/scratch/$USER/trace_att_$SLURM_JOB_ID
mkdir -p $W
t0=$(date +%s)
python -c "from huggingface_hub import snapshot_download as s; s('IUCompPath/TRACE', repo_type='dataset', allow_patterns=['data/TCGA-BRCA/uni_v1/20x/*'], local_dir='$W', max_workers=8)"
echo "download done in $(( $(date +%s) - t0 ))s"
python -u scripts/attention_maps.py --parquet-dir $W/data/TCGA-BRCA/uni_v1/20x --out results/attention
rm -rf $W
echo ATTENTION_OK
