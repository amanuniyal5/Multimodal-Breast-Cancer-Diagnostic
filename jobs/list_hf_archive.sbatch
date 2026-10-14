#!/bin/bash
# List the members of a large .tar.gz on Hugging Face without storing it (streamed through tar -t).
# Used to check whether chadvanderbilt/TCGA_foundation_model_features (uni/BRCA_uni_features.tar.gz, 53.7 GB)
# contains frozen TS/BS slides, which would avoid the slow GDC download for 1c/1e.
#   sbatch --export=ALL,URL=<resolve url>,OUT=results/hf_archive_listing.txt jobs/list_hf_archive.sbatch
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -p u22
#SBATCH -J list_hf_archive
#SBATCH -c 2
#SBATCH --mem-per-cpu=2G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%j.out
set -e
cd $HOME/project
URL=${URL:-https://huggingface.co/datasets/chadvanderbilt/TCGA_foundation_model_features/resolve/main/uni/BRCA_uni_features.tar.gz}
OUT=${OUT:-results/hf_archive_listing.txt}
t0=$(date +%s)
curl -sL --retry 10 --retry-delay 30 "$URL" | tar -tzv > "$OUT"
echo "listed $(wc -l < $OUT) entries in $(( $(date +%s) - t0 ))s"
grep -c -- "-TS" "$OUT" | sed 's/^/TS entries: /'
grep -c -- "-BS" "$OUT" | sed 's/^/BS entries: /'
grep -c -- "-DX" "$OUT" | sed 's/^/DX entries: /'
echo LIST_OK
