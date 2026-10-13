#!/bin/bash
# One-shot (run on the Ada gateway): (re)submit TRIDENT setup after any earlier attempt ends, then the UNI test on
# 3 slides, then all 30 frozen-slide batches (6 array tasks x 5 batches, 2 at a time). Usage: bash jobs/requeue_frozen.sh [after_job]
set -e
cd ~/project
DEP=${1:+--dependency=afterany:$1}
S=$(sbatch --parsable $DEP jobs/trident_setup.sbatch)
T=$(sbatch --parsable -J uni_test --dependency=afterok:$S --export=ALL,ENCODER=uni_v1,LIMIT=3,BATCH=0 jobs/extract_frozen.sbatch)
A=$(sbatch --parsable --array=0-5%2 --dependency=afterok:$T jobs/extract_frozen.sbatch)
echo "trident_setup $S -> uni_test $T -> extract array $A"
