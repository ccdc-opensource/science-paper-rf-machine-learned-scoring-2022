#!/bin/bash
echo "running Scoring workflow.."
# LOAD ENV

echo "Scoring..."
target=$1

file_count=$(ls -d scoring_job_* | wc -l)
out=$(bsub -J rf_scoring[1-$file_count] -q preempt -o lsf_%J_scoring.out -W 01:00 -Zs "/rf_scoring/lsf_input_scoring.sh" $target)
echo "waiting for all HPC jobs to finish..."
bwait -w "ended(${out//[!0-9]/})"

cat scoring_job_*/rescored_ligands.sdf >> rescored_ligands.sdf
