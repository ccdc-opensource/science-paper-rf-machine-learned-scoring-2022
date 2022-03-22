#!/bin/bash
target=$1
flexible_residues=$2
echo $target
echo $flexible_residues

file_count=$(ls ../input_ligands/input_ligands_rdkit_quacpac_moka_split_[0-9]*.sdf | wc -l)

for i in $(seq 1 $file_count)
  do
    mkdir docking_job_$i;
  done

out=$(bsub -J template_dock[1-$file_count] -q preempt -o docking_job_%I/lsf_%J_dock.out -W 01:00 -Zs "/template_docking/lsf_input_gold_docking.sh" $target $flexible_residues)
echo "waiting for all HPC jobs to finish..."
bwait -w "ended(${out//[!0-9]/})"

echo "Collecting results..."
#LOAD ENV

ccdc_roche_scoring/join_docking_results.py
chmod +rwx best_docking_solutions.sdf