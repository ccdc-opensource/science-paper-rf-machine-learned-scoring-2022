#!/bin/bash

# LOAD ENV

target=$1
echo $target
i=$LSB_JOBINDEX
cd scoring_job_$i
ccdc_roche_scoring/stat_potential_inference.py -t $target -d . -l input_ligands.sdf -p ../protein.pdb