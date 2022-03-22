#!/bin/bash

# LOAD ENV
ccdc_roche_scoring/docking.py --input_ligands input_ligands/input_ligands_rdkit_quacpac_moka_split_$i.sdf -t $target -fr=$flexible_residues
