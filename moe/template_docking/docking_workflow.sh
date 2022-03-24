echo "running docking workflow.."
if [[ $1 -eq 1 ]]
then
  echo "preparing input ligands..."
   /template_docking/ligand_preparation_workflow.sh input_ligands.sdf
else
  echo "Not preparing input ligands."
  # LOAD ENV
  mkdir input_ligands
  obabel input_ligands.sdf -osdf -O input_ligands/input_ligands_rdkit_quacpac_moka_split_.sdf -m  
fi

mkdir docking

mkdir tmp_aligned_3d_sdf_sanitized
for file in ligand_*.sdf; do
cat $file >> tmp_aligned_3d_sdf_sanitized/ligand_templates_for_mcs_manual.sdf
echo "\$\$\$\$" >> tmp_aligned_3d_sdf_sanitized/ligand_templates_for_mcs_manual.sdf
done

rm ligand_*.sdf
cd docking
#cp /template_docking/dock.sh .
