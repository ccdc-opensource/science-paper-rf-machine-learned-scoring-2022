########################################################################################################################

import pandas as pd
from ccdc import descriptors, protein, io
from pathlib import Path

########################################################################################################################

ligand_path = 'tmp_aligned_3d_sdf_sanitized/single_files'
protein_path = 'apo_proteins'

distance_dict = {'strucid': [], 'hbond_distance': []}
for protein_file in Path(protein_path).glob('*dry.mol2'):
    ccdc_protein = protein.Protein.from_file(str(protein_file))
    ligand_file = str(list(Path(ligand_path).glob(f'*{ccdc_protein.identifier}*.sdf'))[0])
    ccdc_ligand = io.MoleculeReader(ligand_file)[0]
    ccdc_protein.add_ligand(ccdc_ligand)
    amide_n = [a for a in ccdc_protein.atoms if 'GLN726' in a.residue_label and a.label == 'NE2'][0]
    searcher = descriptors.MolecularDescriptors.AtomDistanceSearch(ccdc_protein)
    close_atoms = searcher.atoms_within_range(amide_n.coordinates, 3.5)
    close_atoms = [a for a in close_atoms if a.protein_atom_type == 'Ligand' and a.atomic_symbol != 'H']
    distances = []
    for close_atom in close_atoms:
        distances.append(descriptors.MolecularDescriptors.atom_distance(amide_n, close_atom))
    if distances:
        shortest_distance = min(distances)
        distance_dict['strucid'].append(ccdc_protein.identifier)
        distance_dict['hbond_distance'].append(shortest_distance)
df = pd.DataFrame(distance_dict)
df.to_csv('gln726_hbond_distances.csv', index=False)
