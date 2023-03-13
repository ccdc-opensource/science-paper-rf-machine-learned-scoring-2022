#!/usr/bin/env python

'''
Generate binding sites by doing an MCS alignment of project data to available crystal structures.
'''

########################################################################################################################

import argparse
import pandas as pd
import subprocess as sp
from ccdc import io, protein, descriptors, entry
from pathlib import Path
from rdkit import Chem


########################################################################################################################


def parse_args():
    '''Define and parse the arguments to the script.'''
    parser = argparse.ArgumentParser(
        description=
        """
        Execute Line of sight contact scripts.
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter  # To display default values in help message.
    )

    parser.add_argument(
        '-t',
        '--target',
        help='Target name.',
        default='pde-10'
    )

    return parser.parse_args()


def _get_alignment_templates(target='pde-10'):
    structure_quality = pd.read_csv(f'{target}_structure_quality.csv')
    high_quality_strucs = structure_quality[(structure_quality['RESOLUTION'] <= 2.5) &
                                            (structure_quality['DENSITYFLAGS'] <= 5) &
                                            (structure_quality['LIGCLASHES'] == 0) &
                                            (structure_quality['LIGSYMCONTACT'] == 0)]['STRUCID'].to_list()
    assayed_structures = list((Path('') / Path(target) / Path('tmp_2d_sdf')).glob('*.sdf'))
    ligand3d_sdf = list((Path('') / Path(target) / Path('tmp_aligned_3d_sdf')).glob('*.sdf'))
    with io.EntryWriter('ligand_templates_for_mcs.sdf') as wri:
        for assayed_sdf in assayed_structures:
            pdb_id = assayed_sdf.stem.split('_')[0]
            if pdb_id not in high_quality_strucs:
                continue
            if len(pdb_id) == 4:
                ligand3D = [f for f in ligand3d_sdf if pdb_id in str(f.stem)][0]
            else:
                with io.EntryReader(str(assayed_sdf)) as rdr:
                    e = rdr[0]
                    ligname = e.identifier.split('-')[0]
                    ligand3D = [f for f in ligand3d_sdf if ligname in str(f)]
                    if len(ligand3D) < 1:
                        continue
                    else:
                        ligand3D = ligand3D[0]
            with io.EntryReader(str(ligand3D)) as rdr:
                for e in rdr:
                    e.attributes['$File'] = str(ligand3D.resolve())
                    e.attributes['STRUCID'] = pdb_id
                    wri.write(e)


def _get_pdb_ligname_df(target='pde-10'):
    pdb_files = list((Path('') / Path(target) / Path('tmp_aligned_for_MOE')).glob('*.pdb'))
    sdf_files = list((Path('') / Path(target) / Path('tmp_2d_sdf')).glob('*.sdf'))
    ligand3d_sdf = list(
        (Path('') / Path(target) / Path('tmp_aligned_3d_sdf')).glob('*.sdf'))

    structure_quality = Path(f'{target}_structure_quality.csv')
    high_quality_strucs = []
    if structure_quality.is_file():
        structure_quality = pd.read_csv(f'{target}_structure_quality.csv')
        high_quality_strucs = structure_quality[(structure_quality['RESOLUTION'] <= 2.5) &
                                                (structure_quality['DENSITYFLAGS'] <= 5) &
                                                (structure_quality['LIGCLASHES'] == 0) &
                                                (structure_quality['LIGSYMCONTACT'] == 0)]['STRUCID'].to_list()

    structure_dicts = []

    for sdf_file in sdf_files:
        strucid = sdf_file.stem.split('_')[0]
        with io.EntryReader(str(sdf_file)) as rdr:
            full_lname = rdr[0].identifier
            if full_lname.startswith('RO'):
                lname = full_lname.split('-')[0]
            else:
                lname = full_lname

        pdb_file = [f for f in pdb_files if f.stem.startswith(strucid)][0]
        template_file = None
        for f in ligand3d_sdf:
            if f.stem.startswith('RO') and f.stem.startswith(lname):
                template_file = f
                break
            elif not f.stem.startswith('RO') and f.stem.startswith(strucid):
                template_file = f
                break

        structure_dicts.append({'ligand': lname, 'full_ligand_name': full_lname, 'pdb_file': str(pdb_file),
                                'template_file': template_file, 'strucid': strucid})

    df = pd.DataFrame(structure_dicts)

    if high_quality_strucs:
        df['is_high_quality'] = df['strucid'].isin(high_quality_strucs)
    else:
        df['is_high_quality'] = True
    df.to_csv('pdb_ligand.csv', index=False)


def _sanitize_pdb_files(target='pde-10'):
    from rdkit import Chem
    pdb_files = list((Path('/rf_scoring') / Path(target) / Path('tmp_aligned_for_MOE')).glob('*.pdb'))
    sanitized_path = Path('/rf_scoring') / Path(target) / Path('tmp_aligned_for_MOE_sanitized')

    structure_quality = pd.read_csv(f'{target}_structure_quality.csv')
    high_quality_strucs = structure_quality[(structure_quality['RESOLUTION'] <= 2.5) &
                                            (structure_quality['DENSITYFLAGS'] <= 5) &
                                            (structure_quality['LIGCLASHES'] == 0) &
                                            (structure_quality['LIGSYMCONTACT'] == 0)
                                            ]['STRUCID'].to_list()

    if not sanitized_path.is_dir():
        sanitized_path.mkdir()
    for pdb_file in pdb_files:
        pdb_id = pdb_file.stem.split('_')[0].split('.')[0]
        if pdb_id not in high_quality_strucs:
            continue

        m = Chem.MolFromPDBFile(str(pdb_file), removeHs=False)
        if not m:
            print('File cannot be read by RDKit: ', pdb_file)
            with io.EntryReader(str(pdb_file)) as rdr:
                p = protein.Protein.from_entry(rdr[0])
                p.remove_unknown_atoms()
                p.kekulize()
                p.standardise_aromatic_bonds()
                for l in p.ligands:
                    if 'PO4' in l.identifier:
                        p.remove_ligand(l.identifier)
                m = Chem.MolFromMol2Block(p.to_string('mol2'))

        sanitized_file = sanitized_path / Path(pdb_file.name.replace('.pdb', '_sanitized.pdb'))
        with open(sanitized_file, 'w') as outfile:
            writer = Chem.rdmolfiles.PDBWriter(outfile)
            writer.write(m)
            writer.close()


def _get_binding_site_residues(target='pde-10'):
    from ccdc_roche.python.los_descriptors import _cut_out_binding_site_by_distance
    sanitized_path = Path('/rf_scoring') / Path(target) / Path(
        'tmp_aligned_for_MOE_sanitized')

    structure_quality = pd.read_csv(f'{target}_structure_quality.csv')
    high_quality_strucs = structure_quality[(structure_quality['RESOLUTION'] <= 2.5) &
                                            (structure_quality['DENSITYFLAGS'] <= 5) &
                                            (structure_quality['LIGCLASHES'] == 0) &
                                            (structure_quality['LIGSYMCONTACT'] == 0)]

    binding_site_residues = []
    strucids = []
    for index, row in high_quality_strucs.iterrows():
        strucid = row['STRUCID']
        pdb_file = list(Path(sanitized_path).glob(f'*{strucid}*_sanitized.pdb'))
        if pdb_file:
            pdb_file = str(pdb_file[0])
            ligname = row['HETMOLRES'][:3]
            with io.EntryReader(pdb_file) as rdr:
                for e in rdr:
                    p = protein.Protein.from_entry(e)
                    ligands = [c for c in p.components if c.atoms[-1].residue_label[:3] == ligname]
                    if not ligands:
                        continue
                    if len(ligands) != 1:
                        continue
                    else:
                        ligand = ligands[0]
                    if not 9 < len(ligand.atoms) < 500:  # kick out fragments and covalently bound ligands
                        continue
                    for l in p.ligands:
                        p.remove_ligand(l.identifier)
                    for a in p.atoms:
                        if a.residue_label[:3] == ligname:
                            p.remove_atoms([a])
                    p = _cut_out_binding_site_by_distance(p, ligand)
                    for residue in p.residues:
                        binding_site_residues.append(residue.identifier)
                        strucids.append(strucid)
    binding_site = pd.DataFrame()
    binding_site['residue_identifier'] = binding_site_residues
    binding_site['strucid'] = strucids
    binding_site = binding_site.drop_duplicates('residue_identifier')
    binding_site = binding_site['residue_identifier'].str.split(':', expand=True).drop_duplicates(subset=1).rename(
        columns={0: 'chain_label', 1: 'residue_identifier'})
    binding_site.to_csv('binding_site_residues.csv', index=False)


def _export_ligands(target='pde-10'):
    protein_out = Path('tmp_aligned_for_MOE_sanitized')
    ligand_out = Path('tmp_aligned_3d_sdf_sanitized/single_files')
    if not protein_out.is_dir():
        protein_out.mkdir()
    pdb_ligand_df = pd.read_csv('pdb_ligand.csv')
    ligand_ids = list(pdb_ligand_df['ligand'].unique())
    protein_files = Path('final_pdb_files').glob('*.pdb')
    ligands = []
    for cnt, protein_file in enumerate(protein_files):

        protein_file_str = str(protein_file)
        strucid = protein_file.stem
        # if strucid not in ['1jnhn'] and cnt != 0:
        #     continue

        ligand_name = pdb_ligand_df[pdb_ligand_df['strucid'] == strucid]['ligand'].values[0]
        ccdc_protein = protein.Protein.from_file(protein_file_str)
        ligands = []
        for l in ccdc_protein.ligands:
            for id in ['LIG', 'L0R', 'UNL', '5M9']:
                if id in l.identifier:
                    ligands.append(l)
        ligand = ligands[0].components[0].copy()
        print('removing...')
        for r in ccdc_protein.residues:
            if r.chain_identifier != 'A' and r.atoms[0].protein_atom_type == 'Amino_acid':
                ccdc_protein.remove_residue(r.identifier)
        for l in ccdc_protein.ligands:
            ccdc_protein.remove_ligand(l.identifier)
        ccdc_protein.add_ligand(ligand)
        ccdc_protein.remove_unknown_atoms()
        ccdc_protein.remove_all_metals()
        ccdc_protein.identifier = strucid

        # if cnt == 0:
        #     template_protein_chain = ccdc_protein.copy().chains[0]
        #
        # else:
        #     transformed_protein = ccdc_protein.copy()
        #     transformation_matrix = protein.Protein.ChainSuperposition().superpose(template_protein_chain, transformed_protein.chains[0])[1]
        #     for cnt, a in enumerate(transformed_protein.atoms):
        #         a.displacement_parameters = ccdc_protein.atoms[cnt].displacement_parameters
        #     ligand.transform(transformation_matrix)
        print('writing...')
        with io.MoleculeWriter(protein_out / protein_file.name) as w:
            w.write(ccdc_protein)

        ligand_entry = entry.Entry.from_molecule(ligand)
        ligand_entry.attributes['strucid'] = strucid
        ligand_entry.attributes['SRN'] = ligand_name
        ligand_entry.identifier = ligand_name
        ligands.append(ligand_entry)

        if not ligand_out.is_dir():
            ligand_out.mkdir()
        with io.EntryWriter(ligand_out / f'{strucid}_{ligand_name}.sdf') as w:
            w.write(ligand_entry)
    return


def main():
    args = parse_args()
    _export_ligands()

    # print('Getting alignment templates...')
    # _get_alignment_templates(args.target)

    # print('Assigning ligand names to PDB files...')
    # _get_pdb_ligname_df(args.target)

    # print('Sanitizing PDB files...')
    # _sanitize_pdb_files(args.target)

    # print('Getting bindingsite residues...')
    # _get_binding_site_residues(args.target)

    print('finished')


if __name__ == '__main__':
    main()
