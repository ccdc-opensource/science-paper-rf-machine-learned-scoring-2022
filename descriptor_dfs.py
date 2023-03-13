#!/usr/bin/env python

'''
Generate binding sites by doing an MCS alignment of project data to available crystal structures.
'''

########################################################################################################################

import pandas as pd
import numpy as np
from pathlib import Path
import re
import argparse
from ccdc import io
from ccdc_roche.python import los_descriptors
from ccdc_roche.python import rf_assignment


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
        '--los_home',
        help='Path to LoS home folder.',
        default='/protonate3d'
    )

    parser.add_argument(
        '-t',
        '--target',
        help='Target name.',
        default='pde-10'
    )

    parser.add_argument(
        '--aligned_ligands',
        help='SDF file with aligned ligands',
        default=False
    )

    parser.add_argument(
        '--docking_dir',
        help='Directory with docked ligands',
        default='.'
    )

    return parser.parse_args()


def _group_by_template(aligned_molecules, project_home):
    '''
    Assign aligned ligands to the template's strucid.
    :param aligned_molecules:
    :param pdb_ligand_file:
    :return:
    '''
    conformer_count = {}
    pdb_ligand_df = pd.read_csv(Path(project_home) / Path('pdb_ligand.csv'))
    with io.EntryReader(aligned_molecules) as rdr:
        group_by_template = {}
        for ligand in rdr.entries():
            ligand_query = ligand.attributes['Query'].split('_')[0]
            if 'RO' in ligand_query:
                ro_number = ligand_query.split('RO')[1]
                ro_number = re.split(r'[A-Z]', ro_number)[0]
                ligand_query = 'RO' + ro_number
            pdb = pdb_ligand_df[pdb_ligand_df['ligand'] == ligand_query]['pdb'].values[0]
            pdb = str(Path(project_home) / Path('tmp_aligned_for_MOE_sanitized') / Path(pdb).name.replace('.pdb',
                                                                                                          '_sanitized.pdb'))

            if ligand.identifier in conformer_count.keys():
                ligand.attributes['conformer'] = conformer_count[ligand.identifier]
                conformer_count[ligand.identifier] += 1
            else:
                ligand.attributes['conformer'] = 0
                conformer_count[ligand.identifier] = 1

            if pdb in group_by_template.keys():
                group_by_template[pdb] = group_by_template[pdb] + [ligand]
            else:
                group_by_template[pdb] = [ligand]
    return group_by_template


def extend_df(df, attr, strucid=''):
    '''
    Add additional columns to Data Frame
    :param df: Dataframe to be extended
    :param attr: Dictionary with additional columns.
    :return:
    '''
    df['identifier'] = strucid
    for a in attr:
        df[a] = attr[a]
    return df


def return_contact_df_from_docking_file(docked_ligand_file, target_home='.', strucid='', pdb_file=None):
    gold_conf = str(list(Path(docked_ligand_file).parents[0].glob(f'gold_*.conf'))[0])

    if pdb_file is None:
        pdb_file = list((Path('../..') / Path('tmp_aligned_for_MOE_sanitized')).glob(f'{strucid}*.pdb'))
        if pdb_file:
            pdb_file = pdb_file[0]
        else:
            pdb_file = list((Path(target_home)).glob(f'{strucid}*prot3d.pdb'))[0]

    assigner = rf_assignment.RfAssigner(str(docked_ligand_file), gold=gold_conf, only_binding_site=False,
                                        intramolecular_protein_contacts=False, bfactor=pdb_file,
                                        interaction_cutoff=1.5, pdb_file=pdb_file)

    assigner.ligand_file = str(docked_ligand_file)
    describer = assigner.describer
    mol_contact_df = assigner.rf_assignments
    mol_contact_df = mol_contact_df.drop(
        columns=[c for c in mol_contact_df.columns if 'Gold.Protein' in c and c != 'Gold.PLP.Chemscore.Protein.Energy'])
    mol_contact_df = mol_contact_df.join(
        assigner.describer.bfactor_df[['normalised_atom_label', 'relative_bfactor']].set_index(
            'normalised_atom_label'), on='los_atom_label')

    rf_count_df = los_descriptors.rf_count_df(mol_contact_df, describer.csd_ligand)
    rf_count_df['rotatable_bonds_num'] = assigner.rotatable_bonds_num
    rf_count_df['frozen_bonds_num'] = assigner.frozen_bonds_num
    rf_score = np.log(
        mol_contact_df[(mol_contact_df['rf_total'] > 0) & (mol_contact_df['is_intramolecular'] == False)][
            'rf_total']).sum()
    rf_count_df['rf_score'] = rf_score
    rf_count_df = extend_df(rf_count_df, describer.csd_ligand_entry.attributes, strucid)

    return rf_count_df, mol_contact_df


def return_contact_df(docking_dir=''):
    if 'decoy' in str(Path(docking_dir).absolute()):
        docked_structures = list(Path(docking_dir).glob(f'*_*/best_decoy_*.sdf'))
    else:
        docked_structures = list(Path(docking_dir).glob(f'*_*/best_soln_*.sdf'))

    contacts_df_list = []
    rf_count_df_list = []

    for docked_ligand_file in docked_structures:
        strucid = docked_ligand_file.stem.split('_')[2]
        rf_count_df, mol_contact_df = return_contact_df_from_docking_file(docked_ligand_file, strucid=strucid)
        rf_count_df_list.append(rf_count_df)
        contacts_df_list.append(mol_contact_df)

    rf_count_df = pd.concat(rf_count_df_list, ignore_index=True)
    contact_df = pd.concat(contacts_df_list, ignore_index=True)

    return contact_df, rf_count_df


def main():
    args = parse_args()
    if args.docking_dir:
        contact_df, rf_count_df = return_contact_df(args.docking_dir)

        rf_count_df.columns = [str(c) for c in rf_count_df.columns]
        rf_count_df = rf_count_df.drop(
            columns=[c for c in rf_count_df.columns if 'Gold.Protein' in c or 'Gold.Chemscore.Hbonds' in c])
        rf_count_df.to_csv(str(Path(args.docking_dir) / Path('docked_rf_count_df.csv')), index=False)

        contact_df.columns = [str(c) for c in contact_df.columns]
        contact_df = contact_df.drop(
            columns=[c for c in contact_df.columns if 'Gold.Protein' in c or 'Gold.Chemscore.Hbonds' in c])
        contact_df.to_csv(str(Path(args.docking_dir) / Path('docked_contact_df.csv')), index=False)

    print('finished')


if __name__ == '__main__':
    main()
