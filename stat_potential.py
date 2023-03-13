#!/usr/bin/env python

########################################################################################################################

import numpy as np
import pandas as pd

np.seterr(all="ignore")
from scipy.optimize import least_squares
from scipy.stats import spearmanr
from pathlib import Path
import itertools
import argparse
import json


########################################################################################################################


def parse_args():
    '''Define and parse the arguments to the script.'''
    parser = argparse.ArgumentParser(
        description=
        """
        Optimize parameters for RF-PLP scoring function.
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter  # To display default values in help message.
    )

    parser.add_argument(
        '-t',
        '--target',
        help='Target name.',
        default='pde-10'
    )

    parser.add_argument(
        '-r',
        '--rescore',
        help='Rescore ligand poses.',
        action='store_true',
        default=False
    )

    parser.add_argument(
        '-b',
        '--benchmark',
        help='Run benchmark test.',
        action='store_true',
        default=False
    )

    parser.add_argument(
        '-df',
        '--dataframe',
        help='DataFrame containing split data',
        default=False
    )

    return parser.parse_args()


def hbond_compensation(ligand_contacts):
    '''
    compensate hbonds for implicit hydrogen treatment
    :param ligand_contacts DataFrame:
    :return ligand_contacts DataFrame:
    '''

    # halogens and sulfur making contact to "mix" or "Water" atom types receive desolvation penalty
    ligand_contacts.loc[
        (ligand_contacts['ligand_atom_symbol']).isin(['S', 'Cl', 'F', 'I']) &
        # (ligand_contacts['rf_total'] < 1 & ligand_contacts['is_hbond']) &
        (ligand_contacts['protein_atom_type'].isin(['O_ali_mix', 'O_pi_mix', 'Water'])),
        'rf_total'] = ligand_contacts['rf_total'] - 0.5

    ligand_contacts.loc[
        (ligand_contacts['los_atom_symbol']).isin(['S']) &
        (ligand_contacts['ligand_atom_type'].str.contains('alcohol')),
        'rf_total'] = ligand_contacts['rf_total'] - 0.5

    ligand_contacts.loc[
        (ligand_contacts['ligand_atom_symbol']).isin(['S', 'Cl', 'F', 'I']) & (
                ligand_contacts['rf_total'] < 1 & ligand_contacts['protein_atom_type'].isin(
            ['O_ali_mix', 'O_pi_mix', 'Water'])), 'interaction_type'] = 'desolvation'

    ligand_contacts.loc[
        (ligand_contacts['ligand_atom_symbol']).isin(['S', 'Cl', 'F', 'I']) | (ligand_contacts['los_atom_symbol']).isin(
            ['S', 'Cl', 'F', 'I']), 'is_hbond'] = 0

    # Classic hbond donors/acceptors that do not make an explicit hbond receive a penalty and are recategorised
    ligand_contacts.loc[(ligand_contacts['interaction_type'] == 'hbond_classic') &
                        (ligand_contacts['is_hbond'] == False) &
                        (ligand_contacts['protein_atom_type'].isin(['O_ali_mix', 'O_pi_mix', 'Water'])), 'rf_total'
    ] = ligand_contacts['rf_total'] - 1
    ligand_contacts.loc[(ligand_contacts['is_hbond'] == False) &
                        (ligand_contacts['interaction_type'] == 'hbond_classic') &
                        (ligand_contacts['is_hbond_mismatch'] == False), 'interaction_type'] = 'desolvation'

    # Remaining hbonds are classic hbonds
    ligand_contacts.loc[(ligand_contacts['is_hbond'] == True) & (
            ligand_contacts['ligand_atom_symbol'] != 'C'), 'interaction_type'] = 'hbond_classic'

    # Polarized C-H bonds make weak hbonds
    ligand_contacts.loc[(ligand_contacts['is_hbond'] == True) &
                        (ligand_contacts['ligand_atom_symbol'] == 'C') &
                        ((ligand_contacts['ligand_atom_type'].str.contains('polarized')) |
                         (ligand_contacts['ligand_atom_type'].str.contains('aromatic'))),
                        'interaction_type'] = 'hbond_weak'
    ligand_contacts.loc[(ligand_contacts['is_hbond'] == True) &
                        (ligand_contacts['los_atom_symbol'] == 'C'), 'interaction_type'] = 'hbond_weak'

    # Polarized C-H...Donor are hbond-mismatch
    ligand_contacts.loc[(ligand_contacts['is_hbond'] == False) &
                        (ligand_contacts['ligand_atom_symbol'] == 'C') &
                        (ligand_contacts['ligand_atom_type'].str.contains('aromatic')) &
                        (ligand_contacts['rf_total'] < 0.8) &
                        (ligand_contacts['protein_atom_type'].str.contains('_don')) &
                        (ligand_contacts['protein_atom_type'].str.contains('_pi')) &
                        (ligand_contacts['protein_h'] < 1) & (ligand_contacts['ligand_h'] < 1),
                        'interaction_type'] = 'hbond_mismatch'

    # Polarized C-H...Donor are hbond-mismatch
    ligand_contacts.loc[(ligand_contacts['is_hbond'] == False) &
                        (ligand_contacts['ligand_atom_symbol'] == 'C') &
                        (ligand_contacts['ligand_atom_type'].str.contains('aromatic')) &
                        (ligand_contacts['rf_total'] < 0.8) &
                        (ligand_contacts['ligand_h'] < 1) &
                        (ligand_contacts['protein_atom_type'].str.contains('_mix')) &
                        (ligand_contacts['los_atom_symbol'] == 'O'),
                        'interaction_type'] = 'hbond_mismatch'

    # classic hbonds are shorter than vdw distance, see: The Geometry of the N-H...O--C Hydrogen Bond.
    # 3.* Hydrogen-Bond Distances and Angles, BY ROBIN TAYLOR, OLGA KENNARD AND WERNER VERSICHEL

    ligand_contacts.loc[ligand_contacts['interaction_type'] == 'hbond_classic', 'rm'] = ligand_contacts['rm'] - 0.3

    ligand_contacts.loc[(ligand_contacts['interaction_type'] == 'hbond_classic') &
                        (ligand_contacts['protein_atom_type'].isin(['O_ali_mix', 'O_pi_mix', 'Water'])),
                        'rf_total'] = ligand_contacts['rf_total'] + 1

    ligand_contacts.loc[(ligand_contacts['interaction_type'] == 'hbond_classic') &
                        (ligand_contacts['ligand_atom_type'].str.contains('alcohol')), 'rf_total'] = ligand_contacts[
                                                                                                         'rf_total'] + 0.8
    ligand_contacts.loc[(ligand_contacts['interaction_type'] == 'hbond_weak') &
                        (ligand_contacts['is_hbond'] == True) &
                        (ligand_contacts['protein_atom_type'].isin(['O_ali_mix', 'O_pi_mix', 'Water'])),
                        'rf_total'] = ligand_contacts['rf_total'] + 0.3

    return ligand_contacts


def turn_off_secondary_clash_contacts(ligand_contacts):
    '''Avoid compensating for clashes by secondary contacts'''
    ligand_contacts['is_clash'] = ((ligand_contacts['interaction_type'] != 'hbond_classic') &
                                   (ligand_contacts['vdw_distance'] < -0.5)) | \
                                  ((ligand_contacts['interaction_type'] == 'hbond_classic') &
                                   (ligand_contacts['vdw_distance'] < -0.7) |
                                   (ligand_contacts['interaction_type'].isin(
                                       ['electrostatic_repulsion', 'hbond_mismatch'])) &
                                   (ligand_contacts['vdw_distance'] < -0.3)
                                   )
    secondary_clash_contact_indices = []
    for ligand_file, df in ligand_contacts.groupby('ligand_file'):
        # Get the indices of the ligand atoms that clash
        clashing_atom_indices = df[df['is_clash']]['ligand_atom_index']
        # clashing_atom_indices = clashing_atom_indices.append(df[(df['is_hbond']==False) & (df['vdw_distance'] < -0.1)]['ligand_atom_index'])

        # Get the df indices of the contacts that clash
        clashing_df_indices = df[df['is_clash']].index
        # clashing_df_indices = clashing_df_indices.append(df[(df['is_hbond'] == False) & (df['vdw_distance'] < -0.1)].index)

        # Get the atom indices of the secondary contacts that do no clash
        clashing_atom_contact_indices = df[(df['ligand_atom_index'].isin(clashing_atom_indices)) &
                                           (df.index.isin(clashing_df_indices) == False)].index

        secondary_clash_contact_indices.extend(clashing_atom_contact_indices)

        # # set water bfactor to binding site minimum
        # water_rows = df[df['protein_atom_type'] == 'Water'].index
        # min_bfactor = df['relative_bfactor'].min()
        # ligand_contacts.loc[water_rows, 'relative_bfactor'] = min_bfactor
    ligand_contacts['clash_count'] = \
        ligand_contacts[ligand_contacts['is_clash'] == True].groupby(['ligand_file', 'ligand_atom_index'])[
            ['ligand_atom_index']].transform('count')
    ligand_contacts['clash_count'] = ligand_contacts['clash_count'].fillna(0)
    ligand_contacts.loc[ligand_contacts['clash_count'] > 3, 'clash_count'] = 3
    ligand_contacts = ligand_contacts.drop(index=secondary_clash_contact_indices)

    return ligand_contacts


def prepare_contact_df(ligand_contacts, count_df):
    ligand_contacts['rf_total'] = ligand_contacts['rf_total'].replace(0, 1)
    ligand_contacts['interaction_type'] = ligand_contacts['interaction_type'].fillna('hydrophobic')

    radii = {'C': 1.7, 'N': 1.55, 'O': 1.52, 'S': 1.8, 'Cl': 1.75, 'F': 1.47, 'Br': 1.85, 'I': 1.98}
    for symbol in radii:
        ligand_contacts.loc[ligand_contacts['los_atom_symbol'] == symbol, 'los_atom_radius'] = radii[symbol]
        ligand_contacts.loc[ligand_contacts['ligand_atom_symbol'] == symbol, 'ligand_atom_radius'] = radii[symbol]
    ligand_contacts['rm'] = ligand_contacts['ligand_atom_radius'] + ligand_contacts['los_atom_radius']

    ligand_contacts = hbond_compensation(ligand_contacts)
    ligand_contacts = turn_off_secondary_clash_contacts(ligand_contacts)
    ligand_contacts.loc[ligand_contacts['interaction_type'].isin(['electrostatic_repulsion', 'hbond_mismatch']),
                        'interaction_type'] = 'repulsive'

    for index, row in count_df[count_df['ligand_file'].isin(ligand_contacts['ligand_file'])].iterrows():
        ligand_file = row['ligand_file']
        pc_frozen_bonds = 0
        if row['rotatable_bonds_num'] > 0:
            pc_frozen_bonds = row['frozen_bonds_num'] / row['rotatable_bonds_num']
        ligand_contacts.loc[ligand_contacts['ligand_file'] == ligand_file, 'pc_frozen_bonds'] = pc_frozen_bonds
    ligand_contacts.loc[ligand_contacts['rf_total'] <= 0, 'rf_total'] = 0.1
    return ligand_contacts


def process_df(strain=['ANI_Strain'], task='pde10_h_pic50'):
    '''
    dG = a*Sum(-ln(RF)) + b*Sum(-ln(RF)) + ... + x * Internal + y * Clash
    rf_score: favorable interactions have negative values, unfavorable have positive values.
    :return:
    '''
    if not Path('best_soln_contacts.csv').is_file():
        contact_df = pd.concat(map(pd.read_csv, Path('..').glob('docking_job_*/docked_contact_df.csv')),
                               ignore_index=True)
        contact_df = contact_df[contact_df['is_intramolecular'] == False]
        contact_df = contact_df.drop(columns=[c for c in contact_df.columns if 'Gold.Protein' in c])
        contact_df = contact_df.drop(columns=[c for c in contact_df.columns if c in ['Gold.Version',
                                                                                     'Gold.Chemscore.Hbonds',
                                                                                     'Gold.Score', 'HAC', 'NETCHARGE',
                                                                                     'ABUNDANCE']])
        # contact_df = contact_df[contact_df['is_primary'] == True]
        contact_df.to_csv('best_soln_contacts.csv', index=False)

    count_df_file = Path('../docked_rf_count_best_docked.csv')
    count_df = pd.read_csv(count_df_file)
    count_df = count_df[['ligand_file', 'frozen_bonds_num', 'rotatable_bonds_num']]

    ligand_contacts = pd.read_csv('best_soln_contacts.csv')
    ligand_contacts = ligand_contacts[ligand_contacts['ligand_file'].isin(count_df['ligand_file'])]
    ligand_contacts = ligand_contacts[ligand_contacts[task].isna() == False]

    if 'ANI_Strain' in strain:
        ani_df = pd.read_csv('../all_docked_ani.csv')
        ani_df = ani_df[['SRN', 'ANI_Strain']]
        ligand_contacts = ligand_contacts.join(ani_df.set_index('SRN'), on='SRN')
        ligand_contacts = ligand_contacts[ligand_contacts['ANI_Strain'].isna() == False]

    ligand_contacts = prepare_contact_df(ligand_contacts, count_df)

    ligand_contacts.to_csv('best_soln_contacts_processed.csv', index=False)
    return ligand_contacts


class PlpScoring(object):
    def __init__(self, task, strain=[], series='', gold=True):
        '''

        :param task:
        :param strain:
        :param series:
        :param gold:
        '''
        self.task = task
        self.gold = gold
        self.series = series
        self.interaction_types = ['pi', 'hydrophobic', 'hbond_weak', 'multipolar', 'repulsive', 'desolvation',
                                  'hbond_classic', 'ionic']
        params = ['w']
        interaction_type_params = [a + '_' + b for a, b in itertools.product(self.interaction_types, params)]

        fav_interaction_type_groups = ['favorable_lip', 'favorable_hb']
        params = ['a', 'b', 'c', 'd', 'e', 'f']
        fav_params = [a + '_' + b for a, b in itertools.product(fav_interaction_type_groups, params)]

        unfav_interaction_type_groups = ['unfavorable']
        params = ['a', 'c', 'f']
        params = interaction_type_params + fav_params + [a + '_' + b for a, b in
                                                         itertools.product(unfav_interaction_type_groups, params)]

        self.interaction_type_groups = fav_interaction_type_groups + unfav_interaction_type_groups

        bounds_dict = {}
        initial_guess = {}
        for p in params:
            if p.endswith('_w'):
                bounds_dict[p] = (-5, 0)
                initial_guess[p] = -1
            if p.endswith('_a'):
                bounds_dict[p] = (0.05, 0.10)
                initial_guess[p] = 0.1
            if p.endswith('unfavorable_a'):
                bounds_dict[p] = (0.05, 0.30)
                initial_guess[p] = 0.1
            if p.endswith('_b'):
                bounds_dict[p] = (0.1, 0.4)
                initial_guess[p] = 0.1
            if p.endswith('_c'):
                bounds_dict[p] = (0.05, 0.8)
                initial_guess[p] = 0.1
            if p.endswith('unfavorable_c'):
                bounds_dict[p] = (0.5, 0.7)
                initial_guess[p] = 0.1
            if p.endswith('_d'):
                bounds_dict[p] = (0.01, 0.5)
                initial_guess[p] = 0.1
            if p.endswith('_e'):
                bounds_dict[p] = (-2, 2)
                initial_guess[p] = 0.1
            if p.endswith('unfavorable_e'):
                bounds_dict[p] = (0, 0.3)
                initial_guess[p] = 0.1

        # Parameters for receptor depths scaling similar as in Boyle, J. Chem. Inf. Model. 2008, 48, 1269–1278
        # bounds_dict['hbond_S'] = (-5, 0)
        bounds_dict['hbond_rho1'] = (0, 70)
        bounds_dict['hbond_rho2'] = (10, 60)
        # initial_guess['hbond_S'] = -0.5
        initial_guess['hbond_rho1'] = 10
        initial_guess['hbond_rho2'] = 10

        # bounds_dict['lipophilic_S'] = (0, 2)
        bounds_dict['lipophilic_rho1'] = (0, 70)
        bounds_dict['lipophilic_rho2'] = (10, 60)
        # initial_guess['lipophilic_S'] = 0.5
        initial_guess['lipophilic_rho1'] = 10
        initial_guess['lipophilic_rho2'] = 10

        # Scale by relative b-factor
        # bounds_dict['bfactor_S'] = (1, 2)
        bounds_dict['bfactor_rho1'] = (-0.6, 0.4)
        bounds_dict['bfactor_rho2'] = (0, 0.3)
        # initial_guess['bfactor_S'] = 1
        initial_guess['bfactor_rho1'] = 10
        initial_guess['bfactor_rho2'] = 10

        for cnt, strain_label in enumerate(strain):
            bounds_dict[f'strain_weight_{cnt}'] = (-2, 0)
            initial_guess[f'strain_weight_{cnt}'] = -0.1
        self.strain = strain

        bounds_dict['rotation_entropy_w'] = (-5, 0)
        initial_guess['rotation_entropy_w'] = -0.1

        bounds_dict['clash_factor'] = (0.2, 1)
        initial_guess['clash_factor'] = 1

        self.initial_guess = initial_guess
        self.bounds_dict = bounds_dict

    def plp_scoring_function(self, x, ligand_contact_df, return_contact_df=False):
        '''
        :param x:
        :param ligand_contact_df:
        :return:
        '''

        params_dict = {p: x[cnt] for cnt, p in enumerate(self.bounds_dict.keys())}
        clash_factor = params_dict['clash_factor']

        params_dict['lipophilic_rho2'] = params_dict['lipophilic_rho2'] + params_dict['lipophilic_rho1']
        params_dict['hbond_rho2'] = params_dict['hbond_rho2'] + params_dict['hbond_rho1']
        params_dict['bfactor_rho2'] = params_dict['bfactor_rho2'] + params_dict['bfactor_rho1']

        int_dict = {'favorable_lip': ['pi', 'hydrophobic', 'hbond_weak', 'multipolar', 'desolvation'],
                    'favorable_hb': ['hbond_classic', 'ionic'],
                    'unfavorable': ['repulsive']}
        ligand_contact_df['plp_energy'] = 0

        A = pd.Series(dtype=float)
        B = pd.Series(dtype=float)
        C = pd.Series(dtype=float)
        D = pd.Series(dtype=float)
        E = pd.Series(dtype=float)

        A_unfav = pd.Series(dtype=float)
        B_unfav = pd.Series(dtype=float)
        F_unfav = pd.Series(dtype=float)

        for interaction_ in int_dict:
            if interaction_ == 'unfavorable':
                a = params_dict[interaction_ + '_a']
                c = params_dict[interaction_ + '_c']

                interaction_indices = ligand_contact_df[
                    ligand_contact_df['interaction_type'].isin(int_dict[interaction_])].index
                rm = ligand_contact_df.loc[interaction_indices]['rm']
                A_unfav = A_unfav.append(rm - a)
                B_unfav = B_unfav.append(rm + c)

                # F_unfav is positive for RF_total < 1
                _c = -(np.log(ligand_contact_df.loc[interaction_indices]['rf_total']))
                F_unfav = F_unfav.append(_c)

            else:
                b = params_dict[interaction_ + '_b']
                a = params_dict[interaction_ + '_a'] + b
                c = params_dict[interaction_ + '_c']
                d = params_dict[interaction_ + '_d'] + c
                e = params_dict[interaction_ + '_e']

                interaction_indices = ligand_contact_df[
                    ligand_contact_df['interaction_type'].isin(int_dict[interaction_])].index
                rm = ligand_contact_df.loc[interaction_indices]['rm']
                A = A.append(rm - a)
                B = B.append(rm - b)
                C = C.append(rm + c)
                D = D.append(rm + d)
                E = E.append(-(np.log(ligand_contact_df.loc[interaction_indices]['rf_total'])) + e)

        A = A.sort_index()
        B = B.sort_index()
        C = C.sort_index()
        D = D.sort_index()
        E = E.sort_index()

        A_unfav = A_unfav.sort_index()
        B_unfav = B_unfav.sort_index()
        F_unfav = F_unfav.sort_index()

        # negative values represent attraction, positive values represent repulsion
        favorable_contacts = ligand_contact_df[
            ligand_contact_df['interaction_type'].isin(int_dict['unfavorable']) == False].copy()
        distances = favorable_contacts['distance'].sort_index()
        clash_counts = favorable_contacts['clash_count'].sort_index()
        favorable_contacts.loc[(distances < A), 'plp_energy'] = clash_factor * (clash_counts + 1) * (
                (A - distances) / A) + clash_factor
        favorable_contacts.loc[(distances >= A) & (distances < B), 'plp_energy'] = E * (distances - A) / (B - A)
        favorable_contacts.loc[(distances >= B) & (distances < C), 'plp_energy'] = E
        favorable_contacts.loc[(distances >= C) & (distances < D), 'plp_energy'] = E * (D - distances) / (D - C)

        unfavorable_contacts = ligand_contact_df[
            ligand_contact_df['interaction_type'].isin(int_dict['unfavorable']) == True].copy()
        clash_counts = unfavorable_contacts['clash_count'].sort_index()
        distances = unfavorable_contacts['distance'].sort_index()
        unfavorable_contacts.loc[(distances < A_unfav), 'plp_energy'] = clash_factor * (clash_counts + 1) * (
                (A_unfav - distances) / A_unfav) + clash_factor
        unfavorable_contacts.loc[(distances >= A_unfav) & (distances < B_unfav), 'plp_energy'] = -F_unfav * (
                distances - A_unfav) / (B_unfav - A_unfav) + F_unfav

        ligand_contact_df = pd.concat([favorable_contacts, unfavorable_contacts])

        for interaction_type, ligand_interaction_df in ligand_contact_df.groupby('interaction_type'):

            if interaction_type in ['hbond_classic', 'ionic']:
                buriedness = ligand_interaction_df['los_atom_buriedness']

                ligand_interaction_df.loc[buriedness < params_dict['hbond_rho1'], 'plp_energy'] = 0

                ligand_interaction_df.loc[(params_dict['hbond_rho1'] <= buriedness) &
                                          (buriedness < params_dict['hbond_rho2']), 'plp_energy'
                ] = ligand_interaction_df['plp_energy'] * (buriedness - params_dict['hbond_rho1']) / (
                        params_dict['hbond_rho2'] - params_dict['hbond_rho1'])

                ligand_interaction_df.loc[buriedness >= params_dict['hbond_rho2'], 'plp_energy'] = \
                    ligand_interaction_df['plp_energy']

            if interaction_type in ['pi', 'hydrophobic', 'multipolar', 'hbond_weak', 'desolvation']:
                buriedness = ligand_interaction_df['los_atom_buriedness']

                ligand_interaction_df.loc[buriedness < params_dict['lipophilic_rho1'], 'plp_energy'] = 0
                ligand_interaction_df.loc[
                    (params_dict['lipophilic_rho1'] <= buriedness) &
                    (buriedness < params_dict['lipophilic_rho2']), 'plp_energy'
                ] = ligand_interaction_df['plp_energy'] * (buriedness - params_dict['lipophilic_rho1']) / (
                        params_dict['lipophilic_rho2'] - params_dict['lipophilic_rho1'])

            ligand_contact_df.loc[ligand_interaction_df.index, 'plp_energy'] = ligand_interaction_df['plp_energy']

        # bfactor scaling
        bfactors = ligand_contact_df['relative_bfactor']

        ligand_contact_df.loc[
            (params_dict['bfactor_rho1'] <= bfactors) &
            (bfactors < params_dict['bfactor_rho2']), 'plp_energy'
        ] = ligand_contact_df['plp_energy'] * (params_dict['bfactor_rho2'] - bfactors) / (
                params_dict['bfactor_rho2'] - params_dict['bfactor_rho1'])

        ligand_contact_df.loc[bfactors >= params_dict['bfactor_rho2'], 'plp_energy'] = 0

        clash_contact_indices = list(unfavorable_contacts[unfavorable_contacts['distance'].lt(A_unfav)].index) + list(
            favorable_contacts[favorable_contacts['distance'].lt(A)].index)

        for interaction_type, interaction_df in ligand_contact_df[
            ligand_contact_df.index.isin(clash_contact_indices) == False].groupby('interaction_type'):
            ligand_contact_df.loc[interaction_df.index, 'plp_energy'] = interaction_df['plp_energy'] * params_dict[
                interaction_type + '_w']
        ligand_contact_df.loc[ligand_contact_df.index.isin(clash_contact_indices), 'plp_energy'] = -ligand_contact_df[
            'plp_energy']
        column_labels = ['ligand_file', self.task, 'pc_frozen_bonds', 'cluster'] + self.strain
        if self.gold:
            column_labels.append('Gold.PLP.Fitness')

        prediction_df = ligand_contact_df[column_labels].drop_duplicates('ligand_file')

        predicted_activities = ligand_contact_df.groupby(['ligand_file'])['plp_energy'].sum()
        prediction_df = prediction_df.join(predicted_activities, on='ligand_file')

        prediction_df['strain_penalty'] = 0
        for cnt, strain_label in enumerate(self.strain):
            prediction_df = prediction_df.astype({strain_label: float})
            if strain_label == 'ANI_Strain':
                prediction_df['strain_penalty'] = params_dict[f'strain_weight_{cnt}'] * np.log(
                    prediction_df[strain_label])
                prediction_df['strain_penalty'] = prediction_df['strain_penalty'].replace([np.inf, -np.inf], 0)
            else:
                prediction_df['strain_penalty'] = prediction_df['strain_penalty'] + params_dict[
                    f'strain_weight_{cnt}'] * prediction_df[strain_label]

        prediction_df['rotation_entropy_penalty'] = params_dict['rotation_entropy_w'] * prediction_df['pc_frozen_bonds']
        prediction_df['score'] = prediction_df['plp_energy'] + prediction_df['strain_penalty'] + prediction_df[
            'rotation_entropy_penalty']
        prediction_df.loc[prediction_df['score'] < 3, 'score'] = 3
        if return_contact_df:
            return prediction_df, ligand_contact_df

        return prediction_df

    def optimization_function(self, x, scores_df=pd.DataFrame()):
        prediction_df = self.plp_scoring_function(x, ligand_contact_df=scores_df)
        if self.minimizer == 'minimize':
            error = (prediction_df[self.task] - prediction_df['score']).abs().median()
            print(error)
            return error
        elif self.minimizer == 'least_squares':
            error = (prediction_df[self.task] - prediction_df['score'])
            print(error.abs().median())
            return error

    def optimize_weights(self, scores_df):

        self.initial_guess = {"pi_w": -0.06329918369899717, "hydrophobic_w": -0.03679486613866424,
                              "hbond_weak_w": -0.043416752578643,
                              "multipolar_w": -0.019966358787951787, "repulsive_w": -0.04166345024299469,
                              "desolvation_w": -0.009318217735359038, "hbond_classic_w": -0.2285049973821565,
                              "ionic_w": -0.1765648797949396,
                              "favorable_lip_a": 0.09999493880055929, "favorable_lip_b": 0.39996346524216725,
                              "favorable_lip_c": 0.7999992895585842, "favorable_lip_d": 0.49999927117641335,
                              "favorable_lip_e": -1.601062326025244, "favorable_hb_a": 0.09351018567138517,
                              "favorable_hb_b": 0.26874110245376265, "favorable_hb_c": 0.4710463226463769,
                              "favorable_hb_d": 0.2947095948102579, "favorable_hb_e": -0.40114233297001367,
                              "unfavorable_a": 0.29990714415021325, "unfavorable_c": 0.5005102385721926,
                              "hbond_rho1": 50.45643860461278,
                              "hbond_rho2": 57.452535372621725, "lipophilic_rho1": 19.205766641846076,
                              "lipophilic_rho2": 21.41581509284962,
                              "bfactor_rho1": 0.1618446393425923, "bfactor_rho2": 0.1639701197217107,
                              "strain_weight_0": -0.0007038247138676942, "strain_weight_1": -0.11715414788931568,
                              "strain_weight_2": -5.154179015906084e-09, "rotation_entropy_w": -3.527789958418665e-06,
                              "clash_factor": 0.5011126435789334}

        self.initial_guess = self.initial_guess.values()

        scores_df = scores_df.drop(
            columns=[c for c in scores_df.columns if 'Gold.' in c and c not in self.strain + ['Gold.PLP.Fitness']])

        from scipy.optimize import minimize
        self.minimizer = 'minimize'
        ret = minimize(self.optimization_function, method='L-BFGS-B', x0=list(self.initial_guess), args=(scores_df,),
                       bounds=(self.bounds_dict.values()), tol=1e-03)  # , tol=1e-05

        print('Least-squares minmization....')
        self.minimizer = 'least_squares'
        ret = least_squares(self.optimization_function, x0=ret.x, args=(scores_df,),
                            bounds=(np.array(list(self.bounds_dict.values()))[:, 0],
                                    np.array(list(self.bounds_dict.values()))[:, 1]), ftol=1e-03)  # , ftol=1e-05

        weight_dict = {p: ret.x[cnt] for cnt, p in enumerate(self.bounds_dict.keys())}
        with open(f'weights_{self.series}.json', "w") as out_file:
            json.dump(weight_dict, out_file)

        print(ret)
        return ret.x, weight_dict


def main():
    args = parse_args()
    if args.target == 'pde-10':
        task = 'pde10_h_pic50'

    strain = ['Gold.PLP.ligand.torsion', 'Gold.PLP.ligand.clash', 'Gold.PLP.Chemscore.Protein.Energy']

    contact_scores_files = Path('best_soln_contacts_processed.csv')
    if not args.rescore:
        print('')
        process_df(strain=strain, task=task)

    elif args.rescore and not contact_scores_files.is_file():
        # RESCORE
        print('Rescoring poses...')
        scores_df = pd.concat([pd.read_csv(f) for f in Path('..').glob('docking_job_*/rescored_rf_contact.csv')],
                              ignore_index=True)
        scores_df = scores_df[scores_df['RMSD_to_mcs'] < 1.5]
        ligand_files_remove = []
        for index, df in scores_df.sort_values(by='rescore', ascending=False).groupby('SRN'):
            ligand_files = df['ligand_file'].unique()
            ligand_files_remove.extend(ligand_files[1:])
        scores_df = scores_df[scores_df['ligand_file'].isin(ligand_files_remove) == False]
        scores_df = scores_df[scores_df['rel_mcs_size'] >= 0.5]
        scores_df = scores_df.drop(columns=[c for c in scores_df.columns if 'Gold.Protein' in c])
        scores_df = scores_df.drop(columns=[c for c in scores_df.columns if c in ['Gold.Version',
                                                                                  'Gold.Chemscore.Hbonds',
                                                                                  'Gold.Score']])
        scores_df.to_csv(contact_scores_files, index=False)

    scores_df = pd.read_csv(contact_scores_files)

    series_file = Path('../../series_assignment.csv')
    if series_file.is_file() and not args.dataframe:
        series_df = pd.read_csv(series_file)
        series_df = series_df[['SRN', 'series']]
        scores_df = scores_df.drop(columns='series').join(series_df.set_index('SRN'), on='SRN')
    elif args.dataframe and Path(args.dataframe).is_file():
        series_df = pd.read_csv(args.dataframe)
        train_val_srns = series_df[series_df['split'] != 'test']['SRN'].to_list()
        scores_df = scores_df[scores_df['SRN'].isin(train_val_srns)].reset_index()
        scores_df['series'] = 'train'
    else:
        scores_df['series'] = 'A'

    for series in scores_df['series'].unique():
        if series != series:
            continue
        print(series)
        series_df = scores_df[scores_df['series'] == series].copy()
        series_df.loc[scores_df['rf_total'] <= 0, 'rf_total'] = 0.1

        series_df['cluster'] = ''

        for cluster in series_df['cluster'].unique():
            optimizer = PlpScoring(task, strain=strain, series=series)
            if args.benchmark:
                optimizer.optimize_weights(series_df[series_df['cluster'] != cluster])
            else:
                optimizer.optimize_weights(series_df)
            with open(f'weights_{series}.json') as f:
                weight_dict = json.load(f)
            weights = list(weight_dict.values())
            prediction_df, contacts_df = optimizer.plp_scoring_function(weights, series_df, return_contact_df=True)

            contacts_df.to_csv(f'contacts_{series}.csv', index=False)

            final_prediction_df = prediction_df.copy()
            final_prediction_df['prediction_error'] = final_prediction_df['score'] - final_prediction_df[task]
            final_prediction_df.to_csv(f'predictions_{series}.csv', index=False)
            if args.target == 'pde-10':
                task = 'pde10_h_pic50'
            print('spearman R RF rescore', spearmanr(final_prediction_df['score'], final_prediction_df[task]))
            print('spearman R Gold.PLP.Fitness',
                  spearmanr(final_prediction_df['Gold.PLP.Fitness'], final_prediction_df[task]))

    return


if __name__ == '__main__':
    main()
