#!/usr/bin/env python

'''
Concatenate all rf_count dfs after assigning to each docking_dir
'''

########################################################################################################################

import sys
import pandas as pd
from pathlib import Path
import argparse
from ccdc import io
from rdkit.Chem import PandasTools

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

    parser.add_argument(
        '-r',
        '--rescore',
        help='Rescore ligand poses.',
        action='store_true',
        default=False
    )

    return parser.parse_args()


def get_best_docking_solutions(df, rescoring=False):

    if rescoring:
        filtered_df = df.sort_values(
            by=['rescore', 'rel_mcs_size', 'tanimoto_similiarity_to_native_ligand'],
            ascending=[False, False, False]).drop_duplicates('ligand_id')

    else:
        df['is_decoy'] = False
        if 'SRN' not in df.columns:
            sys.exit('Ligands do not have SRN, try rescoring mode.')
        filtered_df = df.sort_values(by=['Gold.PLP.Fitness', 'rel_mcs_size', 'tanimoto_similiarity_to_native_ligand'],
                                     ascending=[False, False, False])
        filtered_df['SRN'] = filtered_df['SRN'].str.split('-', expand=True)[0]
        filtered_df = filtered_df.drop_duplicates('SRN')

    filtered_df = filtered_df.drop(columns=[c for c in filtered_df.columns if 'Gold.Protein' in c or 'Gold.Chemscore.Hbonds' in c])
    filtered_df = filtered_df[(filtered_df['RMSD_to_mcs'] <= 1.5) & (filtered_df['rel_mcs_size'] >= 0.5)]

    return filtered_df, df


def concatenate_rf_counts(args):

    if args.rescore:
        dfs = list(Path('.').glob('*docking_job*/rescored_rf_count.csv'))
        scoring_function = 'rescore'
        output = '_rescored'
    else:
        dfs = list(Path('.').glob('*docking_job*/docked_rf_count_df.csv'))
        scoring_function = 'Gold.PLP.Fitness'
        output = ''
    if 'decoy' in str(dfs[0]):

        dfs = [pd.read_csv(df) for df in dfs if Path(df).is_file()]
        edited_dfs = []
        for df in dfs:
            if 'TanimotoCombo_to_scaffold' in df.columns:
                if 'TanimotoCombo_to_reference' in df.columns:
                    df = df.drop('TanimotoCombo_to_scaffold', axis=1)
                else:
                    df = df.rename(columns={'TanimotoCombo_to_scaffold': 'TanimotoCombo_to_reference'})
            edited_dfs.append(df)
        dfs = edited_dfs
        df = pd.concat(dfs, ignore_index=True)
        df['is_decoy'] = True
        df.to_csv('docked_rf_count_df.csv', index=False)
        filtered_df = df.sort_values(by=['rel_mcs_size', scoring_function], ascending=[False, True])
        active_srn = pd.read_csv('../docked_rf_count_df.csv')['SRN']
        filtered_df = filtered_df[filtered_df['SRN'].isin(active_srn)]
        filtered_df = filtered_df[filtered_df['TanimotoCombo_to_reference'] < 0.4]
        output_file = f'docked{output}_rf_count_best_decoys.csv'
        filtered_df.to_csv(output_file, index=False)

    else:
        dfs = [pd.read_csv(df) for df in dfs if Path(df).is_file()]
        df = pd.concat(dfs, ignore_index=True)
        filtered_df, df = get_best_docking_solutions(df, rescoring=args.rescore)
        df = df.drop(columns=[c for c in filtered_df.columns if 'Gold.Protein' in c or 'Gold.Chemscore.Hbonds' in c])
        df.to_csv(f'docked{output}_rf_count_df.csv', index=False)
        output_file = f'docked{output}_rf_count_best_docked.csv'
        filtered_df.to_csv(output_file, index=False)


def concatenate_docked_strucs():
    ligand_files = Path('.').glob('*docking_job*/*_*/gold_soln_docking_input_*.sdf')

    if 'decoy' in ligand_files[0]:
        ligand_dicts = []
        for ligand_file in ligand_files:
            with io.EntryReader(ligand_file) as rdr:
                e = rdr[0]
            attributes = e.attributes
            if 'TanimotoCombo_to_reference' in attributes:
                if float(attributes['TanimotoCombo_to_reference']) < 0.4:
                    attributes['ligand_file'] = ligand_file
                    ligand_dicts.append(attributes)
        df = pd.DataFrame.from_dict(ligand_dicts)
        df.to_csv('all_good_decoys.csv', index=False)
    else:
        ligand_dicts = []
        for ligand_file in ligand_files:
            with io.EntryReader(ligand_file) as rdr:
                e = rdr[0]
            attributes = e.attributes
            if float(attributes['RMSD_to_mcs']) <= 1.5:
                attributes['ligand_file'] = ligand_file
                ligand_dicts.append(attributes)
        df = pd.DataFrame.from_dict(ligand_dicts)
        df.to_csv('all_good_dockings.csv', index=False)
    return


def main():
    args = parse_args()
    concatenate_rf_counts(args)


if __name__ == '__main__':
    main()
