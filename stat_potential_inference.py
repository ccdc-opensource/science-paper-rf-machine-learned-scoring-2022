#!/usr/bin/env python

########################################################################################################################

import argparse
import json
import numpy as np
import pandas as pd
from ccdc import io, search
from ccdc_roche_scoring import descriptor_dfs, stat_potential, scoring_parameters, docking
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

    parser.add_argument(
        '-s',
        '--series',
        help='Series name.',
        default=''
    )

    parser.add_argument(
        '-l',
        '--ligand',
        help='Ligand file.',
        default=''
    )

    parser.add_argument(
        '-p',
        '--protein',
        help='Protein file.',
        default=''
    )

    parser.add_argument('-g', '--gold', action='store_true',
                        help='Input is Gold solution file')

    parser.add_argument('-df', '--dataframe', default=False,
                        help='Input is CSV file')

    parser.add_argument(
        '-d',
        '--docking_path',
        help='docking_path.',
        default=''
    )

    return parser.parse_args()


class StatPotInferer():
    def __init__(self, target, series, ligand_file, gold=True, protein=False, mol_contact_df=pd.DataFrame(),
                 rf_count_df=pd.DataFrame(), scoring_parameter_file=Path('.'), task='pde10_h_pic50'):
        self.target = target
        self.target_home = f'/pstore/db/cadd/moe/moe_projects/{target}_bindingsites'
        if not Path(self.target_home).is_dir():
            self.target_home = '.'
        self.protein = protein
        self.task = task
        self.series = series
        if scoring_parameter_file == Path('.'):
            scoring_parameter_file = Path('').resolve().parent / f'second_iteration/weights_{series}.json'
            if not scoring_parameter_file.is_file():
                scoring_parameter_file = Path('').resolve().parent / f'first_iteration/weights_{series}.json'
            if not scoring_parameter_file.is_file():
                scoring_parameter_path = Path(scoring_parameters.__path__._path[0])
                scoring_parameter_file = scoring_parameter_path / self.target / f'weights_{series}.json'
        self.ligand_file = ligand_file

        if mol_contact_df.shape[0] == 0:
            self.gold = gold
            self._return_dfs()
        else:
            self.mol_contact_df = mol_contact_df
            self.rf_count_df = rf_count_df

        if scoring_parameter_file.is_file() and self.mol_contact_df.shape[0] > 0:
            with open(scoring_parameter_file) as f:
                weight_dict = json.load(f)
            self.weights = list(weight_dict.values())
            self.score = self._return_stat_pot_score()
        else:
            self.score = False

    def _return_dfs(self):
        if self.gold:
            rf_count_df, mol_contact_df = descriptor_dfs.return_contact_df_from_docking_file(
                self.ligand_file, target_home=self.target_home,
                strucid=self.ligand_file.parent.resolve().name.split('_')[0])

        else:
            rdkit_protein = Chem.MolFromPDBFile(self.protein, removeHs=False)
            self.protein = './sanitized.pdb'
            print(Path('.').resolve())
            Chem.MolToPDBFile(rdkit_protein, self.protein)
            docking.gold_rescoring(str(self.ligand_file), self.protein)

            rf_count_df, mol_contact_df = pd.DataFrame(), pd.DataFrame()
            if Path('rescored_ligand.sdf').is_file():
                rf_count_df, mol_contact_df = descriptor_dfs.return_contact_df_from_docking_file(
                    'rescored_ligand.sdf', pdb_file=self.protein)
        self.rf_count_df = rf_count_df
        self.mol_contact_df = mol_contact_df

    def _return_stat_pot_score(self):
        self.mol_contact_df = stat_potential.prepare_contact_df(self.mol_contact_df, self.rf_count_df)
        self.mol_contact_df = self.mol_contact_df.drop(
            columns=[c for c in self.mol_contact_df.columns if 'Gold.PLP.part' in c or c in ['Gold.Version']])
        self.mol_contact_df['cluster'] = ''
        if self.task not in self.mol_contact_df.columns:
            self.mol_contact_df[self.task] = np.nan

        optimizer = stat_potential.PlpScoring(
            self.task, strain=['Gold.PLP.ligand.torsion', 'Gold.PLP.ligand.clash', 'Gold.PLP.Chemscore.Protein.Energy'],
            series=self.series, gold=False)
        prediction_df, self.mol_contact_df = optimizer.plp_scoring_function(self.weights, self.mol_contact_df,
                                                                            return_contact_df=True)
        self.rf_count_df['rescore'] = prediction_df['score'].values[0]
        self.mol_contact_df['rescore'] = prediction_df['score'].values[0]
        if self.task in self.rf_count_df.columns:
            self.rf_count_df['error'] = self.rf_count_df['rescore'] - self.rf_count_df[self.task]
            self.mol_contact_df['error'] = self.mol_contact_df['rescore'] - self.mol_contact_df[self.task]
        self.rf_count_df['ligand_file'] = Path(self.ligand_file).resolve()
        self.mol_contact_df['ligand_file'] = Path(self.ligand_file).resolve()
        print(round((self.mol_contact_df['rescore'].values[0]), 2))
        return self.rf_count_df, self.mol_contact_df


def main():
    args = parse_args()
    print(args.target)
    # scorer = StatPotInferer(args.target, args.series, args.ligand)
    rf_count_df_list = []
    mol_contact_df_list = []

    # score from ligands
    if not args.dataframe:
        if Path(args.ligand).is_file():
            ligand_files = [Path(args.ligand)]
        else:
            ligand_files = Path(args.docking_path).glob('*_*/best_soln*.sdf')
        with io.EntryWriter('rescored_ligands.sdf') as output_ligand_writer:
            for ligand_file in ligand_files:
                ligand_file = Path(ligand_file)
                ligands = io.EntryReader(str(ligand_file))
                for ligand in ligands:
                    temp_ligand_file = ligand_file.parent / Path('temp_ligand.sdf')
                    with io.EntryWriter(temp_ligand_file) as w:
                        w.write(ligand)
                    ligand_series = ''
                    scoring_parameter_file = Path('../weights_train.json')
                    if 'series' in ligand.attributes.keys():
                        ligand_series = ligand.attributes['series']
                        scoring_parameter_file = Path(
                            __file__).parent / f'scoring_parameters/{args.target}/weights_{ligand_series}.json'

                    else:
                        with open(
                                Path(__file__).parent / Path(f'series_definitions/{args.target}.json')) as series_defs:
                            series_dict = json.load(series_defs)
                        for series in series_dict:
                            smarts = search.SMARTSSubstructure(series_dict[series])
                            searcher = search.SubstructureSearch()
                            searcher.add_substructure(smarts)
                            hits = searcher.search(ligand.molecule)
                            if hits:
                                ligand_series = series
                                scoring_parameter_file = Path(
                                    __file__).parent / f'scoring_parameters/{args.target}/weights_{ligand_series}.json'
                                if not scoring_parameter_file.is_file():
                                    scoring_parameter_file = Path(
                                        __file__).parent / f'scoring_parameters/{args.target}/weights_train.json'
                                break
                    if ligand_series:
                        print('scoring...')
                        scorer = StatPotInferer(args.target, ligand_series, temp_ligand_file, gold=args.gold,
                                                protein=args.protein,
                                                scoring_parameter_file=scoring_parameter_file)
                        if scorer.score:
                            rf_count_df, mol_contact_df = scorer.score
                            rf_count_df['series'] = ligand_series
                            mol_contact_df['series'] = ligand_series
                            rf_count_df_list.append(rf_count_df)
                            mol_contact_df_list.append(mol_contact_df)
                            ligand.attributes['RF Score'] = round(rf_count_df['rescore'].values[0], 1)
                            output_ligand_writer.write(ligand)
                    else:
                        print('No scoring function for this series.')

    # Score from Contact DF
    elif args.dataframe:
        mol_contact_dfs = pd.read_csv('best_soln_contacts_processed.csv')
        mol_contact_dfs['SRN'] = mol_contact_dfs['SRN'].str.split('-', expand=True)[0]

        rf_count_dfs = pd.read_csv('../docked_rf_count_best_docked.csv')
        split_df = pd.read_csv(args.dataframe)
        split_df['SRN'] = split_df['SRN'].str.split('-', expand=True)[0]
        test_srns = split_df[split_df['split'] == 'test']['SRN'].to_list()

        mol_contact_dfs = mol_contact_dfs[mol_contact_dfs['SRN'].isin(test_srns)]
        rf_count_dfs = rf_count_dfs[rf_count_dfs['SRN'].isin(test_srns)]
        ligand_series = 'train'
        for ligand_file, mol_contact_df in mol_contact_dfs.groupby('ligand_file'):
            rf_count_df = rf_count_dfs[rf_count_dfs['ligand_file'] == ligand_file]
            scorer = StatPotInferer(args.target, ligand_series, ligand_file, mol_contact_df=mol_contact_df,
                                    rf_count_df=rf_count_df,
                                    scoring_parameter_file=Path('weights_train.json'))
            if scorer.score:
                rf_count_df, mol_contact_df = scorer.score
                rf_count_df['series'] = ligand_series
                mol_contact_df['series'] = ligand_series
                rf_count_df_list.append(rf_count_df)
                mol_contact_df_list.append(mol_contact_df)

    if rf_count_df_list:
        rf_count_df = pd.concat(rf_count_df_list, ignore_index=True)
        rf_count_df = rf_count_df.drop(columns=[c for c in rf_count_df.columns if 'Gold.Protein' in str(c)])
        rf_count_df = rf_count_df.drop(columns=[c for c in rf_count_df.columns if 'Chemscore' in str(c)])
        rf_count_df.to_csv('rescored_rf_count.csv', index=False)
        mol_contact_df = pd.concat(mol_contact_df_list, ignore_index=True)
        mol_contact_df = mol_contact_df.drop(columns=[c for c in mol_contact_df.columns if 'Gold.Protein' in str(c)])
        mol_contact_df = mol_contact_df.drop(columns=[c for c in mol_contact_df.columns if 'Chemscore' in str(c)])
        mol_contact_df.to_csv('rescored_rf_contact.csv', index=False)
    else:
        pd.DataFrame({'rescore': [np.nan]}).to_csv('rescored_rf_count.csv', index=False)
        with io.EntryWriter('rescored_ligands.sdf') as w, io.EntryReader(args.ligand) as rdr:
            for e in rdr:
                e.attributes['RF Score'] = ''
                w.write(e)


if __name__ == '__main__':
    main()
