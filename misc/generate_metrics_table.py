#!/usr/bin/env python

import pandas as pd
from rdkit.Chem import PandasTools
from scipy.stats import spearmanr, pearsonr, bootstrap
from sklearn.metrics import mean_squared_error
from pathlib import Path
import numpy as np
# from sklearn.preprocessing import MinMaxScaler
from sklearn.linear_model import LinearRegression


def add_attentive_fp(attentive_fp_df, predictions_df, dgl_path, extension=''):
    afp_column_name = 'AttentiveFP' + extension
    hybrid_column_name = '2D3D' + extension
    smiles_df = pd.read_csv(dgl_path / 'attentive_fp/2d_dataset_df.csv')
    smiles_df = smiles_df[['SMILES', 'SRN']]
    attentive_fp_df = attentive_fp_df.join(smiles_df.set_index('SMILES'), on='smiles')
    attentive_fp_df = attentive_fp_df[['SRN', 'prediction']]
    attentive_fp_df = attentive_fp_df.rename(columns={'prediction': afp_column_name})
    attentive_fp_df['SRN'] = attentive_fp_df['SRN'].str.split('-', expand=True)[0]
    predictions_df = predictions_df.join(attentive_fp_df.set_index('SRN'), on='SRN')
    predictions_df[hybrid_column_name] = (predictions_df['RF-PLP'] + predictions_df[afp_column_name]) / 2
    return predictions_df


def return_pearsonr(x, y):
    return pearsonr(x, y)[0]


def return_rmse(x, y):
    return mean_squared_error(x, y, squared=False)


def return_r2(x, y):
    return pearsonr(x, y)[0]**2


def return_spearmanr(x, y):
    return spearmanr(x, y)[0]


def return_ci(x, y, fn):
    ci = bootstrap((x, y), fn, vectorized=False, paired=True).confidence_interval
    ci = (ci[1]-ci[0]) / 2
    return ci


class TableMaker(object):
    def __init__(self):
        vanilla_df = PandasTools.LoadSDF('vanilla_gold/vanilla_gold.sdf')
        vanilla_df = vanilla_df.drop(columns='ROMol')
        vanilla_df['SRN'] = vanilla_df['ID'].str.split('|', expand=True)[1].str.split('-', expand=True)[0]
        vanilla_df = vanilla_df.astype({'Gold.PLP.Fitness': float})
        vanilla_df = vanilla_df.sort_values('Gold.PLP.Fitness', ascending=False).drop_duplicates('SRN')
        vanilla_df = vanilla_df[vanilla_df['Gold.PLP.Fitness'] > 0]
        self.vanilla_df = vanilla_df

    def make_table(self, dgl_path=Path('dgl_models_random'), rf_score_path=Path('docking_template_date/first_iteration_random_split'), output=Path('table1_with_avg_scaled.csv')):
        '''Collect data for manuscript Table 1: Random split'''

        dataset_df = pd.read_csv(dgl_path / 'dataset_df.csv')
        dataset_df['SRN'] = dataset_df['SRN'].str.split('-', expand=True)[0]
        test_srns = dataset_df[dataset_df['split'] == 'test']['SRN'].unique()
        train_srns = dataset_df[dataset_df['split'].isin(['val', 'train'])]['SRN'].unique()

        predictions_df = dataset_df[['SRN', 'split', 'pde10_h_pic50']]
        predictions_df = predictions_df[predictions_df['SRN'].isin(test_srns)]

        # training set pIC50 average as null modell

        train_set_mean = dataset_df[dataset_df['split'].isin(['val', 'train'])]['pde10_h_pic50'].mean()
        predictions_df['train_set_mean'] = train_set_mean

        # Vanilla Gold
        vanilla_test_df = self.vanilla_df[self.vanilla_df['SRN'].isin(test_srns)]
        vanilla_test_df = vanilla_test_df[['SRN', 'Gold.PLP.Fitness']]
        vanilla_test_df = vanilla_test_df.rename(columns={'Gold.PLP.Fitness': 'vanilla_Gold.PLP.Fitness'})

        # Fit docking score to pIC50
        vanilla_train_df = self.vanilla_df[self.vanilla_df['SRN'].isin(train_srns)]
        vanilla_train_df = vanilla_train_df.join(dataset_df[dataset_df['SRN'].isin(train_srns)][['SRN', 'pde10_h_pic50']].set_index('SRN'), on='SRN')
        linear_fit = LinearRegression().fit(vanilla_train_df['Gold.PLP.Fitness'].to_numpy().reshape(-1, 1), vanilla_train_df['pde10_h_pic50'].to_numpy().reshape(-1, 1))
        vanilla_test_df['vanilla_Gold.PLP.Fitness_scaled'] = linear_fit.predict(
            vanilla_test_df['vanilla_Gold.PLP.Fitness'].to_numpy().reshape(-1, 1))

        predictions_df = predictions_df.join(vanilla_test_df.set_index('SRN'), on='SRN')

        # Template Docking
        template_df = pd.read_csv(rf_score_path / 'rescored_rf_count.csv')
        template_df = template_df[['SRN', 'Gold.PLP.Fitness', 'rescore']]
        template_df = template_df[template_df['SRN'].isin(test_srns)]
        template_df = template_df.rename(columns={'Gold.PLP.Fitness': 'template_Gold.PLP.Fitness',
                                                  'rescore': 'RF-PLP'})

        predictions_df = predictions_df.join(template_df.set_index('SRN'), on='SRN')

        attentive_fp_df = pd.read_csv(dgl_path / 'attentive_fp/predictions.csv')
        predictions_df = add_attentive_fp(attentive_fp_df, predictions_df, dgl_path)

        # ACNN
        acnn_df = pd.read_csv(dgl_path / 'ACNN_predictions.csv')
        acnn_df = acnn_df[['SRN', 'pred_pic50']]
        acnn_df = acnn_df.rename(columns={'pred_pic50': 'ACNN'})
        acnn_df['SRN'] = acnn_df['SRN'].str.split('-', expand=True)[0]
        predictions_df = predictions_df.join(acnn_df.set_index('SRN'), on='SRN')

        # PotentialNet
        pnet_df = pd.read_csv(dgl_path / 'PotentialNet_predictions.csv')
        pnet_df = pnet_df[['SRN', 'pred_pic50']]
        pnet_df = pnet_df.rename(columns={'pred_pic50': 'PotentialNet'})
        pnet_df['SRN'] = pnet_df['SRN'].str.split('-', expand=True)[0]
        predictions_df = predictions_df.join(pnet_df.set_index('SRN'), on='SRN')

        extended_df_path = (dgl_path / 'attentive_fp_extended_training_set/predictions.csv')
        if extended_df_path.is_file():
            attentive_fp_df = pd.read_csv(extended_df_path)
            predictions_df = add_attentive_fp(attentive_fp_df, predictions_df, dgl_path, extension='_ed')

        table_dict = {'rmse': [], 'rmse_ci': [], 'spearman': [], 'spearman_ci': [], 'r2': [], 'r2_ci': [], 'model': []}

        for model in predictions_df.columns:
            if model in ['SRN', 'split', 'pde10_h_pic50']:
                continue
            if 'scaled' in model and 'Gold' in model:
                continue
            label = predictions_df['pde10_h_pic50'].to_numpy()
            predictions = predictions_df[model].to_numpy()

            # spearmanr and pearsonr not defined for constant arrays
            if model == 'train_set_mean':
                table_dict['spearman'].append(np.nan)
                table_dict['spearman_ci'].append(np.nan)
                table_dict['r2'].append(np.nan)
                table_dict['r2_ci'].append(np.nan)
            else:
                table_dict['spearman'].append(round(return_spearmanr(label, predictions), 2))
                table_dict['spearman_ci'].append(return_ci(label, predictions, return_spearmanr))
                table_dict['r2'].append(round(return_r2(label, predictions), 2))
                table_dict['r2_ci'].append(return_ci(label, predictions, return_r2))

            if 'Gold' in model:
                _predictions = predictions_df[model + '_scaled'].to_numpy()
            else:
                _predictions = predictions
            table_dict['rmse'].append(round(return_rmse(label, _predictions), 2))
            table_dict['rmse_ci'].append(return_ci(label, _predictions, return_rmse))

            table_dict['model'].append(model)
        table_df = pd.DataFrame(table_dict)
        table_df = table_df[['model', 'spearman', 'spearman_ci', 'rmse', 'rmse_ci', 'r2', 'r2_ci']]
        table_df.to_csv(output, index=False)


def main():
    table_maker = TableMaker()
    # Table 1
    table_maker.make_table()

    table_maker.make_table(dgl_path=Path('dgl_models_temporal_2011'),
                           rf_score_path=Path('docking_template_date/first_iteration_temporal_2011_split'),
                           output=Path('table1_2011_scaled.csv'))

    table_maker.make_table(dgl_path=Path('dgl_models_temporal_2012'),
                           rf_score_path=Path('docking_template_date/first_iteration_temporal_2012_split'),
                           output=Path('table1_2012_scaled.csv'))

    table_maker.make_table(dgl_path=Path('dgl_models_temporal_2013'),
                           rf_score_path=Path('docking_template_date/first_iteration_temporal_2013_split'),
                           output=Path('table1_2013_scaled.csv'))

    # Table 2
    table_maker.make_table(dgl_path=Path('dgl_models_aminohetaryl_c1_amide'),
               rf_score_path=Path('docking_template_date/first_iteration_aminohetaryl_c1_amide'),
               output=Path('table3_aminohetaryl_c1_amide_with_avg_scaled.csv'))


    table_maker.make_table(dgl_path=Path('dgl_models_c1_hetaryl_alkyl_c2_hetaryl'),
               rf_score_path=Path('docking_template_date/first_iteration_c1_hetaryl_alkyl_c2_hetaryl'),
               output=Path('table4_c1_hetaryl_alkyl_c2_hetaryl_with_avg_scaled.csv'))

    table_maker.make_table(dgl_path=Path('dgl_models_aryl_c1_amide_c2_hetaryl'),
               rf_score_path=Path('docking_template_date/first_iteration_aryl_c1_amide_c2_hetaryl'),
               output=Path('table5_aryl_c1_amide_c2_hetaryl_with_avg_scaled.csv'))

    return


if __name__ == "__main__":
    main()
