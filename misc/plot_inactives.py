import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path


########################################################################################################################


def main():
    print('plotting...')
    dfs = []
    for df_file in Path('.').glob('rf_score_random_model/*/rescored_rf_count.csv'):
        strucid = str(df_file.parent.stem)
        df = pd.read_csv(df_file)
        df['strucid'] = strucid
        dfs.append(df)
    dfs = pd.concat(dfs, ignore_index=True)
    dfs['scoring_function'] = 'RF-PLP'
    dfs = dfs.rename(columns={'rescore': 'prediction'})
    df2d = pd.read_csv('attentive_fp_random/ligands.csv')
    df2d['scoring_function'] = 'AttentiveFP'
    inference_df = pd.read_csv('attentive_fp_random/regression_inference_results/prediction.csv')
    df2d = df2d.join(inference_df['task_1'])
    df2d = df2d.rename(columns={'task_1': 'prediction'})
    df = pd.concat([df2d, dfs])
    df['is_parent'] = False
    df.loc[df['ligand_id'].str.len() < 8, 'is_parent'] = True
    df = df.reset_index(drop=True)
    for t, strucid_df in df.groupby(['strucid', 'scoring_function']):
        strucid, scoring_function = t
        parent_pred = strucid_df[strucid_df['is_parent'] is True]['prediction'].values[0]
        df.loc[strucid_df.index, 'parent_prediction'] = parent_pred

    df['pIC50 change'] = df['prediction'] - df['parent_prediction']
    df.loc[df['ligand_id'].str.contains('intra'), 'category'] = 'strain'
    df.loc[df['ligand_id'].str.contains('inter'), 'category'] = 'clash'
    df.loc[(df['ligand_id'].str.contains('inter') is False) & (
            df['ligand_id'].str.contains('intra') is False), 'category'] = 'contact'
    df = df[df['is_parent'] is False]

    proasis_pdb_df = pd.read_csv('../proasis_pdb.csv')
    df = df.join(proasis_pdb_df.set_index('Label'), on='strucid')
    df['IDs'] = df['IDs'].fillna('5edi').str.lower()
    df = df.rename(columns={'strucid': 'proasis_strucid', 'IDs': 'strucid'})

    df.to_csv('inactive_compounds_pic50.csv', index=False)

    sns.set_theme(style="white", font_scale=1.1)
    f, ax = plt.subplots()
    markers = ['o', '^', 'D']
    it_markers = iter(markers)
    categories = []
    strucids = df['strucid'].unique()
    for category, type_df in df.groupby('category'):
        marker = next(it_markers)
        categories.append(category)
        for strucid in strucids:
            if strucid not in type_df['strucid'].to_list():
                type_df = type_df.append({'strucid': strucid, 'pIC50 change': np.nan}, ignore_index=True)
        type_df = type_df.sort_values('strucid')
        g = sns.stripplot(data=type_df, x='strucid', y='pIC50 change', hue='scoring_function', marker=marker,
                          dodge=True, linewidth=1, ax=ax)
        g.legend_.remove()

    g.axhline(0, color='black', linewidth=0.5)
    handles, labels = g.get_legend_handles_labels()
    leg = g.legend(handles[:2], labels[:2], bbox_to_anchor=(0, -0.4), loc='lower left', ncol=len(categories))
    #
    # add legend for marker
    ax.add_artist(leg)
    h = [plt.plot([], [], color="gray", marker=marker, ls="")[0] for marker in markers]
    plt.legend(handles=h, labels=categories, bbox_to_anchor=(0, -0.3), loc='lower left', ncol=len(categories))

    f.savefig('comparison.png', dpi=300, bbox_inches='tight')
    f.clf()


if __name__ == '__main__':
    main()
