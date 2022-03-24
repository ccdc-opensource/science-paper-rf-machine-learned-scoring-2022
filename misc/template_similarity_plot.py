
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

docking_df = pd.read_csv('docking_template_date/docked_rf_count_best_docked.csv')
series_df = pd.read_csv('series_assignment.csv')
series_df = series_df[['SRN', 'series']]
series_df['SRN'] = series_df['SRN'].str.split('-', expand=True)[0]
docking_df = docking_df.drop(columns='series').join(series_df.set_index('SRN'), on='SRN')
docking_df = docking_df.rename(columns={'tanimoto_similiarity_to_native_ligand': 'Tanimoto Similarity', 'rel_mcs_size': 'Relative MCS size'})
ax = sns.scatterplot(x='Tanimoto Similarity', y='Relative MCS size', color='grey', data=docking_df,
                              hue='series', palette='colorblind')
sns.move_legend(ax, "lower center", bbox_to_anchor=(.5, -0.35), ncol=1, title=None, frameon=False)
plt.savefig('manuscript_data/template_sim_vs_mcs_size.png', dpi=300, bbox_inches='tight')
plt.clf()
print('nice')