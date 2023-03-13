########################################################################################################################

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from pathlib import Path
from rdkit.Chem import AllChem
from rdkit.Chem import PandasTools
from sklearn.manifold import TSNE


##########################################################################################################


class TsnePlotter():
    def __init__(self):
        self.sdf = PandasTools.LoadSDF('../pde-10/DataView_PDE10_profiling_1__export.sdf')

        fps = [AllChem.GetMorganFingerprintAsBitVect(x, radius=4, nBits=1024) for x in self.sdf['ROMol']]
        fp_df = pd.DataFrame(np.array(fps), columns=['bit_' + str(i) for i in range(len(fps[0]))])
        self.embedded_df = TSNE(n_components=2, metric='cosine').fit_transform(fp_df)
        fp_df['SRN'] = self.sdf['SRN']
        self.fp_df = fp_df
        return

    def make_tsne(self, dgl_path=Path('dgl_models_random'), output='tsne_random.png', palette='colorblind'):
        '''Collect data for manuscript Table 1: Random split'''

        dataset_df = pd.read_csv(dgl_path / 'dataset_df.csv')
        dataset_df = dataset_df[['SRN', 'split']].drop_duplicates('SRN')
        fp_df = self.fp_df.join(dataset_df.set_index('SRN'), on='SRN')
        fp_df.loc[fp_df['split'].isna(), 'split'] = 'extended dataset'
        plot_df = pd.DataFrame({'x': self.embedded_df[:, 0], 'y': self.embedded_df[:, 1], 'split': fp_df['split']})
        ax1 = sns.scatterplot(x='x', y='y', color='grey', data=plot_df[plot_df['split'] == 'extended dataset'],
                              hue='split', palette='Greys')
        # l1_handles, l1_lables = ax1.get_legend_handles_labels()
        ax2 = sns.scatterplot(x='x', y='y', hue='split', hue_order=['train', 'val', 'test'],
                              data=plot_df[plot_df['split'] != 'extended dataset'], palette=palette, edgecolor='black')
        # l2_handles, l2_lables = ax2.get_legend_handles_labels()
        # ax2.legend(l1_handles + l2_handles, l1_lables + l2_lables)

        plt.savefig(output, dpi=300, bbox_inches='tight')
        plt.clf()
        return


def main():
    TSNE = TsnePlotter()
    TSNE.make_tsne(palette='colorblind')
    TSNE.make_tsne(dgl_path=Path('dgl_models_temporal_2011'), output='tsne_temporal_2011.png', palette='colorblind')
    TSNE.make_tsne(dgl_path=Path('dgl_models_temporal_2012'), output='tsne_temporal_2012.png', palette='colorblind')
    TSNE.make_tsne(dgl_path=Path('dgl_models_temporal_2013'), output='tsne_temporal_2013.png', palette='colorblind')
    TSNE.make_tsne(dgl_path=Path('dgl_models_aminohetaryl_c1_amide'), output='tsne_aminohetaryl_c1_amide.png',
                   palette='colorblind')
    TSNE.make_tsne(dgl_path=Path('dgl_models_aryl_c1_amide_c2_hetaryl'), output='tsne_aryl_c1_amide_c2_hetaryl.png',
                   palette='colorblind')
    TSNE.make_tsne(dgl_path=Path('dgl_models_c1_hetaryl_alkyl_c2_hetaryl'),
                   output='tsne_c1_hetaryl_alkyl_c2_hetaryl.png', palette='colorblind')


if __name__ == '__main__':
    main()
