
########################################################################################################################

from ccdc import io
import pandas as pd
from pathlib import Path
import shutil
########################################################################################################################


def blind(mol, pdb_complex_id, extension):
    mol.identifier = pdb_complex_id
    with io.MoleculeWriter(Path(pdb_complex_id) / (pdb_complex_id + extension)) as w:
        w.write(mol)
    return


def main():
    # proasis_pdb_df = pd.read_csv('../proasis_pdb.csv')
    # for complex_folder in Path('../pde-10_pdb_bind_format').glob('*_*'):
    #     if complex_folder.is_dir():
    #         proasis_complex_id = complex_folder.stem
    #         proasis_id = proasis_complex_id.split('_')[0]
    #         pdb_id = proasis_pdb_df[proasis_pdb_df['Label'] == proasis_id]['IDs'].str.lower().values[0]
    #         compound_id = complex_folder.stem.split('_')[1]
    #         pdb_complex_id = pdb_id + '_' + compound_id
    #         blinded_complex_path = Path(pdb_complex_id)
    #         if not blinded_complex_path.is_dir():
    #             blinded_complex_path.mkdir()
    #
    #         for extension in ['_protein.pdb', '_pocket.pdb', '_ligand.sdf', '_ligand.mol2']:
    #             mol = io.MoleculeReader(str(complex_folder / (proasis_complex_id + extension)))[0]
    #             blind(mol, pdb_complex_id, extension)

    table_s1_df = pd.read_csv('../table_s1.csv')
    table_s1_df = table_s1_df[table_s1_df['temporal_2013_split'] == 'test']
    for index, docking_folder in table_s1_df['docking_folder'].iteritems():
        shutil.copytree(docking_folder, '../aqemia_2013_test_set/' + docking_folder)
    return


if __name__ == '__main__':
    main()

