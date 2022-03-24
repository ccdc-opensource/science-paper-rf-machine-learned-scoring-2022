#!/usr/bin/env python

########################################################################################################################

from rdkit.Chem import PandasTools
from pathlib import Path
import pandas as pd
from ccdc import io

########################################################################################################################


def main():
    print('Collecting results...')
    docked_solns = [str(f) for f in Path('.').glob('docking_job_*/*/rescored_ligands.sdf')]
    if not docked_solns:
        docked_solns = [str(f) for f in Path('.').glob('docking_job_*/*/best_soln_*.sdf')]

    dfs = []
    for f in docked_solns:
        mol_df = PandasTools.LoadSDF(f)
        mol_df['docking_folder'] = str(Path(f).parent)
        dfs.append(mol_df)


    df = pd.concat(dfs, ignore_index=True, sort=False)
    df['lig_id'] = df['ID'].str.split('|', expand=True)[1]
    df = df.astype({'RMSD_to_mcs': float})
    df = df[df['RMSD_to_mcs'] <= 2]
    df = df.sort_values('Gold.PLP.Fitness', ascending=False).drop_duplicates('lig_id')
    df = df.drop(columns=[c for c in df.columns if 'Gold.PLP.C' in c])
    df = df.drop(columns=[c for c in df.columns if 'Gold.PLP.p' in c])
    df = df.drop(columns=[c for c in df.columns if 'Gold.Protein' in c])
    df = df.drop(columns=[c for c in df.columns if 'Gold.Chemscore' in c])
    df = df.drop(columns=[c for c in df.columns if 'Gold.Id.' in c])
    df = df.drop(columns=[c for c in df.columns if 'Gold.PLP.S' in c])
    PandasTools.WriteSDF(df, 'best_docking_solutions.sdf', properties=list(df.columns))

    docking_folders = df['docking_folder'].to_list()
    docked_pockets = [str(f) for f in Path('.').glob('docking_job_*/*/best_soln_pocket.mol2') if str(f.parent) in docking_folders]

    rdr = io.EntryReader(docked_pockets)
    with io.EntryWriter('pockets.mol2') as w:
        for docked_pocket in rdr:
            w.write(docked_pocket)


if __name__ == '__main__':
    main()
