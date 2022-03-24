from pymol import cmd
from pathlib import Path


def main():
    protein_files = Path('.').glob('*.pdb')
    for protein_file in protein_files:
        strucid = protein_file.stem[:5]
        cmd.delete('all')
        cmd.remove('all')

        # Setup protein
        cmd.load(protein_file)
        cmd.select('protein', 'not ((resn GLN or resn TYR) and (resi 726 or resi 693))')
        cmd.remove('protein')
        cmd.select('bb', 'name c+o+n')
        cmd.hide('all', 'bb')
        cmd.select('protein', 'all')
        cmd.show('sticks', 'all')
        cmd.hide('(hydro)')
        cmd.hide('lines', 'all')
        cmd.select('polar', 'ele N+O')
        cmd.show('sticks', 'ele h and neighbor(ele n+o)')
        cmd.show('sticks', 'polar')
        cmd.color('green', 'ele C')

        # load ligand
        for ligand_file in Path('.').glob(f'*{strucid}*.sdf'):
            ligname = ligand_file.stem
            ligand_color = 'cyan'
            if 'native' not in ligname:
                ligand_color = 'purple'
            cmd.load(ligand_file)
            cmd.select('ligand', ligname)
            cmd.set('stick_radius', 0.1, 'ligand')
            cmd.color(ligand_color, ligname + ' and ele C')
            cmd.hide('sticks', 'ligand and ele H')

            # Finalize scene
            cmd.set('depth_cue', 0)
            cmd.set('ray_shadows', 0)
            cmd.set('valence', 0)
            cmd.zoom('all')
            cmd.set('label_size', 30)
            cmd.set('dash_width', 5)
            cmd.bg_color('white')

            ligid = ligname.split('_')[1]
            cmd.save(f'test_scene_{strucid}_{ligid}.pse')
            cmd.remove(ligname)


if __name__ == '__main__':
    main()
