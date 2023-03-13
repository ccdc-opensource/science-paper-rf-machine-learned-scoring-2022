#!/usr/bin/env python

########################################################################################################################

import argparse
import itertools
import os
import pandas as pd
import subprocess as sp
import sys
from ccdc import io, protein, descriptors, entry, search
from ccdc.docking import Docker
from ccdc_roche.python.los_descriptors import _cut_out_binding_site_by_distance
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
        '--los_home',
        help='Path to LoS home folder.',
        default=''
    )

    parser.add_argument(
        '--input_ligands',
        help='SDF file with aligned ligands',
        default=False
    )

    parser.add_argument(
        '--docking_dir',
        help='Directory with docked ligands',
        default=False
    )

    parser.add_argument(
        '-t',
        '--target',
        help='Target name.',
        default=None
    )

    parser.add_argument(
        '-fr',
        '--flexible_residues',
        nargs='+',
        help='Target name.',
        default=[]
    )

    return parser.parse_args()


class MCS(object):
    '''Calculate Maximum Common Substructure between ligand and template.'''

    def __init__(self, ligand, template):
        '''

        :param ligand:
        :param template:
        :param strict_mcs Set to True if a stricter MCS should be calculated. :
        '''
        self.ligand = ligand
        self.template = template
        self.mcs_scaffold, self.mcs_atoms, self.mcs_bonds = self.return_mcs_scaffold()

    # def return_strict_mcs(self):
    #     strict_mcs_searcher = descriptors.MolecularDescriptors.MaximumCommonSubstructure()
    #     strict_mcs_searcher.settings.check_bond_count = True
    #     strict_mcs_searcher.settings.check_hydrogen_count = True
    #     strict_scaffold = self.template.copy()
    #     _ligand = self.ligand.copy()
    #     strict_scaffold.standardise_aromatic_bonds()
    #     _ligand.standardise_aromatic_bonds()
    #     strict_scaffold.standardise_delocalised_bonds()
    #     _ligand.standardise_delocalised_bonds()
    #     strict_mcs_atoms = strict_mcs_searcher.search(strict_scaffold, _ligand, search_step_limit=1000000)[0]
    #     strict_mcs_atoms = [a[0] for a in strict_mcs_atoms]
    #     remove_atoms = [a for a in strict_scaffold.atoms if a not in strict_mcs_atoms]
    #
    #     # remove partial ring matches
    #     strict_mcs_searcher = descriptors.MolecularDescriptors.MaximumCommonSubstructure()
    #     strict_mcs_searcher.settings.check_bond_count = False
    #     strict_mcs_searcher.settings.check_hydrogen_count = True
    #     while remove_atoms:
    #         strict_scaffold.remove_atoms(remove_atoms)
    #         max_component_size = 0
    #         if strict_scaffold.atoms:
    #             for c in strict_scaffold.components:
    #                 c_size = len(c.atoms)
    #                 if c_size > max_component_size:
    #                     max_component = c
    #                     max_component_size = c_size
    #             strict_scaffold = max_component
    #         else:
    #             break
    #         mcs_atoms, mcs_bonds = strict_mcs_searcher.search(strict_scaffold, _ligand)
    #         remove_atoms = self._mcs_partial_ring_match_from_bond(mcs_bonds, strict=True)
    #         if not remove_atoms:
    #             strict_scaffold.add_hydrogens(mode='missing')
    #             remove_atoms = remove_atoms + _compare_mcs_stereo_chemistry(strict_scaffold, _ligand, mcs_atoms)
    #
    #         mcs_atoms = [a[0] for a in mcs_atoms]
    #         remove_atoms = remove_atoms + [a for a in strict_scaffold.atoms if a not in mcs_atoms]
    #         remove_atoms = list(set(remove_atoms))
    #         if not [a for a in remove_atoms if a.atomic_symbol != 'H']:
    #             strict_scaffold.remove_atoms(remove_atoms)
    #             break
    #
    #     mcs_atoms, mcs_bonds = strict_mcs_searcher.search(strict_scaffold, _ligand)
    #     remove_atoms = SubstitutentComparer().compare(strict_scaffold, _ligand, mcs_atoms)
    #     # remove_atoms = _compare_aromatic_rings(strict_scaffold, _ligand, mcs_atoms)
    #     strict_scaffold.remove_atoms(remove_atoms)
    #     max_component_size = 0
    #     if strict_scaffold.atoms:
    #         for c in strict_scaffold.components:
    #             c_size = len(c.atoms)
    #             if c_size > max_component_size:
    #                 max_component = c
    #                 max_component_size = c_size
    #         strict_scaffold = max_component
    #
    #     if len(strict_scaffold.heavy_atoms) < 7:
    #         strict_scaffold = False
    #
    #     return strict_scaffold

    def return_mcs_scaffold(self, partial_ring_matches_allowed=True, ignore_hydrogens=True):
        mcs_searcher = descriptors.MolecularDescriptors.MaximumCommonSubstructure()
        mcs_searcher.settings.ignore_hydrogens = ignore_hydrogens
        scaffold = self.template.copy()

        mcs_searcher.settings.check_bond_type = False
        mcs_atoms, mcs_bonds = mcs_searcher.search(scaffold, self.ligand, search_step_limit=1000000)

        if partial_ring_matches_allowed:
            mcs_bonds = [mcs_pair for mcs_pair in mcs_bonds if mcs_pair[0].is_cyclic == mcs_pair[1].is_cyclic]
            mcs_scaffold_atoms_from_bonds = list(itertools.chain.from_iterable([b[0].atoms for b in mcs_bonds]))
            mcs_scaffold_atoms_from_bonds = [a.label for a in mcs_scaffold_atoms_from_bonds]
            mcs_atoms = [mcs_pair for mcs_pair in mcs_atoms if mcs_pair[0].label in mcs_scaffold_atoms_from_bonds]
        else:
            mcs_atoms = [mcs_pair for mcs_pair in mcs_atoms if mcs_pair[0].is_cyclic == mcs_pair[1].is_cyclic]

        mcs_template_atoms = [a[0] for a in mcs_atoms]
        scaffold_atoms = scaffold.atoms
        if ignore_hydrogens:
            remove_atoms = [a for a in scaffold_atoms if a not in mcs_template_atoms and a.atomic_symbol != 'H']
        else:
            remove_atoms = [a for a in scaffold_atoms if a not in mcs_template_atoms]

        if partial_ring_matches_allowed:
            remove_atoms = remove_atoms + self._mcs_partial_ring_match_from_bond(mcs_bonds)
        else:
            remove_atoms = remove_atoms + self._mcs_partial_ring_match(mcs_atoms)
            remove_atoms = remove_atoms + SubstitutentComparer().compare(self.ligand, mcs_atoms)

        remove_atoms = remove_atoms + _compare_mcs_stereo_chemistry(scaffold, self.ligand, mcs_atoms)
        remove_atoms = list(set(remove_atoms))

        # remove partial ring matches
        while remove_atoms:
            scaffold.remove_atoms(remove_atoms)
            max_component_size = 0
            if scaffold.atoms:
                for c in scaffold.components:
                    c_size = len(c.atoms)
                    if c_size > max_component_size:
                        max_component = c
                        max_component_size = c_size
                scaffold = max_component
            else:
                break
            mcs_atoms, mcs_bonds = mcs_searcher.search(scaffold, self.ligand, search_step_limit=1000000)

            # mcs_atoms_ = mcs_searcher.search(scaffold, self.ligand, search_step_limit=1000000)[0]
            # remove_atoms = [mcs_pair[0] for mcs_pair in mcs_atoms if mcs_pair not in mcs_atoms_]
            remove_atoms = []

            if partial_ring_matches_allowed:
                remove_atoms = remove_atoms + self._mcs_partial_ring_match_from_bond(mcs_bonds)
            else:
                remove_atoms = self._mcs_partial_ring_match(mcs_atoms) + remove_atoms
                remove_atoms = remove_atoms + SubstitutentComparer().compare(self.ligand, mcs_atoms)
            remove_atoms = remove_atoms + _compare_mcs_stereo_chemistry(scaffold, self.ligand, mcs_atoms)
            remove_atoms = list(set(remove_atoms))
        return scaffold, mcs_atoms, mcs_bonds

    def _mcs_partial_ring_match(self, mcs_atoms):
        '''
        :param mcs_atoms:
        :return: List of atoms that have different ring attributes.
        '''
        remove_atoms = []
        for mcs_pair in mcs_atoms:
            atom0 = mcs_pair[0]
            atom1 = mcs_pair[1]
            if atom0.is_cyclic != atom1.is_cyclic:
                remove_atoms.append(atom0)
            elif atom0.is_spiro != atom1.is_spiro:
                remove_atoms.append(atom0)
            elif atom0.is_cyclic:
                if len(atom0.rings) != len(atom1.rings):
                    remove_atoms.append(atom0)

        return remove_atoms

    def _mcs_partial_ring_match_from_bond(self, mcs_bonds, strict=False):
        remove_atoms = []
        safe_atoms = []
        for mcs_pair in mcs_bonds:
            bond0 = mcs_pair[0]
            bond1 = mcs_pair[1]
            if bond0.is_cyclic != bond1.is_cyclic:
                if strict:
                    atom0 = mcs_pair[0].atoms[0]
                    if atom_has_open_valency(atom0):
                        remove_atoms.append(atom0)
                    atom1 = mcs_pair[0].atoms[1]
                    if atom_has_open_valency(atom1):
                        remove_atoms.append(atom1)
                else:
                    remove_atoms.append(bond0.atoms[0])
                    remove_atoms.append(bond0.atoms[1])
                    if [a for a in bond0.atoms if a.atomic_symbol == 'H']:
                        safe_atoms.append(bond0.atoms[0].label)
                        safe_atoms.append(bond0.atoms[1].label)
        remove_atoms = [a for a in remove_atoms if a.label not in safe_atoms]

        for mcs_pair in mcs_bonds:
            if not mcs_pair[0].is_cyclic and 'ar' in mcs_pair[0].sybyl_type:
                atom0 = mcs_pair[0].atoms[0]
                if atom_has_open_valency(atom0):
                    remove_atoms.append(atom0)
                atom1 = mcs_pair[0].atoms[1]
                if atom_has_open_valency(atom1):
                    remove_atoms.append(atom1)
        return remove_atoms


class SubstitutentComparer(object):

    def __init__(self):
        self.smarts_dict = {'amide': 'C(=O)N'}

    def compare(self, ligand, mcs_atoms):
        '''
        Return atoms that should be removed from strict MCS because their number of heavy atom neighbours is
        different from ligand. Also removes partially matched amide carbonyls.
        :param ligand:
        :param mcs_atoms:
        :return: list of atoms
        '''
        remove_atoms = []
        for substructure in self.smarts_dict:
            smarts = self.smarts_dict[substructure]
            searcher = search.SubstructureSearch()
            searcher.add_substructure(search.SMARTSSubstructure(smarts))
            ligand_hits = searcher.search(ligand)

            for hit in ligand_hits:
                ligand_match_atoms = hit.match_atoms()
                ligand_mcs_pairs = [a for a in mcs_atoms if a[1] in ligand_match_atoms]
                for mcs_pair in ligand_mcs_pairs:
                    scaffold_atom = mcs_pair[0]
                    heavy_scaffold_neighbours = [a for a in scaffold_atom.neighbours if a.atomic_symbol != 'H']
                    heavy_ligand_neighbours = [a for a in mcs_pair[1].neighbours if a.atomic_symbol != 'H']

                    if scaffold_atom.atomic_symbol == 'N':
                        # for partial amide matches remove entire amide, to enable its rotation
                        if len(heavy_scaffold_neighbours) != len(heavy_ligand_neighbours):
                            remove_atoms.append(scaffold_atom)
                            if len(ligand_mcs_pairs) == len(ligand_match_atoms):
                                remove_atoms.append(ligand_mcs_pairs[1][0])
                        # N(-C)-C will lead to symmetric matches. If both C have open valency, strict scaffold match can be off
                        elif len(heavy_scaffold_neighbours) == 3:
                            open_valency_neigbours = [n for n in heavy_scaffold_neighbours if atom_has_open_valency(n)]
                            if len(open_valency_neigbours) == 2:
                                remove_atoms.extend(open_valency_neigbours)
        return remove_atoms


# def _compare_aromatic_rings(scaffold, ligand, mcs_atoms):
#     smarts = '[aD3]1[aD2][aD2][aD2][aD2][aD2]1'
#     searcher = search.SubstructureSearch()
#     searcher.add_substructure(search.SMARTSSubstructure(smarts))
#     scaffold_hits = searcher.search(scaffold)
#     ligand_hits = searcher.search(ligand)
#
#     remove_atoms = []
#     for hit in scaffold_hits:
#         scaffold_match_atoms = hit.match_atoms()[1:]
#         scaffold_mcs_pairs = [a for a in mcs_atoms if a[0] in scaffold_match_atoms]
#
#         ligand_mcs_pairs = []
#         for ligand_hit in ligand_hits:
#             ligand_match_atoms = ligand_hit.match_atoms()[1:]
#             ligand_mcs_pairs = [a for a in mcs_atoms if a[1] in ligand_match_atoms]
#
#         if scaffold_mcs_pairs == ligand_mcs_pairs:
#             continue
#         elif scaffold_mcs_pairs:
#             remove_atoms = remove_atoms + [a[0] for a in scaffold_mcs_pairs]
#
#     return remove_atoms


class GeometryFlag(object):
    def __init__(self, smarts, atom_index_1, atom_index_2, min_distance, has_to_be_los=True):
        self.smarts = smarts
        self.substructure = search.SMARTSSubstructure(smarts)
        self.atom_index_1 = atom_index_1
        self.atom_index_2 = atom_index_2
        self.min_distance = min_distance
        self.has_to_be_los = has_to_be_los


class ConformerPruner(object):
    def __init__(self):
        self.flagged_geometries = self._flagged_geometries()

    def _flagged_geometries(self):
        flagged_geometries = []
        flagged_geometries.append(GeometryFlag('N(H)-C(=O)-!@[#6]@[#6]-!@N(H)', 1, 7, 2.2))
        flagged_geometries.append(GeometryFlag('N(H)-C(=O)-!@[#6](~[#7D2H0])@[#6]-!@N(H)', 1, 8, 2.5))
        flagged_geometries.append(GeometryFlag('N(H)-C(=O)-!@[#6](~[#7D2H0])@[#6]-!@N(H)', 3, 6, 2.9))
        flagged_geometries.append(GeometryFlag('N(H)-C(=O)-!@[#6](~[#7D2H0])@[#6]-!@N(H)', 1, 7, 3.0))
        flagged_geometries.append(GeometryFlag('N(H)-C(=O)-!@[#6]@[#6]-!@N(H)', 0, 6, 3.0))
        flagged_geometries.append(GeometryFlag('N(H)-C(=O)-!@[#6](~[#6H1])~[#6]~[#7D2H0]', 3, 7, 3.8))
        flagged_geometries.append(GeometryFlag("N([H])-C(=O)-!@[#6]@[#7]([H])", 1, 6, 2.1))
        flagged_geometries.append(GeometryFlag('C(=O)-!@[#6]@[#6]-!@N-C(=O)', 1, 6, 2.5))
        flagged_geometries.append(GeometryFlag('[#6](=O)-!@[c]([cD2])@[c]([cD2])-O-H', 1, 6, 3))
        flagged_geometries.append(GeometryFlag('[#6](=O)-[c]([cD3][#6])@[c]([cD2])-O-H', 1, 7, 3))
        flagged_geometries.append(GeometryFlag('[NH1][C](=O)!@[#6]~[#6][N-]', 2, 5, 3.5))
        flagged_geometries.append(GeometryFlag('[#7D2H0]~[#6](~[#6])-!@[#6](~[#6])~[#7D2H0]', 0, 5, 3))
        flagged_geometries.append(GeometryFlag('[#7D2H0]~[#6](~[#6][H])-!@[#6](~[#6][H])~[#7D2H0]', 3, 6, 2.5))
        flagged_geometries.append(GeometryFlag('C(=O)-[#6]@[#7D2H0]', 1, 3, 3.0, False))
        flagged_geometries.append(GeometryFlag('C(=O)-N-[#6]@[#7D2H0]', 1, 4, 3.0))

        return flagged_geometries

    def is_bad_conformer(self, ligand):
        for flagged_geometry in self.flagged_geometries:
            searcher = search.SubstructureSearch()
            searcher.add_substructure(flagged_geometry.substructure)
            hits = searcher.search(ligand)
            for hit in hits:
                match_atoms = hit.match_atoms()
                is_los = match_atoms[flagged_geometry.atom_index_1].is_in_line_of_sight(
                    match_atoms[flagged_geometry.atom_index_2])
                distance = descriptors.MolecularDescriptors.atom_distance(match_atoms[flagged_geometry.atom_index_1],
                                                                          match_atoms[flagged_geometry.atom_index_2])
                if distance < flagged_geometry.min_distance:
                    if flagged_geometry.has_to_be_los:
                        if is_los:
                            return True
                    else:
                        return True
        return False


def update_gold_conf(gold_conf, water_paths=False, fixed_bonds=None, args=None):
    f = open(gold_conf, 'r')
    newdata = f.read()
    f.close()

    newdata = newdata.replace("internal_ligand_h_bonds = 0", "internal_ligand_h_bonds = 1")
    newdata = newdata.replace("rms_tolerance = 1.5", "rms_tolerance = 0.5")
    newdata = newdata.replace('save_lone_pairs = 1', 'save_lone_pairs = 0')

    newdata = newdata.replace('pt_crosswt = 0', 'pt_crosswt = 95')
    newdata = newdata.replace('allele_mutatewt = 0', 'allele_mutatewt = 95')
    newdata = newdata.replace('migratewt = 0', 'migratewt = 10')

    newdata = newdata.replace('solvate_all = 1',
                              f'solvate_all = 1\n' + fixed_bonds)
    newdata = newdata.replace('param_file = DEFAULT',
                              f'param_file = /rf_scoring/gold_custom.params')
    newdata = newdata.replace('rescore_param_file = DEFAULT',
                              f'rescore_param_file = /rf_scoring/gold_custom.params')

    if water_paths:
        newdata = newdata + '  WATER DATA\n'
        for water_file in water_paths:
            newdata = newdata + f'water 1 toggle trans_spin 0.5 {water_file}\n'

    f = open(gold_conf, 'w')
    f.write(newdata)
    f.close()


def is_in_ring_with_open_valency(atom, mol):
    '''
    Check if an atom belongs to a "flat" ring system that has an open valency
    :param atom: 
    :param mol: 
    :return: 
    '''

    if atom not in mol.atoms:
        return Exception('Atom is not part of Molecule.')
    if atom.is_cyclic:
        for ring in mol.rings:
            # Don't rotate saturate rings larger than 5, as their conformation is not sampled
            if ring.is_aromatic or ring.is_fully_conjugated or len(ring.atoms) == 5:
                if atom in ring.atoms:
                    for ring_atom in ring.atoms:
                        if atom_has_open_valency(ring_atom):
                            print('open valency')
                            return True
                        if ring.is_fused:
                            for ring_atom in ring.atoms:
                                fused_rings = [fused_ring for fused_ring in ring_atom.rings if fused_ring != ring]
                                if len(fused_rings) > 1:
                                    for fused_ring in fused_rings:
                                        for fused_ring_atom in fused_ring.atoms:
                                            if atom_has_open_valency(fused_ring_atom):
                                                print('Fused ring open valency')
                                                return True
    return False


def atom_has_open_valency(atom):
    num_neighbours = len(atom.neighbours)
    if atom.sybyl_type == 'C.3' and num_neighbours < 4:
        return True
    if atom.sybyl_type == 'C.2' and num_neighbours < 3:
        return True
    if atom.sybyl_type == 'N.pl3' and num_neighbours < 3:
        return True
    if atom.sybyl_type == 'N.3' and atom.formal_charge == 0 and num_neighbours < 3:
        return True
    if atom.sybyl_type == 'N.am' and atom.formal_charge == 0 and num_neighbours < 3:
        return True
    if atom.sybyl_type == 'N.3' and atom.formal_charge == 1 and num_neighbours < 4:
        return True
    if atom.sybyl_type == 'N.4' and num_neighbours < 4:
        return True
    return False


def _mcs_metrics(docked_ligand_entry, scaffold):
    docked_ligand = docked_ligand_entry.molecule
    docked_ligand.set_formal_charges()
    rdkit_ligand = Chem.MolFromMolBlock(docked_ligand.to_string('sdf'), removeHs=False)
    docked_ligand = docked_ligand.from_string(Chem.MolToMolBlock(rdkit_ligand))

    docked_ligand.remove_unknown_atoms()
    mcs_searcher = descriptors.MolecularDescriptors.MaximumCommonSubstructure()
    mcs_searcher.settings.check_bond_type = False
    mcs_searcher.settings.ignore_hydrogens = True
    mcs_searcher.settings.check_bond_count = False
    temp_docked_ligand = docked_ligand.copy()
    temp_scaffold = scaffold.copy()
    mcs_atoms_docked_template = mcs_searcher.search(temp_docked_ligand, temp_scaffold, search_step_limit=1000000)[0]

    # Keep only MCS part to allow RMSD calculation to consider symmetry
    ligand_mcs_atoms = [mcs_pair[0] for mcs_pair in mcs_atoms_docked_template]
    remove_atoms = [a for a in temp_docked_ligand.atoms if a not in ligand_mcs_atoms]
    temp_docked_ligand.remove_atoms(remove_atoms)
    temp_docked_ligand.remove_hydrogens()

    scaffold_mcs_atoms = [mcs_pair[1] for mcs_pair in mcs_atoms_docked_template]
    remove_atoms = [a for a in temp_scaffold.atoms if a not in scaffold_mcs_atoms]
    temp_scaffold.remove_atoms(remove_atoms)
    temp_scaffold.remove_hydrogens()

    rmsd = descriptors.MolecularDescriptors.rmsd(temp_docked_ligand, temp_scaffold,
                                                 # atoms=mcs_atoms_docked_template,
                                                 exclude_hydrogens=True, with_symmetry=True
                                                 )
    docked_ligand_entry.attributes['RMSD_to_mcs'] = rmsd
    abs_mcs_size = len(scaffold.heavy_atoms)
    rel_mcs_size = abs_mcs_size / len(docked_ligand.heavy_atoms)
    docked_ligand_entry.attributes['rel_mcs_size'] = rel_mcs_size
    docked_ligand_entry.attributes['abs_mcs_size'] = abs_mcs_size
    return docked_ligand_entry, docked_ligand


def _select_best_soln(docking_folder, scaffold, row, srn, pdb_id, attributes, scaffold_filename):
    docked_ligand_files = docking_folder.glob(f'gold_soln_docking_input_*_*_*.sdf')

    best_soln = entry.Entry()
    best_soln_fitness = -999
    for docked_ligand_file in docked_ligand_files:
        docked_ligand_entries = []
        docked_ligand_file = docked_ligand_file.resolve()

        with io.EntryReader(str(docked_ligand_file)) as rdr:
            for docked_ligand_entry in rdr:

                docked_ligand_entry, docked_ligand = _mcs_metrics(docked_ligand_entry, scaffold, )
                docked_ligand_entry.attributes['template_strucid'] = row['template_strucid']
                docked_ligand_entry.attributes['rel_mcs_size_to_native_ligand'] = row['rel_mcs_size_to_native_ligand']
                docked_ligand_entry.attributes['tanimoto_similiarity_to_native_ligand'] = row['similarity']
                docked_ligand_entry.attributes['is_decoy'] = False
                docked_ligand_entry.attributes.update(attributes)
                docked_ligand_entries.append(docked_ligand_entry)
                gold_rescore_fitness = float(docked_ligand_entry.attributes['Gold.PLP.Fitness'])
                if gold_rescore_fitness > best_soln_fitness and not ConformerPruner().is_bad_conformer(
                        docked_ligand_entry.molecule):
                    best_soln = docked_ligand_entry
                    best_soln_file = docked_ligand_file
                    best_soln_fitness = gold_rescore_fitness

        with io.EntryWriter(docked_ligand_file) as wr:
            for docked_ligand_entry in docked_ligand_entries:
                wr.write(docked_ligand_entry)

    shape_similarity = sp.check_output(
        [sys.executable, '/ccdc_roche_scoring/shape_similarity.py',
         '--template_ligand', scaffold_filename,
         '--docked_ligand', str(best_soln_file)], env=os.environ)

    shape_similarity = float(shape_similarity.decode('utf-8').split('\n')[0])
    print('shape', shape_similarity)

    best_soln.attributes['is_decoy'] = False
    best_soln.attributes['TanimotoCombo'] = shape_similarity
    best_soln.attributes.update(attributes)

    best_soln_mol = best_soln.molecule
    best_soln_mol.set_formal_charges()
    rdkit_ligand = Chem.MolFromMolBlock(best_soln_mol.to_string('sdf'), removeHs=False)
    best_soln_mol = best_soln_mol.from_string(Chem.MolToMolBlock(rdkit_ligand))
    best_soln_mol.remove_unknown_atoms()

    new_best_soln = entry.Entry.from_molecule(best_soln_mol)
    new_best_soln.attributes = best_soln.attributes

    best_soln_file = str(Path(docking_folder) / Path(f'best_soln_{pdb_id}_{srn}.sdf'))
    with io.EntryWriter(best_soln_file) as w:
        w.write(new_best_soln)

    gold_conf = str(list(Path(docking_folder).glob('*.conf'))[0])
    _write_best_soln_pocket(gold_conf, best_soln_file)
    return True


def _select_best_decoy(docking_folder, mcs_searcher, scaffold, srn, pdb_id, reference_ligand_file):
    docked_ligand_files = docking_folder.glob(f'gold_soln_docking_input_*_*_*.mol2')

    for docked_ligand_file in docked_ligand_files:
        docked_ligand_entries = []
        docked_ligand_file = docked_ligand_file.resolve()
        with io.EntryReader(str(docked_ligand_file)) as rdr:
            for docked_ligand_entry in rdr:
                docked_ligand_entry, docked_ligand = _mcs_metrics(docked_ligand_entry, scaffold)
                rmsd_to_mcs = docked_ligand_entry.attributes['RMSD_to_mcs']
                docked_ligand_entry.attributes['template_strucid'] = pdb_id
                docked_ligand_entries.append(docked_ligand_entry)

        shape_similarity = sp.check_output(
            [sys.executable, '/ccdc_roche_scoring/shape_similarity.py',
             '--template_ligand', str(reference_ligand_file),
             '--docked_ligand', str(docked_ligand_file)], env=os.environ)

        shape_similarity = float(shape_similarity.decode('utf-8').split('\n')[0])
        print('shape', shape_similarity)

        # set attributes for decoy
        with io.EntryReader(str(reference_ligand_file)) as rdr:
            attributes = rdr[0].attributes
        docked_ligand_entry.attributes.update(attributes)
        docked_ligand_entry.attributes['RMSD_to_mcs'] = rmsd_to_mcs
        docked_ligand_entry.attributes['TanimotoCombo_to_reference'] = shape_similarity
        docked_ligand_entry.attributes['is_decoy'] = True

        with io.EntryWriter(docked_ligand_file) as wr:
            for docked_ligand_entry in docked_ligand_entries:
                wr.write(docked_ligand_entry)
    return True


def _protein_preparation(pdb, dry_receptor_file, pdb_id, target):
    protein_res_labels = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET',
                          'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HOH', '']

    target_protein = protein.Protein.from_file(str(pdb))

    water_atoms = [a for a in target_protein.atoms if a.protein_atom_type == 'Water' and
                   is_important_water(a, target_protein, target) == False]
    ligand_atoms = [a for a in target_protein.atoms if
                    ('LIG1' in a.residue_label or a.residue_label[:3] not in protein_res_labels)]

    target_protein.remove_atoms(ligand_atoms + water_atoms)
    target_protein.standardise_aromatic_bonds()
    target_protein.standardise_delocalised_bonds()

    if not dry_receptor_file.parent.is_dir():
        dry_receptor_file.parent.mkdir()

    water_paths = list(dry_receptor_file.parent.glob('*water*.mol2'))

    if len(target_protein.waters) > 0 and len(water_paths) == 0:
        # water_paths = []
        for cnt, water in enumerate(target_protein.waters):
            water_path = (dry_receptor_file.parent / Path(f'{pdb_id}_water_{cnt}.mol2')).resolve()
            water_paths.append(water_path)
            with io.MoleculeWriter(water_path) as w_wr:
                water.add_hydrogens()
                w_wr.write(water)
    target_protein.remove_all_waters()

    # save as pdb first to ensure correct format
    temp_pdb = str(dry_receptor_file.parent / Path(dry_receptor_file.stem + '.pdb'))
    with io.MoleculeWriter(temp_pdb) as p_wr:
        p_wr.write(target_protein)

    target_protein = protein.Protein.from_file(temp_pdb)
    target_protein.add_hydrogens()
    with io.MoleculeWriter(str(dry_receptor_file)) as p_wr:
        p_wr.write(target_protein)

    return water_paths, target_protein


def is_important_water(atom, protein, target='pde-10'):
    if target == 'pde-10':
        distance_searcher = descriptors.MolecularDescriptors.AtomDistanceSearch(protein)
        close_atoms = distance_searcher.atoms_within_range(atom.coordinates, 3.7)

    else:
        return False

    if target == 'pde-10':
        tyr_contact = False
        gln_contact = False
        trp_contact = False
        asp1_contact = False
        asp2_contact = False

        for at in close_atoms:
            res_label = at.residue_label[0:3]
            atom_label = at.label
            if res_label == 'TYR' and atom_label == 'OH':
                tyr_contact = True
                continue

            if res_label == 'TRP' and atom_label == 'NE1':
                trp_contact = True
                continue

            if res_label == 'GLN' and atom_label == 'OE1':
                gln_contact = True
                continue

            if res_label == 'ASP' and atom_label == 'O':
                asp1_contact = True
                continue

            if res_label == 'ASP' and atom_label == 'OD1':
                asp2_contact = True
                continue

        if gln_contact and trp_contact and tyr_contact:
            return True
        elif tyr_contact and asp1_contact and asp2_contact:
            return True

        else:
            return False


def _fix_rotatable_bond(mcs_scaffold_bond, strict_scaffold):
    '''
    Returns True if rotatable bond should be fixed in docking.
    :param mcs_scaffold_bond:
    :param strict_scaffold:
    :return:
    '''
    if atom_has_open_valency(mcs_scaffold_bond.atoms[0]) or atom_has_open_valency(mcs_scaffold_bond.atoms[1]):
        return False
    if is_in_ring_with_open_valency(mcs_scaffold_bond.atoms[0], strict_scaffold):
        return False
    if is_in_ring_with_open_valency(mcs_scaffold_bond.atoms[1], strict_scaffold):
        return False
    return True


def _dock(docker, dry_receptor_file, target_protein: protein.Protein, native_ligand, ligand_filename, docking_folder,
          pdb_id, srn, water_paths, scaffold=False, strict_scaffold=False, reference_ligand_file=False,
          diverse_solutions=True, args=None) -> None:
    settings = docker.settings
    settings.fitness_function = 'plp'

    if diverse_solutions:
        settings.diverse_solutions = True, 1, 0.5

    settings.add_protein_file(str(dry_receptor_file.resolve()))
    native_ligand.add_hydrogens(mode='missing')
    settings.binding_site = settings.BindingSiteFromLigand(target_protein, native_ligand, 10.0)

    flexible_residues = []
    for flexible_residue in args.flexible_residues:
        res = [r for r in settings.proteins[0].residues if flexible_residue in r.identifier][0]
        flexible_residues.append(res)
    for flexible_residue in flexible_residues:
        rl = settings.RotamerLibrary(settings.protein_files[0], flexible_residue)
        rl.add_default_rotamers()
        settings.add_rotamer_library(settings.proteins[0], rl)

    if scaffold:
        scaffold.add_hydrogens(mode='missing')
        scaffold.assign_bond_types()
        settings.add_constraint(settings.TemplateSimilarityConstraint('all', scaffold, weight=40))

    fixed_bonds = ''

    if strict_scaffold:
        settings.add_constraint(settings.ScaffoldMatchConstraint(strict_scaffold, weight=1000))
        ligand_mol = io.MoleculeReader(str(ligand_filename))[0]
        mcs_searcher = descriptors.MolecularDescriptors.MaximumCommonSubstructure()
        mcs_searcher.settings.ignore_hydrogens = True
        mcs_searcher.settings.check_bond_type = False
        mcs_atoms, mcs_bonds = mcs_searcher.search(ligand_mol, strict_scaffold)

        mcs_ligand_bonds = [b[0] for b in mcs_bonds]
        ligand_rotatable_bonds = [b for b in ligand_mol.bonds if b.is_rotatable and b not in mcs_ligand_bonds]

        # fix bonds
        for mcs_ligand_bond, mcs_scaffold_bond in mcs_bonds:
            if mcs_ligand_bond.is_rotatable:
                if _fix_rotatable_bond(mcs_scaffold_bond, strict_scaffold):
                    a0 = mcs_ligand_bond.atoms[0].index + 1
                    a1 = mcs_ligand_bond.atoms[1].index + 1
                    fixed_bonds = fixed_bonds + f'fix_rotatable_bond = {a0} {a1}\n'
                else:
                    ligand_rotatable_bonds.append(mcs_ligand_bond)

    autoscale = AutoScaleSetter(ligand_rotatable_bonds).autoscale
    settings.autoscale = autoscale
    ndocks = 12
    if autoscale == 0:
        settings._settings.set_niche_size(100)
        settings._settings.set_population_size(100)
        settings._settings.set_maxops(100000)
        settings._settings.set_n_islands(5)
        settings._settings.set_selection_pressure(1.1)
        ndocks = 6

    if reference_ligand_file:
        settings.reference_ligand_file = str(reference_ligand_file.resolve())

    settings.add_ligand_file(str(ligand_filename.resolve()), ndocks=ndocks)
    settings.output_file = str(docking_folder / Path(f'./{pdb_id}_docked_ligands_{srn}.sdf'))

    gold_conf = docking_folder / Path(f'gold_{pdb_id}_{srn}.conf')
    settings.make_absolute_file_names(str(gold_conf.resolve()))
    settings.write_options = ['NO_LINK_FILES', 'NO_RNK_FILES', 'NO_PLP_MOL2_FILES', 'NO_BESTRANKING_LST_FILE',
                              'NO_GOLD_LIGAND_MOL2_FILE', 'NO_LOG_FILES', 'NO_FIT_PTS_FILES']
    settings._settings.set_flip_planar_N(False)

    settings.write(str(gold_conf.resolve()))

    update_gold_conf(gold_conf, water_paths, fixed_bonds, args=args)

    docker = Docker()
    docker.settings = docker.settings.from_file(str(gold_conf.resolve()))
    docker.dock(file_name=str(gold_conf.resolve()))
    print('GOLD docking finished.')

    return


def _return_cis_trans_bicycle(atom):
    atom_hydrogen = [a for a in atom.neighbours if a.is_cyclic == False][0]
    for neigh in atom.neighbours:
        if (len(neigh.rings) == len(atom.rings)):
            if len(neigh.neighbours) < 4:
                return False
            neigh_hydrogen = [a for a in neigh.neighbours if a.is_cyclic == False][0]
            if descriptors.MolecularDescriptors.atom_distance(neigh_hydrogen, atom_hydrogen) < 2.8:
                return 'cis'
            else:
                return 'trans'


def _bicycle_cis_trans_mismatch(atom1, atom2):
    if atom1.sybyl_type.split('.')[1] in ['2', 'ar']:
        return False
    if atom2.sybyl_type.split('.')[1] in ['2', 'ar']:
        return False
    if atom1.atomic_symbol != 'C':
        return False
    if atom2.atomic_symbol != 'C':
        return False
    atom1_stereo = _return_cis_trans_bicycle(atom1)
    atom2_stereo = _return_cis_trans_bicycle(atom2)

    if atom1_stereo != atom2_stereo:
        return True
    else:
        return False


def _bicyclic_cis_trans_chirality(atom1, atom2):
    '''Check if the bicyclic system has cis-trans chirality.'''
    if 1 < len(atom1.rings) == len(atom2.rings) and atom1.is_spiro == False and atom2.is_spiro == False:
        atom1_second_bridge_atom = [n for n in atom1.neighbours if len(n.rings) > 1][0]
        atom2_second_bridge_atom = [n for n in atom2.neighbours if len(n.rings) > 1][0]
        if len(atom1_second_bridge_atom.neighbours) == len(atom2_second_bridge_atom.neighbours) == 4:
            return True
        else:
            return False
    else:
        return False


def _return_stereocenter_hybridization_dict(rdkit_mol):
    stereo_center_dict = {}
    hybridization_dict = {}
    for a in rdkit_mol.GetAtoms():
        hybridization_dict[a.GetProp('_TriposAtomName')] = a.GetHybridization()
        if '_CIPCode' in a.GetPropsAsDict().keys():
            stereo_center_dict[a.GetProp('_TriposAtomName')] = a.GetProp('_CIPCode')
    return stereo_center_dict, hybridization_dict


def _compare_mcs_stereo_chemistry(scaffold, ligand, mcs_atoms) -> bool:
    '''
    :return: True if MCS have the same stereochemistry, else return False
    '''
    # rdkit_scaffold
    rdkit_scaffold = Chem.MolFromMol2Block(scaffold.to_string('mol2'), removeHs=False)
    params = Chem.RemoveHsParameters()
    params.removeDegreeZero = True
    rdkit_scaffold = Chem.RemoveHs(rdkit_scaffold, params)

    # rdkit_ligand
    rdkit_ligand = Chem.MolFromMol2Block(ligand.to_string('mol2'))

    Chem.rdmolops.FindPotentialStereo(rdkit_scaffold)
    Chem.rdmolops.FindPotentialStereo(rdkit_ligand)

    scaffold_atom_stereo, scaffold_atom_hybridization = _return_stereocenter_hybridization_dict(rdkit_scaffold)
    ligand_atom_stereo, ligand_atom_hybridization = _return_stereocenter_hybridization_dict(rdkit_ligand)

    scaffold_atoms_stereo_mismatch = []
    hybridization_mismatch = []
    for mcs_atom_pair in mcs_atoms:
        if mcs_atom_pair[0].atomic_symbol == 'H':
            continue
        ligand_stereo = None
        scaffold_stereo = None
        scaffold_atom = mcs_atom_pair[0]
        ligand_atom = mcs_atom_pair[1]
        scaffold_atom_label = scaffold_atom.label
        ligand_atom_label = ligand_atom.label

        if scaffold_atom_label in scaffold_atom_stereo.keys():
            scaffold_stereo = scaffold_atom_stereo[scaffold_atom_label]
        if ligand_atom_label in ligand_atom_stereo.keys():
            ligand_stereo = ligand_atom_stereo[ligand_atom_label]
        if scaffold_stereo is not None and ligand_stereo is not None and scaffold_stereo != ligand_stereo:
            ligand_neighbours = sorted(a.atomic_symbol for a in list(
                itertools.chain.from_iterable(n.neighbours for n in ligand_atom.neighbours)))
            scaffold_neighbours = sorted(a.atomic_symbol for a in
                                         list(itertools.chain.from_iterable(
                                             n.neighbours for n in scaffold_atom.neighbours)))
            if len(ligand_neighbours) > len(scaffold_neighbours):
                scaffold_atoms_stereo_mismatch.append(scaffold_atom)
            elif ligand_neighbours == scaffold_neighbours:
                scaffold_atoms_stereo_mismatch.append(scaffold_atom)
        if ligand_atom_hybridization[ligand_atom_label] != scaffold_atom_hybridization[scaffold_atom_label]:
            hybridization_mismatch.append(scaffold.atom(scaffold_atom_label))
    return list(set(scaffold_atoms_stereo_mismatch + hybridization_mismatch))


def _mcs_templates_df(native_ligand_entries, ligand_mol, series_template_strucids=None):
    '''
    Find maximum common substructure. Eliminate ring and hybridization mismatch atoms.
    :param native_ligand_entries:
    :param mcs_searcher:
    :param ligand_mol:
    :return:
    Example: Ligand RO6898508, template: 1qgx
    #>>> abs_mcs_size
    #27
    '''
    templates_df = pd.DataFrame()
    templates = {'template_strucid': [], 'abs_mcs_size': [], 'rel_mcs_size_to_ligand': [],
                 'rel_mcs_size_to_native_ligand': [], 'native_ligand': [], 'scaffold': [],
                 'mcs_atom_labels': [], 'mcs_object': []}
    for cnt, native_ligand_entry in enumerate(native_ligand_entries):

        if series_template_strucids is not None and native_ligand_entry.attributes[
            'STRUCID'] not in series_template_strucids:
            continue

        # if native_ligand_entry.attributes['STRUCID'] == '1qhdw':
        #     with io.MoleculeWriter('test_def.sdf') as w:
        #         w.write(native_ligand_entry.molecule)
        #     print('nice')

        mcs = MCS(ligand_mol, native_ligand_entry.molecule)
        scaffold = mcs.mcs_scaffold
        mcs_atoms = mcs.mcs_atoms
        mcs_atom_labels = [(mcs_pair[0].label, mcs_pair[1].label) for mcs_pair in mcs_atoms if
                           mcs_pair[0] in scaffold.atoms]

        mcs_size = len(scaffold.heavy_atoms)
        templates['template_strucid'].append(native_ligand_entry.attributes['STRUCID'])
        templates['abs_mcs_size'].append(mcs_size)
        temp_ligand = ligand_mol.copy()
        templates['rel_mcs_size_to_ligand'].append(mcs_size / len(temp_ligand.heavy_atoms))
        templates['rel_mcs_size_to_native_ligand'].append(mcs_size / len(native_ligand_entry.molecule.heavy_atoms))
        templates['native_ligand'].append(native_ligand_entry)
        templates['scaffold'].append(scaffold)
        templates['mcs_object'].append(mcs)
        templates['mcs_atom_labels'].append(mcs_atom_labels)
        templates_df = pd.DataFrame(templates).sort_values(by='abs_mcs_size', ascending=False).reset_index(drop=True)

    # Tanimoto similarity for highest MCS compounds
    searcher = search.SimilaritySearch(ligand_mol)
    for index, row in templates_df.iterrows():
        similarity = searcher.search_molecule(row['native_ligand'].molecule).similarity
        templates_df.loc[index, 'similarity'] = similarity

    templates_df.loc[templates_df['rel_mcs_size_to_ligand'] < 0.5, 'similarity'] = 0
    templates_df = templates_df.sort_values(by=['abs_mcs_size', 'similarity'], ascending=False).reset_index(drop=True)

    return templates_df


def _write_starting_ligand(ligand_mol, ligand_filename, docking_folder, scaffold, scaffold_filename):
    rdkit_scaffold = Chem.MolFromMol2Block(scaffold.to_string(), removeHs=False)
    w = Chem.SDWriter(str(scaffold_filename))
    w.write(rdkit_scaffold)
    w.close()

    tmp_ligand_file = 'tmp_ligand_mol.sdf'
    tmp_ligand_file = docking_folder / tmp_ligand_file
    ligand_mol.add_hydrogens(mode='missing')
    rdkit_ligand = Chem.MolFromMol2Block(ligand_mol.to_string('mol2'), removeHs=False)
    w = Chem.SDWriter(str(tmp_ligand_file))
    w.write(rdkit_ligand)
    w.close()

    # oechem is imcompatible with ccdc
    sp.check_output([sys.executable, '/ccdc_roche_scoring/template_alignment.py',
                     '--ligand',
                     str(Path(tmp_ligand_file).resolve()),
                     '--scaffold',
                     str(Path(scaffold_filename).resolve()),
                     '--output',
                     str(Path(ligand_filename).resolve())], env=os.environ)

    # Ensure stereo centers were preserved by omega:
    rdkit_ligand_stereo_dict = _return_stereocenter_hybridization_dict(rdkit_ligand)[0]
    input_ligand = Chem.MolFromMol2File(str(Path(ligand_filename).resolve()))
    input_ligand_stereo_dict = _return_stereocenter_hybridization_dict(input_ligand)[0]
    ccdc_input_ligand = io.MoleculeReader(str(Path(ligand_filename).resolve()))[0]
    mcs_searcher = descriptors.MolecularDescriptors.MaximumCommonSubstructure()
    mcs_searcher.settings.ignore_hydrogens = True
    mcs_atoms = mcs_searcher.search(ligand_mol, ccdc_input_ligand, search_step_limit=1000000)[0]
    stereo_center_pairs = [a for a in mcs_atoms if a[0].label in rdkit_ligand_stereo_dict.keys()]
    for stereo_center_pair in stereo_center_pairs:
        if rdkit_ligand_stereo_dict[stereo_center_pair[0].label] != input_ligand_stereo_dict[
            stereo_center_pair[1].label]:
            print('Omega returned wrong stereo center')
            return False
    return True


def gold_decoy_docking(reference_docking_job, args=None):
    from ccdc_roche_scoring import join_docked_rf_counts
    mcs_searcher = descriptors.MolecularDescriptors.MaximumCommonSubstructure()
    mcs_searcher.settings.ignore_hydrogens = True

    df = pd.read_csv(Path(reference_docking_job) / Path('docked_rf_count_df.csv'))
    best_solns = join_docked_rf_counts.get_best_docking_solutions(df, args.target)[0]
    best_solns = [Path(ligand_file).parent.stem for ligand_file in best_solns['ligand_file'].values]
    reference_docking_dirs = [reference_dir for reference_dir in Path(reference_docking_job).glob('*_*') if
                              reference_dir.stem in best_solns]
    for reference_docking_dir in reference_docking_dirs:
        if reference_docking_dir.is_dir():

            docking_folder = Path(reference_docking_dir.name)
            if not docking_folder.is_dir():
                docking_folder.mkdir()
            pdb_id = docking_folder.name.split('_')[0].lower()
            srn = docking_folder.name.split('_')[1]
            # if 'RO6858069' not in srn:
            #     continue

            reference_ligand_file = list(reference_docking_dir.glob('best_soln_*.sdf'))[0].resolve()

            # load scaffold file
            scaffold_file = list(reference_docking_dir.glob('scaffold_*.sdf'))[0].resolve()
            with io.MoleculeReader(str(scaffold_file)) as rdr:
                scaffold = rdr[0]

            # load cavity file
            cavity_ligand_filepath = list(reference_docking_dir.glob('cavity_*.mol2'))[0].resolve()
            with io.MoleculeReader(str(cavity_ligand_filepath)) as rdr:
                cavity_ligand = rdr[0]

            # load target file
            dry_receptor_file = (reference_docking_dir / Path('gold_protein.mol2')).resolve()
            with io.EntryReader(str(dry_receptor_file)) as rdr:
                target_protein = protein.Protein.from_entry(rdr[0])
                # target_protein = rdr[0]

            # run docking
            docker = Docker()
            ligand_filepath = list(reference_docking_dir.glob('docking_input_*.mol2'))[0].resolve()
            _dock(docker, dry_receptor_file, target_protein, cavity_ligand, ligand_filepath, docking_folder, pdb_id,
                  srn, water_paths=False, scaffold=False, reference_ligand_file=reference_ligand_file, autoscale=25,
                  args=args)

            _select_best_decoy(docking_folder, mcs_searcher, scaffold, srn, pdb_id, reference_ligand_file)


def _setup_docking(docking_folder, targets, pdb_id, project_home, ligand_mol, srn, scaffold, scaffold_filename,
                   native_ligand, args, strict_scaffold):
    if not docking_folder.is_dir():
        docking_folder.mkdir()

    # get target structure and water
    dry_receptor_file = Path(project_home) / Path('apo_proteins') / Path(f'{pdb_id}_dry.mol2')
    if not Path(str(dry_receptor_file)).is_file():
        pdb = [p for p in targets if pdb_id in str(p)][0]
        print('Preparing protein MOL2 files...')
        _protein_preparation(pdb, dry_receptor_file, pdb_id, target=args.target)
    with io.EntryReader(str(dry_receptor_file)) as t_rdr:
        target_protein = protein.Protein.from_entry(t_rdr[0])
        water_paths = [str(water_path.resolve()) for water_path in
                       dry_receptor_file.parent.glob(f'{pdb_id}_water*.mol2')]
        if len(water_paths) == 0:
            print('No water molecules selected for docking...')
            water_paths = False
    target_protein.identifier = pdb_id

    ligand_filename = Path(docking_folder) / Path(f'docking_input_{srn}.mol2')  # mol2 to preserve atom types

    if _write_starting_ligand(ligand_mol, ligand_filename, docking_folder,
                              strict_scaffold,
                              scaffold_filename):
        docker = Docker()
        diverse_solutions = True

        _dock(docker, dry_receptor_file, target_protein, native_ligand, ligand_filename,
              docking_folder, pdb_id, srn, water_paths, scaffold=scaffold,
              strict_scaffold=strict_scaffold, diverse_solutions=diverse_solutions, args=args)


class AutoScaleSetter(object):
    '''
    Set autoscale to 0 if there are very few rotatable bonds.
    '''

    def __init__(self, ligand_rotatable_bonds):
        self.ligand_rotatable_bonds = [b for b in ligand_rotatable_bonds if b.sybyl_type != 'am']
        self.autoscale = self._set_autoscale()

    def _bond_rotates_heavy_atoms(self, bond):
        for atom in bond.atoms:
            num_small_neighbours = len([n for n in atom.neighbours if n.atomic_symbol in ['H', 'F']])
            if len(atom.neighbours) - num_small_neighbours == 1:
                return False
        else:
            return True

    def _set_autoscale(self):

        ligand_rotatable_bonds = [b for b in self.ligand_rotatable_bonds if self._bond_rotates_heavy_atoms(b)]

        if len(ligand_rotatable_bonds) > 4:
            autoscale = 80
        else:
            autoscale = 0

        return autoscale


def _write_best_soln_pocket(gold_conf, ligand_file):
    docking_settings = Docker.Settings().from_file(gold_conf)
    docking_results = Docker.Results(docking_settings)
    csd_ligand_entry = docking_results.DockedLigandReader(ligand_file, docking_settings)[0]

    # setup protein-ligand-complex
    protein = docking_results.make_complex(csd_ligand_entry)
    apo_protein = protein.copy()
    apo_protein.remove_ligand(apo_protein.ligands[0].identifier)
    apo_protein.remove_unknown_atoms()
    apo_protein.identifier = csd_ligand_entry.identifier

    apo_pocket = protein.copy()
    apo_pocket = _cut_out_binding_site_by_distance(apo_pocket, csd_ligand_entry.molecule)
    apo_pocket.remove_unknown_atoms()
    apo_pocket.identifier = csd_ligand_entry.identifier

    pocket_file = Path(ligand_file).parent / 'best_soln_pocket.mol2'
    protein_soln_file = Path(ligand_file).parent / 'best_soln_protein.pdb'
    with io.EntryWriter(pocket_file) as w:
        w.write(apo_pocket)

    apo_protein.kekulize()
    with io.EntryWriter(protein_soln_file) as w:
        w.write(apo_protein)
    return


def gold_scaffold_docking(ligands_for_alignment, project_home, args=None):
    rel_mcs_size_to_ligand_threshold = 0.5
    rel_mcs_size_to_native_ligand_threshold = 0.3
    similarity_threshold = 0.7
    pdb_ligand_df_file = Path(project_home, 'pdb_ligand.csv')

    if args.target == 'default':
        targets = list(((Path(project_home) / Path('tmp_aligned_for_MOE_sanitized')).glob('*.pdb')))
        template_ligands_sdf = str(
            Path(project_home) / 'tmp_aligned_3d_sdf_sanitized/ligand_templates_for_mcs_manual.sdf')
        if pdb_ligand_df_file.is_file():
            pdb_ligand_df = pd.read_csv(pdb_ligand_df_file)
            pdb_ligand_df = pdb_ligand_df[pdb_ligand_df['is_high_quality'] == True]
            pdb_ligand_df = pdb_ligand_df[pdb_ligand_df['template_file'].isna() == False]

    series_df = Path('../../series_assignment.csv')
    if series_df.is_file():
        series_df = pd.read_csv(series_df)
        series_df = series_df[[c for c in series_df.columns if c in ['Proasis ID', 'SRN', 'series']]].rename(
            {'Proasis ID': 'strucid'}, axis=1)
        pdb_ligand_df = pdb_ligand_df.join(series_df.set_index('strucid'), on='strucid')
    else:
        series_df = pd.DataFrame()

    native_ligand_entries = []
    with io.EntryReader(template_ligands_sdf) as sdf_rdr:
        for cnt, native_ligand_entry in enumerate(sdf_rdr):
            if 'strucid' in native_ligand_entry.attributes:
                strucid = native_ligand_entry.attributes['strucid']
            elif args.target != 'default':
                strucid = native_ligand_entry.identifier.split('.')[0]
            elif pdb_ligand_df_file.is_file() and args.target == 'default':
                ligand_file = str(Path(native_ligand_entry.attributes['$File']).resolve())
                strucid = pdb_ligand_df[pdb_ligand_df['template_file'] == ligand_file]['strucid']
                if strucid.shape[0] == 1:
                    strucid = strucid.to_list()[0]
                else:
                    strucid = False
            else:
                strucid = native_ligand_entry.identifier.split('.')[0]

            native_ligand = native_ligand_entry.molecule
            if args.target != 'default':
                max_comp_size = 0
                max_comp = None
                for native_ligand_comp in native_ligand.components:
                    comp_size = len(native_ligand_comp.heavy_atoms)
                    if comp_size > max_comp_size:
                        max_comp_size = comp_size
                        max_comp = native_ligand_comp
                native_ligand = max_comp
                native_ligand.add_hydrogens(mode='missing')

            native_ligand.normalise_labels()
            new_native_ligand_entry = entry.Entry.from_molecule(native_ligand)
            new_native_ligand_entry.attributes = native_ligand_entry.attributes
            native_ligand_entry = new_native_ligand_entry
            if strucid:
                native_ligand_entry.attributes['STRUCID'] = strucid
                native_ligand_entries.append(native_ligand_entry)

    ligands_for_alignment_sdf = Path(project_home) / Path(ligands_for_alignment)
    with io.EntryReader(str(ligands_for_alignment_sdf)) as lig_sd_rdr:
        for lig_cnt, ligand in enumerate(lig_sd_rdr):
            # try:
            ligand_entry = lig_sd_rdr[lig_cnt]
            if 'SRN' in ligand_entry.attributes:
                srn = ligand_entry.attributes['SRN'].replace(' ', '')
            else:
                srn = ligand_entry.identifier.replace(' ', '')

            print(srn)
            series_template_strucids = None
            if series_df.shape[0] > 1:
                series = series_df[series_df['SRN'] == srn]['series'].values[0]
                if series == series:
                    series_template_strucids = pdb_ligand_df[pdb_ligand_df['series'] == series]['strucid'].values
                else:
                    continue

            ligand_mol = ligand_entry.molecule.components[0]

            # kekulize with RDKit because it handles aromaticity better
            rdkit_ligand = Chem.MolFromMolBlock(ligand_mol.to_string('sdf'), removeHs=False)
            ligand_mol = ligand_mol.from_string(Chem.MolToMolBlock(rdkit_ligand))

            ligand_mol.normalise_labels()
            templates_df = _mcs_templates_df(native_ligand_entries, ligand_mol, series_template_strucids)

            # Remove template structures that were solved after the assay date
            if args.target == 'pde-10':
                proasis_df = pd.read_csv('../../proasis_date.csv').astype({'DATESOLVED': 'datetime64'})
                templates_df = templates_df.join(proasis_df.set_index('STRUCID'), on='template_strucid')
                assay_date = ligand_entry.attributes[
                    'PDE10_FULL_LEN_CGMP_SPA_IC50_h-PDE10A(14-779)-E.Coli-c: AP005128;Min;EXP Test Date']
                assay_date = pd.to_datetime(assay_date)
                templates_df = templates_df[templates_df['DATESOLVED'] < assay_date]

            dock = False
            for index, row in templates_df.iterrows():
                criterion1 = row['rel_mcs_size_to_ligand'] >= rel_mcs_size_to_ligand_threshold and \
                             row['rel_mcs_size_to_native_ligand'] >= rel_mcs_size_to_native_ligand_threshold
                criterion2 = row['abs_mcs_size'] >= 10 and \
                             row['similarity'] >= similarity_threshold
                if index < 1 and row['abs_mcs_size'] > 5 and (criterion1 or criterion2):
                    pdb_id = row['template_strucid'].split('.')[0].lower()
                    native_ligand = row['native_ligand'].molecule
                    docking_folder = Path(f'{pdb_id}_{srn}_{lig_cnt}')
                    scaffold = row['scaffold']
                    mcs = row['mcs_object']
                    scaffold_filename = docking_folder / Path(f'scaffold_{srn}.sdf')
                    strict_scaffold = \
                        mcs.return_mcs_scaffold(partial_ring_matches_allowed=False, ignore_hydrogens=False)[0]
                    _setup_docking(docking_folder, targets, pdb_id, project_home, ligand_mol, srn, scaffold,
                                   scaffold_filename, native_ligand, args, strict_scaffold)
                    _select_best_soln(docking_folder, scaffold, row, srn, pdb_id, ligand_entry.attributes,
                                      scaffold_filename)
                    dock = True
            if not dock:
                print('No suitable scaffold for ', srn)

        # except Exception as e:
        #     print(e)
        #     continue


def gold_rescoring(ligand_file: str, protein_file: str):
    print('rescoring......')

    with io.MoleculeReader(ligand_file) as rdr:
        ligand = rdr[0]
    ligand.add_hydrogens(mode='missing')
    target_protein = protein.Protein.from_file(protein_file)

    docker = Docker()
    settings = docker.settings
    settings.fitness_function = None
    settings.rescore_function = 'plp'
    settings.write_options = ['NO_LINK_FILES', 'NO_RNK_FILES', 'NO_PLP_MOL2_FILES', 'NO_BESTRANKING_LST_FILE',
                              'NO_GOLD_LIGAND_MOL2_FILE', 'NO_LOG_FILES', 'NO_FIT_PTS_FILES']
    settings._settings.set_rescore_with_simplex(False)
    settings._settings.set_fix_protein_rotatable_bonds(True)
    settings.add_protein_file(str(protein_file))

    settings.binding_site = settings.BindingSiteFromLigand(target_protein, ligand, 10.0)
    settings.add_ligand_file(ligand_file)
    settings.output_file = str(Path(f'./rescored_ligand.sdf'))

    gold_conf = 'gold_rescoring.conf'
    settings.make_absolute_file_names(gold_conf)
    settings.write(gold_conf)

    docker = Docker()
    docker.settings = docker.settings.from_file(gold_conf)
    docker.dock(file_name=gold_conf)

    ligand_attributes = {}
    if Path('rescored_ligand.sdf').is_file():
        ligand_entry = io.EntryReader('rescored_ligand.sdf')[0]
        ligand_attributes = ligand_entry.attributes
        ligand_attributes['Gold.PLP.Chemscore.Protein.Energy'] = 0
        ligand_entry.attributes = ligand_attributes
        with io.EntryWriter('rescored_ligand.sdf') as w:
            w.write(ligand_entry)

    return ligand_attributes


def main():
    args = parse_args()
    args.flexible_residues = [s.upper() for s in args.flexible_residues if s.strip()]
    project_home = '../..'

    print(project_home)
    if args.input_ligands:
        gold_scaffold_docking(args.input_ligands, project_home=project_home, args=args)
        print('finished scaffold docking')
    if args.docking_dir:
        gold_decoy_docking(args.docking_dir, args=args)
        print('finished decoy docking')


if __name__ == '__main__':
    main()
