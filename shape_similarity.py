#!/usr/bin/env python

import argparse
import openeye
import time

from openeye.oechem import OESetSDData, OEGetSDDataPairs, OEGraphMol, OESuppressHydrogens, OEAddExplicitHydrogens, \
    OEWriteMolecule, OEReadMolecule
# from openeye.oeomega import *
from openeye.oeshape import OEColorForceField, OEColorFFType_ExplicitMillsDean, OEColorOptions, OEOverlapPrep, \
    OEOverlapResults, OERemoveColorAtoms, OEExactColorFunc, OEOverlapFunc


def parse_args():
    """Define and parse the arguments to the script."""
    parser = argparse.ArgumentParser(
        description='Execute Line of sight contact scripts.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter  # To display default values in help message.
    )

    parser.add_argument(
        '--template_ligand',
        help='Path to template_ligand.',
        default=''
    )

    parser.add_argument(
        '--docked_ligand',
        help='Path to docking_ligand.',
    )

    return parser.parse_args()


class Rescorer(object):

    def __init__(self):
        cff = OEColorForceField()
        cff.Clear()
        # cff.Init(OEColorFFType_ImplicitMillsDean)
        cff.Init(OEColorFFType_ExplicitMillsDean)
        self.prep = OEOverlapPrep()
        self.prep.SetColorForceField(cff)
        self.prep.SetUseHydrogens(True)

        options = OEColorOptions()
        options.SetColorForceField(cff)
        self.func = OEOverlapFunc(OEExactColorFunc(options))

    def is_similar(self, m, _diverse_confs):
        """
        Make sure that conformation is sufficiently different from those already saved
        :param m:
        :param _diverse_confs:
        :return:
        """
        res = OEOverlapResults()
        conf_is_similar = False
        for diverse_conf in _diverse_confs:
            self.func.SetupRef(diverse_conf)
            self.func.Overlap(m, res)
            print('TanimotCombo', res.GetTanimotoCombo())
            if res.GetTanimotoCombo() > 1.9:
                conf_is_similar = True
        return conf_is_similar

    def rescore_and_write(self, templates, dmols, nmax, ofs):
        """
        molecules are rescored with ROCS to make sure the non-matched parts overlap as well as possible
        :param templates:
        :param dmols:
        :param nmax:
        :param ofs:
        :return:
        """

        # colorFunc = OEAnalyticColorFunc()
        t1 = time.time()
        print('start template scoring')
        res = OEOverlapResults()
        res_val = []

        for i, dmol in enumerate(dmols):
            dmol.SweepConfs()
            self.prep.Prep(dmol)

            confs = list(dmol.GetConfs())
            for t in templates:
                # --> templates already prepared with "prep.Prep(t) "
                self.func.SetupRef(t)
                for k, conf in enumerate(confs):
                    self.func.Overlap(conf, res)
                    res_val.append([i, t.GetTitle(), k, res.GetTanimotoCombo(), conf])
                    # print("SHAPE %f COLOR %f COMBO %f" %(res.GetTanimoto(), res.GetColorTanimoto(), res.GetTanimotoCombo())) # REMOVE

        res_val = sorted(res_val, key=lambda x: x[3], reverse=True)
        t2 = time.time()
        print(f'finish: time {t2 - t1}')
        diverse_confs = []
        for j, lx in enumerate(res_val):
            if j >= nmax:
                break
            if (lx[0], lx[2]) in diverse_confs:
                continue
            else:
                diverse_confs.append((lx[0], lx[2]))
            tmpmol = OEGraphMol(lx[4])

            # re-generate hydrogen positions as they can be distorted in rare occasions
            OESuppressHydrogens(tmpmol)
            OEAddExplicitHydrogens(tmpmol)

            # --> add all SD-tags from first molecule (are all the same anyway)
            for dp in OEGetSDDataPairs(dmols[0]):
                OESetSDData(tmpmol, dp.GetTag(), dp.GetValue())
            OESetSDData(tmpmol, "TanimotoCombo", str(lx[3]))
            OESetSDData(tmpmol, "Query", str(lx[1]))
            OESetSDData(tmpmol, "Alignment_Type", "MCS")

            # --> remove added pseudo atoms from ROCS calc.
            OERemoveColorAtoms(tmpmol)
            OEWriteMolecule(ofs, tmpmol)
        return True


def shape_similarity(template_ligand_path, docked_ligand_path):
    # shape similarity to template
    oe_native_ligand_stream = openeye.oechem.oemolistream(template_ligand_path)
    oe_native_ligand = openeye.oechem.OEGraphMol()
    openeye.oechem.OEReadMolecule(oe_native_ligand_stream, oe_native_ligand)

    oe_docked_ligand_stream = openeye.oechem.oemolistream(docked_ligand_path)
    max_tanimoto_combo = 0
    for oe_docked_ligand in oe_docked_ligand_stream.GetOEMols():
        # oe_docked_ligand = openeye.oechem.OEGraphMol()
        # openeye.oechem.OEReadMolecule(oe_docked_ligand_stream, oe_docked_ligand)

        rescorer = Rescorer()
        rescorer.prep.Prep(oe_native_ligand)
        rescorer.func.SetupRef(oe_native_ligand)

        res = openeye.oeshape.OEOverlapResults()
        rescorer.prep.Prep(oe_docked_ligand)
        rescorer.func.Overlap(oe_docked_ligand, res)
        TanimotoCombo = res.GetTanimotoCombo()
        if TanimotoCombo > max_tanimoto_combo:
            max_tanimoto_combo = TanimotoCombo
    print(max_tanimoto_combo)
    return max_tanimoto_combo


def main():
    args = parse_args()
    shape_similarity(args.template_ligand, args.docked_ligand)


if __name__ == '__main__':
    main()
