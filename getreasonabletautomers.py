#!/usr/bin/env python

# (C) 2017 OpenEye Scientific Software Inc. All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of OpenEye products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. OpenEye claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable OpenEye offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall OpenEye be
# liable for any damages or liability in connection with the Sample Code
# or its use.


import sys
from openeye import oechem
from openeye import oequacpac


def is_bad_tautomer(tautomer):
    '''
    https://docs.eyesopen.com/toolkits/python/oechemtk/patternmatch.html#section-patternmatch-sss
    :param tautomer:
    :return:
    '''

    bad_smarts = {'triazole_4H': 'c1[nD2][nD2]c[nH]1',
                  'pyridone1': 'C1C=CC(=O)N=C1',
                  'pyridone2': 'C1C=CC=NC1=O',
                  'iminol': 'C(-O)=!@N',
                  'enol': 'C(-O)=!@[#6]',
                  'imine1': '[#6]-[#7]=[#6]1-[*]=[*]-[*]-[*]=[*]1',
                  'imine2': '[#6]-[#7]=[#6]1[*]:[*]:[*]:[*][*]1'}
    for smarts in bad_smarts.values():
        # create a substructure search object
        ss = oechem.OESubSearch(smarts)
        oechem.OEPrepareSearch(tautomer, ss)
        if ss.SingleMatch(tautomer):
            print(smarts)
            return True
        else:
            continue
    return False


def extended_enumeration(mol):
    extended_smarts = {'pyridone': 'n@~[#6]=O',
                       'triazole': 'nnncc'}
    for smarts in extended_smarts.values():
        # create a substructure search object
        ss = oechem.OESubSearch(smarts)
        oechem.OEPrepareSearch(mol, ss)
        if ss.SingleMatch(mol):
            return True
    return False


def main(argv=[__name__]):
    if len(argv) != 3:
        oechem.OEThrow.Usage("%s <mol-infile> <mol-outfile>" % argv[0])

    ifs = oechem.oemolistream()
    if not ifs.open(argv[1]):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % argv[1])

    ofs = oechem.oemolostream()
    if not ofs.open(argv[2]):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % argv[2])

    tautomerOptions = oequacpac.OETautomerOptions()
    pKaNorm = True

    for mol in ifs.GetOEGraphMols():
        if extended_enumeration(mol):
            tautomerOptions.SetCarbonHybridization(False)
            tautomerOptions.SetSaveStereo(True)
            tautomers = oequacpac.OEEnumerateTautomers(mol, tautomerOptions)
        else:
            tautomers = oequacpac.OEGetReasonableTautomers(mol, tautomerOptions, pKaNorm)
        for tautomer in tautomers:
            if is_bad_tautomer(tautomer):
                continue
            oechem.OEAssignAromaticFlags(tautomer)
            oechem.OEWriteMolecule(ofs, tautomer)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
