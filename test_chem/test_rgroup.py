"""
Copyright (C) 2017 Anodyne Informatics, LLC
"""

import os
from typing import List
from unittest import TestCase
from unittest import main

from rdkit.Chem.rdchem import Mol

from rdkit import Chem
from ruse.rdkit.rgroup import Rgroup, RgroupDecomposer, mol_to_cores, Core
from rdkit.Chem.Descriptors import HeavyAtomMolWt

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


class TestRgroup(TestCase):

    def test_decomp(self) -> None:
        cores = ["*1CCCC1", "C1CCC1"]
        # cores = ["*1CCCC1"]
        smiles = ["C1CCCC1-C2CCC2Cl", "N1CCCC1-C2CCC2Cl", "O1CCCC1-C2CCC2Cl", "N1OCCC1-C2CCC2Cl", "N1OCSC1-C2CCC2Cl"]

        group = Rgroup()
        group.decompose(cores, smiles)

        self.assertTrue(len(group.columns), 4)

    def test_core_matching(self) -> None:
        self.core_matching('C1CC([*:2])C1[*:1]', 8, 2)

    def test_core_matching2(self) -> None:
        self.core_matching('*1C([*:1])C([*:2])CC1', 10, 4)

    def core_matching(self, core_smi: str, n_expected_substructure_matches, n_expected_core_matches) -> None:
        core = Core(0, Chem.MolFromSmarts(core_smi))

        smi = 'ClC1CCC1C1CCCC1'
        mol = Chem.MolFromSmiles(smi)
        (n_substructure_matches, n_core_matches) = (0, 0)

        for match in mol.GetSubstructMatches(core.matching_mol, uniquify=False):
            rest = Chem.ReplaceCore(mol, core.matching_mol, match, labelByIndex=True, replaceDummies=True)
            sidechains = Chem.GetMolFrags(rest, asMols=True, sanitizeFrags=False)
            sidechain_matches = core.match_sidechains(mol, match, sidechains)
            n_substructure_matches += 1
            if sidechain_matches:
                n_core_matches += 1

        self.assertEqual(n_core_matches, n_expected_core_matches)
        self.assertEqual(n_substructure_matches, n_expected_substructure_matches)

    def test_cyclic_match(self) -> None:
        core_smi = "[*:1]C1C([*:2])CCC1"
        mol_smi = 'C1(CCC2)C2CCCC1'

        core_mol = Chem.MolFromSmarts(core_smi)
        mol = Chem.MolFromSmiles(mol_smi)
        core = Core(0, core_mol)

        (n_substructure_matches, n_core_matches) = (0, 0)

        for match in mol.GetSubstructMatches(core.matching_mol, uniquify=False):
            rest = Chem.ReplaceCore(mol, core.matching_mol, match, labelByIndex=True, replaceDummies=True)
            sidechains = Chem.GetMolFrags(rest, asMols=True, sanitizeFrags=False)
            sidechain_matches = core.match_sidechains(mol, match, sidechains)
            n_substructure_matches += 1
            if sidechain_matches:
                n_core_matches += 1

        self.assertEqual(n_core_matches, 2)
        self.assertEqual(n_substructure_matches, 10)

    def test_decomposer(self) -> None:
        cores = ["*1CCCC1", "C1CCC1"]
        smiles = ["C1CCCC1-C2CCC2Cl", "N1CCCC1-C2CCC2Cl", "O1CCCC1-C2CCC2Cl", "N1OCCC1-C2CCC2Cl", "N1OCSC1-C2CCC2Cl"]

        group = RgroupDecomposer()
        group.decompose(cores, smiles)

        no_molecules = len(smiles)
        self.assertEqual(no_molecules, len(group.molecules))
        self.assertEqual(no_molecules, len(group.decomposition))
        # no of molecules decomposed
        self.assertEqual(5, len([d for d in group.decomposition if d]))
        # total number of decompositions
        self.assertEqual(8, len([i for d in group.decomposition if d for i in d]))

    def test_az_decomposer(self) -> None:
        cores = self.read_query_mols()
        molecules = self.read_molecules()

        group = RgroupDecomposer()
        group.decompose_molecules(cores, molecules)
        mollist = group.to_molecule_grid(True)
        self.assertEqual(952, len(mollist))

        self.assertEqual(len(cores), len(group.cores))
        no_molecules = len(molecules)

        self.assertEqual(no_molecules, len(group.molecules))
        self.assertEqual(no_molecules, len(group.decomposition))
        # no of molecules decomposed
        self.assertEqual(952, len([d for d in group.decomposition if d]))
        # total number of decompositions
        self.assertEqual(952, len([i for d in group.decomposition if d for i in d]))

        # r  group counts
        self.assertEqual(952, len([m[2] for m in mollist if m[2]]))
        self.assertEqual(952, len([m[2] for m in mollist if m[2] and HeavyAtomMolWt(m[2])]))
        self.assertEqual(935, len([m[3] for m in mollist if m[3]]))
        self.assertEqual(582, len([m[3] for m in mollist if m[3] and HeavyAtomMolWt(m[3])]))
        self.assertEqual(952, len([m[4] for m in mollist if m[4]]))
        self.assertEqual(952, len([m[4] for m in mollist if m[4] and HeavyAtomMolWt(m[4])]))

    def test_az_decomposer_match_first_core_only(self) -> None:
        cores = self.read_query_mols()
        molecules = self.read_molecules()

        group = RgroupDecomposer()
        group.match_first_core_only = True
        group.decompose_molecules(cores, molecules)
        mollist = group.to_molecule_grid(True)
        self.assertEqual(952, len(mollist))

        self.assertEqual(len(cores), len(group.cores))
        no_molecules = len(molecules)

        self.assertEqual(no_molecules, len(group.molecules))
        self.assertEqual(no_molecules, len(group.decomposition))
        # no of molecules decomposed
        self.assertEqual(952, len([d for d in group.decomposition if d]))
        # total number of decompositions
        self.assertEqual(952, len([i for d in group.decomposition if d for i in d]))

        # r  group counts
        self.assertEqual(952, len([m[2] for m in mollist if m[2]]))
        self.assertEqual(952, len([m[2] for m in mollist if m[2] and HeavyAtomMolWt(m[2])]))
        self.assertEqual(935, len([m[3] for m in mollist if m[3]]))
        self.assertEqual(582, len([m[3] for m in mollist if m[3] and HeavyAtomMolWt(m[3])]))
        self.assertEqual(952, len([m[4] for m in mollist if m[4]]))
        self.assertEqual(952, len([m[4] for m in mollist if m[4] and HeavyAtomMolWt(m[4])]))

    def test_az_rdkit_decomp(self) -> None:
        """
        Test of RDkit's decomp code on Matt's example compounds
        :return:
        """
        cores = self.read_query_mols()
        molecules = self.read_molecules()

        group = Rgroup()
        group.decompose_molecules(cores, molecules)
        mollist = group.to_molecule_grid()
        no_molecules = len(molecules)
        self.assertEqual(no_molecules, len(mollist))
        self.assertEqual(len(cores), len(group.cores))
        self.assertEqual(no_molecules, len(group.molecules))

        # r  group counts
        self.assertEqual(913, len([m[1] for m in mollist if m[1]]))
        self.assertEqual(913, len([m[2] for m in mollist if m[2]]))
        self.assertEqual(913, len([m[2] for m in mollist if m[2] and HeavyAtomMolWt(m[2])]))
        self.assertEqual(913, len([m[3] for m in mollist if m[3]]))
        self.assertEqual(544, len([m[3] for m in mollist if m[3] and HeavyAtomMolWt(m[3])]))
        self.assertEqual(913, len([m[4] for m in mollist if m[4]]))
        self.assertEqual(913, len([m[4] for m in mollist if m[4] and HeavyAtomMolWt(m[4])]))

    def test_structure1(self) -> None:
        molecules = [Chem.MolFromSmiles('COC(=O)c1cccc(N2CC(NC(=O)c3[nH]c(C)c(Cl)c3Cl)C2=O)c1')]
        cores = self.read_query_mols()[1:2]
        group = RgroupDecomposer()
        group.decompose_molecules(cores, molecules)
        mollist = group.to_molecule_grid(True)

        self.assertEqual(len(mollist), 1)
        self.assertEqual(len([m for m in mollist[0] if m]), 5)

    def test_symmetric_decomp(self) -> None:
        core_smarts = ["*1C([*:1])C([*:2])CC1", "C([*:1])1C([*:2])CC1"]
        smiles = ["C1CCCC1-C2CCC2Cl", "N1CCCC1-C2CCC2Cl", "O1CCCC1-C2CCC2Cl", "N1OCCC1-C2CCC2Cl", "N1OCSC1-C2CCC2Cl"]
        cores = [Chem.MolFromSmarts(s) for s in core_smarts]
        mols = [Chem.MolFromSmiles(s) for s in smiles]

        group = RgroupDecomposer()
        group.decompose_molecules(cores, mols, process_symmetric_groups=True)
        mol_lists = group.to_molecule_grid(True)

    @classmethod
    def read_query_mols(cls, file='pyrrolamide_cores.mol') -> List[Mol]:
        file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources", file))
        with open(file, 'r') as fh:
            mol_block = fh.read()
        return mol_to_cores(mol_block)

    @classmethod
    def read_molecules(cls, file='AZ_Pyrrolamides_4OE.csv') -> List[Mol]:
        if file.endswith('.mol') or file.endswith('sdf'):
            return cls.read_molecules_from_sdf(file)
        file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources", file))
        with open(file, 'r') as fh:
            lines = fh.readlines()
        return [Chem.MolFromSmiles(line.split(',')[1]) for line in lines[1:]]

    @classmethod
    def read_molecules_from_sdf(cls) -> List[Mol]:
        file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources", file))
        supplier = Chem.SDMolSupplier(file)
        return [m for m in supplier]


if __name__ == "__main__":
    main()
