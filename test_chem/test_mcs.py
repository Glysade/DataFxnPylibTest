import unittest
from unittest import TestCase, main

from rdkit import Chem

from ruse.rdkit.mcs import MCS, mcs_with_constant_smiles, mcs_match_fragment_smiles
from ruse.rdkit.mmp import create_combined_structure


class TestMcs(TestCase):

    @unittest.skip("Test is currently intractable")
    def test1(self):
        query = 'CCn1nnc(C2OC(n3cnc4c(N)nc(NC(CO)[CH2:1][c:1]5cc[c:2]([OH:2])cc5)nc43)C(O)C2O)n1'
        product = 'CCn1nnc(C2OC(n3cnc4c(N)nc(NC(CO)[CH2:1][c:1]5ccc(O)[c:2]([OH:2])c5)nc43)C(O)C2O)n1'
        mol_a = Chem.MolFromSmiles(query)
        mol_b = Chem.MolFromSmiles(product)
        mcs = MCS(mol_a, mol_b)
        mcs.find_mcs()

    def test2(self):
        query = 'c1ccccc1CCO'
        product = 'c1ccccc1CO'
        mol_a = Chem.MolFromSmiles(query)
        mol_b = Chem.MolFromSmiles(product)
        mcs = MCS(mol_a, mol_b, check_degree=False)
        mappings = mcs.find_mcs()
        self.assertEqual(len(mappings), 26)
        mcs = MCS(mol_a, mol_b, check_degree=True)
        mappings = mcs.find_mcs()
        self.assertEqual(len(mappings), 6)

    def test3(self):
        query = 'CCn1nnc(C2OC(n3cnc4c(N)nc(NC(CO)[CH2:1][c:1]5cc[c:2]([OH:2])cc5)nc43)C(O)C2O)n1'
        product = 'CCn1nnc(C2OC(n3cnc4c(N)nc(NC(CO)[CH2:1][c:1]5ccc(O)[c:2]([OH:2])c5)nc43)C(O)C2O)n1'
        constant = '[*]CC(CO)Nc1nc(N)c2ncn(C3OC(c4nnn(CC)n4)C(O)C3O)c2n1.[*]O'

        mol_a = Chem.MolFromSmiles(query)
        mol_b = Chem.MolFromSmiles(product)
        mapping = mcs_with_constant_smiles(mol_a, mol_b, constant, check_degree=False)
        a_mapping, b_mapping = list(zip(*mapping))

        a_missing = [a for a in range(mol_a.GetNumAtoms()) if a not in a_mapping]
        b_missing = [b for b in range(mol_b.GetNumAtoms()) if b not in b_mapping]

        self.assertEqual(len(a_missing), 1)
        self.assertEqual(len(b_missing), 2)

    def test4(self):
        query = 'CCn1nnc(C2OC(n3cnc4c(N)nc(N[CH:1](CO)[CH2:12][c:2]5ccc(O)cc5)nc43)C(O)C2O)n1'
        product = 'CCn1nnc(C2OC(n3cnc4c(N)nc(N[CH:1](CO)[CH2:1][CH2:2][c:2]5ccc(O)cc5)nc43)C(O)C2O)n1'
        constant = '[*]C(CO)Nc1nc(N)c2ncn(C3OC(c4nnn(CC)n4)C(O)C3O)c2n1.[*]c1ccc(O)cc1'

        mol_a = Chem.MolFromSmiles(query)
        mol_b = Chem.MolFromSmiles(product)
        mapping = mcs_with_constant_smiles(mol_a, mol_b, constant)
        a_mapping, b_mapping = list(zip(*mapping))

        a_missing = [a for a in range(mol_a.GetNumAtoms()) if a not in a_mapping]
        b_missing = [b for b in range(mol_b.GetNumAtoms()) if b not in b_mapping]

        self.assertEqual(len(a_missing), 1)
        self.assertEqual(len(b_missing), 2)

    def test5(self):
        query = 'Nc1c(N(C(=O)CCc2ccccc2Cl)[CH2:3][CH:123]([CH3:1])[CH3:2])c(=O)[nH]c(=O)n1Cc1ccccc1'
        product = 'C=[C:3](C[N:2](C[CH2:1][CH3:1])[CH3:2])[CH2:3]N(C(=O)CCc1ccccc1Cl)c1c(N)n(Cc2ccccc2)c(=O)[nH]c1=O'
        constant = '*C.*C.*CN(C(=O)CCc1ccccc1Cl)c1c(N)n(Cc2ccccc2)c(=O)[nH]c1=O'

        mol_a = Chem.MolFromSmiles(query)
        mol_b = Chem.MolFromSmiles(product)
        mapping = mcs_with_constant_smiles(mol_a, mol_b, constant)
        a_mapping, b_mapping = list(zip(*mapping))

        a_missing = [a for a in range(mol_a.GetNumAtoms()) if a not in a_mapping]
        b_missing = [b for b in range(mol_b.GetNumAtoms()) if b not in b_mapping]

        self.assertEqual(len(a_missing), 1)
        self.assertEqual(len(b_missing), 6)



class TestFragmentMcs(TestCase):

    def test_overlap(self):
        from_smiles = 'c1cc([*:2])ccc1[*:1]'
        to_smiles = 'c1cc([*:2])ccc1C[*:1]'
        atoms = mcs_match_fragment_smiles(from_smiles, to_smiles)
        self.assertEqual(len(atoms), 7)

    def test_overlap2(self):
        from_smiles = 'c1cc([*:2])ccc1[*:1]'
        to_smiles = 'Oc1cc([*:1])ccc1[*:2]'
        atoms = mcs_match_fragment_smiles(from_smiles, to_smiles, check_degree=False)
        self.assertEqual(len(atoms), 8)
        atoms = mcs_match_fragment_smiles(from_smiles, to_smiles, check_degree=True)
        self.assertEqual(len(atoms), 2)

    def test_overlap3(self):
        from_smiles = 'C([*:1])[*:2]'
        to_smiles = 'C(C[*:2])[*:1]'
        atoms = mcs_match_fragment_smiles(from_smiles, to_smiles)
        self.assertEqual(len(atoms), 2)

    def test_overlap4(self):
        from_smiles = 'c1cc([*:2])ccc1[*:1]'
        to_smiles = 'Oc1ccc([*:1])cc1[*:2]'
        atoms = mcs_match_fragment_smiles(from_smiles, to_smiles, check_degree=False)
        self.assertEqual(len(atoms), 7)
        atoms = mcs_match_fragment_smiles(from_smiles, to_smiles, check_degree=True)
        self.assertEqual(len(atoms), 2)

    def test_overlap5(self):
        from_smiles = 'Oc1ccc(C[*:1])cc1'
        to_smiles = 'O=[N+]([O-])c1ccc(C[*:1])cc1'
        atoms = mcs_match_fragment_smiles(from_smiles, to_smiles)
        self.assertEqual(len(atoms), 8)

    def test_overlap6(self):
        from_smiles = 'c1cc([*:2])ccc1C[*:1]'
        to_smiles = 'c1cc([*:2])ccc1[*:1]'
        atoms = mcs_match_fragment_smiles(from_smiles, to_smiles)
        self.assertEqual(len(atoms), 7)

    def test_fragment_mapping(self) -> None:
        from_smiles = 'c1cc([*:2])ccc1C[*:1]'
        to_smiles = 'c1cc([*:2])ccc1[*:1]'
        constant_smiles = '*C(CO)Nc1nc(N)c2ncn(C3OC(c4nnn(CC)n4)C(O)C3O)c2n1.*O'
        mapped_query_mol = create_combined_structure(constant_smiles, from_smiles)
        fragment_mapping = mcs_match_fragment_smiles(from_smiles, to_smiles)
        query_mapping, _ = list(zip(*fragment_mapping))
        query_missing = [atom.GetIdx() for atom in mapped_query_mol.GetAtoms()
                         if atom.HasProp('variableNo') and atom.GetIntProp('variableNo') not in query_mapping]
        self.assertEqual(len(query_mapping), 7)
        self.assertEqual(len(query_missing), 1)

    def test_fragment_mapping2(self) -> None:
        from_smiles = 'c1cc([*:2])ccc1C[*:1]'
        to_smiles = '[*:1]c1ccc([*:2])cc1'
        constant_smiles = '*C(CO)Nc1nc(N)c2ncn(C3OC(c4nnn(CC)n4)C(O)C3O)c2n1.*O'
        mapped_query_mol = create_combined_structure(constant_smiles, from_smiles)
        fragment_mapping = mcs_match_fragment_smiles(from_smiles, to_smiles)
        query_mapping, _ = list(zip(*fragment_mapping))
        query_missing = [atom.GetIdx() for atom in mapped_query_mol.GetAtoms()
                         if atom.HasProp('variableNo') and atom.GetIntProp('variableNo') not in query_mapping]
        self.assertEqual(len(query_mapping), 7)
        self.assertEqual(len(query_missing), 1)

    def test_correspondence_graph(self) -> None:
        from_smiles = 'c1cc([*:2])ccc1C[*:1]'
        to_smiles = '[*:1]c1ccc([*:2])cc1'
        mol_a = Chem.MolFromSmarts(from_smiles)
        mol_b = Chem.MolFromSmarts(to_smiles)
        mcs = MCS(mol_a, mol_b)
        graph = mcs.build_correspondence_graph()
        mapping = [(3, 5), (6, 1), (5, 7), (2, 4), (0, 2), (4, 6), (1, 3)]
        clique = mcs.mapping_to_clique(graph, mapping)
        try:
            mcs.check_clique(graph, clique)
        except RuntimeError:
            self.fail("Clique check failed")
        mapping.append((8, 0))
        clique = mcs.mapping_to_clique(graph, mapping)
        self.assertRaises(RuntimeError, mcs.check_clique, graph, mapping)


if __name__ == "__main__":
    main()
