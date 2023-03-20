#!/usr/bin/env python

# Test script for test_linkers.py.

import csv
import unittest
from pathlib import Path
import df.find_linkers as fl

from rdkit import Chem


def linker_from_dict(linker_dict) -> fl.Linker:
    # only the linker_smiles needs to be in the dict,
    # the rest is optional for testing purposes.
    linker_smiles = linker_dict['linker_smiles']
    linker = Chem.MolFromSmiles(linker_smiles)
    linker_length = len(fl.shortest_path_between_dummies(linker))

    try:
        name = linker_dict['name']
    except KeyError:
        name = 'Str_0'

    try:
        num_atoms = linker_dict['num_atoms']
    except KeyError:
        num_atoms = 0

    try:
        mol_smiles = linker_dict['mol_smiles']
    except KeyError:
        mol_smiles = ''

    try:
        left_smiles = linker_dict['left_smiles']
    except KeyError:
        left_smiles = ''

    try:
        right_smiles = linker_dict['right_smiles']
    except KeyError:
        right_smiles = ''

    try:
        num_donors = linker_dict['num_donors']
    except KeyError:
        num_donors = 0

    try:
        num_acceptors = linker_dict['num_donors']
    except KeyError:
        num_acceptors = 0

    return fl.Linker(name=name, linker_smiles=linker_smiles,
                     num_atoms=num_atoms,
                     mol_smiles=mol_smiles, left_smiles=left_smiles,
                     right_smiles=right_smiles,
                     linker_length=linker_length, num_donors=num_donors,
                     num_acceptors=num_acceptors)


class TestLinkers(unittest.TestCase):
    def test_shortest_path(self) -> None:
        test_cases = [('*CCC*', [0, 1, 2, 3, 4]),
                      ('*C1CCCC(*)C1', [0, 1, 7, 5, 6])]
        for tc in test_cases:
            linker = Chem.MolFromSmiles(tc[0])
            sp = fl.shortest_path_between_dummies(linker)
            self.assertEqual(tc[1], sp)

    def test_linker_subs_ok(self) -> None:
        test_cases = [('*CCC*', True), ('*CC(C)C*', True),
                      ('*CC(CC)*', False), ('*CCCCC*', False), ('*CS(=O)(=O)*', True),
                      ('*C(OC)*', False), ('*C1COC(*)CC1', True),
                      ('*C1COC(*)C(CC)C1', False), ('*C1OC(*)C(CC)C1', False),
                      ('*C(C1CC1)*', False), ('*C(C)(CC)C*', False),
                      ('*C(C)(C)C*', True), ('*C1OC(C)(C)CC1*', True),
                      ('*C(C1C(F)C1)*', False), ('*CCCCCCCCCC*', False)]
        for tc in test_cases[0:]:
            linker = Chem.MolFromSmiles(tc[0])
            self.assertEqual(tc[1], fl.linker_is_ok(linker, 8), tc[0])

    def test_linker_bonds_ok(self) -> None:
        test_cases = [('*CCC*', True), ('*CCCC*', False),
                      ('*CC=C*', True), ('*C1CCCC(*)C1', True)]
        for tc in test_cases:
            linker = Chem.MolFromSmiles(tc[0])
            self.assertEqual(tc[1], fl.linker_is_ok(linker, 8), tc[0])

    def test_linker_atoms_ok(self) -> None:
        test_cases = [('*CCC*', True), ('*CCCC*', False),
                      ('*CNCC*', True)]
        for tc in test_cases:
            linker = Chem.MolFromSmiles(tc[0])
            self.assertEqual(tc[1], fl.linker_is_ok(linker, 8), tc[0])

    def test_make_subs(self) -> None:
        test_cases = [['c1ccccc1CCc1cccnc1',
                       [{'name': 'Str_0',
                         'linker_smiles': 'C(C[*:2])[*:1]',
                         'left_smiles': 'c1ccc([*:1])cc1',
                         'right_smiles': 'c1cncc([*:2])c1',
                         'mol_smiles': 'c1ccc(CCc2cccnc2)cc1',
                         'symmetrical': True}]],
                      ['c12ccccc2cncc1', []],
                      ['c1ccccc1', []],
                      ['c1ccccc1C2CCC(OC2)c1cccnc1', [{'name': 'Str_3',
                                                       'linker_smiles': 'C1CC([*:2])OCC1[*:1]',
                                                       'left_smiles': 'c1ccc([*:1])cc1',
                                                       'right_smiles': 'c1cncc([*:2])c1',
                                                       'mol_smiles': 'c1ccc(C2CCC(c3cccnc3)OC2)cc1',
                                                       'symmetrical': False},
                                                      {'name': 'Str_3',
                                                       'linker_smiles': 'C1CC([*:2])COC1[*:1]',
                                                       'left_smiles': 'c1cncc([*:1])c1',
                                                       'right_smiles': 'c1ccc([*:2])cc1',
                                                       'mol_smiles': 'c1ccc(C2CCC(c3cccnc3)OC2)cc1',
                                                       'symmetrical': False}]],
                      ['c1ccccc1C2CC(CCCCCCCCCC)C(OC2)c1cccnc1', []],
                      ['C12COCCC2CC3CCCNC3C1', []],
                      ['c1ccccc1CCc2ccccc2OCOc3cccnc3', [{'name': 'Str_6',
                                                          'linker_smiles': 'C(C[*:2])[*:1]',
                                                          'left_smiles': 'c1ccc([*:1])cc1',
                                                          'right_smiles': 'c1cncc(OCOc2ccccc2[*:2])c1',
                                                          'mol_smiles': 'c1ccc(CCc2ccccc2OCOc2cccnc2)cc1',
                                                          'symmetrical': True},
                                                         {'name': 'Str_6',
                                                          'linker_smiles': 'C(O[*:1])O[*:2]',
                                                          'left_smiles': 'c1ccc(CCc2ccccc2[*:1])cc1',
                                                          'right_smiles': 'c1cncc([*:2])c1',
                                                          'mol_smiles': 'c1ccc(CCc2ccccc2OCOc2cccnc2)cc1',
                                                          'symmetrical': True}]],
                      ['c1ccccc1CCc2[nH]ccc2CC(CCCC)c2ncncn2', [{'name': 'Str_7',
                                                                 'linker_smiles': 'C(C[*:2])[*:1]',
                                                                 'left_smiles': 'c1ccc([*:1])cc1',
                                                                 'right_smiles': 'CCCCC(Cc1cc[nH]c1[*:2])c1ncncn1',
                                                                 'mol_smiles': 'CCCCC(Cc1cc[nH]c1CCc1ccccc1)c1ncncn1',
                                                                 'symmetrical': True}]],
                      ['c1ccccc1CCCCc2[nH]ccc2CC=CCc2ncncn2', [{'name': 'Str_8',
                                                                'linker_smiles': 'C(=CC[*:2])C[*:1]',
                                                                'left_smiles': 'c1ccc(CCCCc2[nH]ccc2[*:1])cc1',
                                                                'right_smiles': 'c1ncnc([*:2])n1',
                                                                'mol_smiles': 'C(=CCc1cc[nH]c1CCCCc1ccccc1)Cc1ncncn1',
                                                                'symmetrical': True}]],
                      ['c1ccccc1CCCCc2[nH]ccc2C=CCCc2ncncn2', [{'name': 'Str_9',
                                                                'linker_smiles': 'C(=C[*:1])CC[*:2]',
                                                                'left_smiles': 'c1ccc(CCCCc2[nH]ccc2[*:1])cc1',
                                                                'right_smiles': 'c1ncnc([*:2])n1',
                                                                'mol_smiles': 'C(=Cc1cc[nH]c1CCCCc1ccccc1)CCc1ncncn1',
                                                                'symmetrical': False},
                                                               {'name': 'Str_9',
                                                                'linker_smiles': 'C(=C[*:2])CC[*:1]',
                                                                'left_smiles': 'c1ncnc([*:1])n1',
                                                                'right_smiles': 'c1ccc(CCCCc2[nH]ccc2[*:2])cc1',
                                                                'mol_smiles': 'C(=Cc1cc[nH]c1CCCCc1ccccc1)CCc1ncncn1',
                                                                'symmetrical': False}]]]
        for i, tc in enumerate(test_cases[0:]):
            _, linkers = fl.find_linkers((tc[0], f'Str_{i}'))
            test_linkers = [linker_from_dict(ld) for ld in tc[1]]
            self.assertEqual(len(test_linkers), len(linkers), f'Test {i}')
            self.assertEqual(test_linkers, linkers, f'Test {i}')

    def test_different_phenyls(self) -> None:
        # There was a time when these different environments for the
        # pbenyl linker gave different SMILES strings due to the
        # dummies being left with aromatic state.
        smis = ['C1CCC(CC1)c1ccc(cc1)C1CCCCC1', 'c1ncccc1c1ccc(cc1)C1CCCCC1',
                'c1ncccc1c1ccc(cc1)c1cccnc1']
        for i, smi in enumerate(smis):
            _, linkers = fl.find_linkers((smi, f'Str_{i + 1}'))
            self.assertEqual('c1cc([*:2])ccc1[*:1]',
                             linkers[0].linker_smiles, f'Test {i}')

    def test_split_linker1(self) -> None:
        # This one threw a Kekulization error at one point.  There
        # should be 2 linkers as it is asymmetrical.
        _, linkers = fl.find_linkers(('Cc1ccccc1-c1nn(-c2c(Cl)cc(Cl)cc2Cl)c(=N)s1', 'Str_1'))
        self.assertEqual(2, len(linkers))
        self.assertEqual('N=c1sc([*:1])nn1[*:2]', linkers[0].linker_smiles)

    def test_split_linker2(self) -> None:
        # There should only be 1 linker, the one from the two end rings
        # is too long.
        _, linkers = fl.find_linkers(('Cc1ccccc1CCc1nn(-c2c(Cl)cc(Cl)cc2Cl)c(=N)s1', 'Str_1'))
        self.assertEqual(1, len(linkers))

    def test_remove_atom_map_num(self) -> None:
        mol = Chem.MolFromSmiles('[CH3:1][C:2](=[O:3])[OH:4]')
        fl.remove_atom_maps(mol)
        self.assertEqual('CC(=O)O', Chem.MolToSmiles(mol))

    def test_is_linker_symmetric(self) -> None:
        test_cases = [('O([*:1])[*:2]', True), ('C([*:1])[*:2]', True),
                      ('O=C(N[*:1])[*:2]', False), ('c1cc([*:1])cc([*:2])c1', True),
                      ('c1cc([*:1])nc([*:2])n1', False)]
        for tc in test_cases:
            linker = fl.Linker(name='Str_0', linker_smiles=tc[0],
                               num_atoms=0, mol_smiles='',
                               left_smiles='', right_smiles='',
                               linker_length=0, num_donors=0, num_acceptors=0)
            self.assertEqual(tc[1], linker.symmetrical)

    def test_linkers_match(self) -> None:
        test_cases = [({'linker_smiles': 'O([*:1])[*:2]',
                        'num_atoms': 0,
                        'left_smiles': 'c1cncc([*:1])c1',
                        'right_smiles': 'c1ccc([*:2])cc1'},
                       {'linker_smiles': 'O([*:1])[*:2]',
                        'num_atoms': 0,
                        'left_smiles': 'c1cncc([*:1])c1',
                        'right_smiles': 'c1ccc([*:2])cc1'},
                       True),
                      ({'linker_smiles': 'C([*:1])[*:2]',
                        'num_atoms': 0,
                        'left_smiles': 'c1cncc([*:1])c1',
                        'right_smiles': 'c1ccc([*:2])cc1'},
                       {'linker_smiles': 'O([*:1])[*:2]',
                        'num_atoms': 0,
                        'left_smiles': 'c1cncc([*:1])c1',
                        'right_smiles': 'c1ccc([*:2])cc1'},
                       False),
                      ({'linker_smiles': 'O([*:1])[*:2]',
                        'num_atoms': 0,
                        'left_smiles': 'c1cncc([*:1])c1',
                        'right_smiles': 'c1ccc([*:2])cc1'},
                       {'linker_smiles': 'O([*:1])[*:2]',
                        'num_atoms': 0,
                        'left_smiles': 'c1ccc([*:1])cc1',
                        'right_smiles': 'c1cncc([*:2])c1'},
                       True),
                      ({'linker_smiles': 'O([*:1])[*:2]',
                        'num_atoms': 0,
                        'left_smiles': 'c1ccc([*:1])cc1',
                        'right_smiles': 'O=C1C2CCCCC2C(=O)N1CC(=O)N1CCCSC1=[*:2]'},
                       {'linker_smiles': 'O([*:1])[*:2]',
                        'num_atoms': 0,
                        'left_smiles': 'O=C(CN1C(=O)C2CCCCC2C1=O)N1CCCSC1=[*:1]',
                        'right_smiles': 'c1ccc([*:2])cc1'},
                       True),
                      ]
        for i, tc in enumerate(test_cases):
            l1 = linker_from_dict(tc[0])
            l2 = linker_from_dict(tc[1])
            self.assertEqual(tc[2], l1 == l2, f'Test {i}')

    def test_write_csv(self) -> None:
        smis = ['c1ccccc1CCc1cccnc1', 'c1ccccc1CCCCc2[nH]ccc2C=CCCc2ncncn2',
                'c1ccccc1C2CCC(OC2)c1cccnc1', 'c1ccccc1CCc1ccccc1']
        mols = [(smi, f'Str_{i + 1}') for i, smi in enumerate(smis)]
        all_linkers = fl.find_all_linkers(mols, 8, 5)
        sorted_linkers = fl.sort_linkers(all_linkers)

        temp_path = Path('test_write_csv.csv')
        res = fl.write_csv_linkers(sorted_linkers, temp_path)
        self.assertTrue(res)
        with open(temp_path, 'r', newline='') as fi:
            csvreader = csv.reader(fi)
            rows = [row for row in csvreader]
        self.assertEqual(6, len(rows))
        self.assertEqual(3, len(rows[0]))
        self.assertEqual(3, len(rows[-1]))
        self.assertEqual('2', rows[1][2])
        self.assertEqual('Linker_5', rows[-1][0])
        self.assertEqual('C(=C[*:2])CC[*:1]', rows[3][1])

        # temp_path.unlink()

    def test_write_sdf(self) -> None:
        smis = ['c1ccccc1CCc1cccnc1', 'c1ccccc1CCCCc2[nH]ccc2C=CCCc2ncncn2',
                'c1ccccc1C2CCC(OC2)c1cccnc1', 'c1ccccc1CCc1ccccc1']
        mols = [(smi, f'Str_{i + 1}') for i, smi in enumerate(smis)]
        all_linkers = fl.find_all_linkers(mols, 8, 5)
        sorted_linkers = fl.sort_linkers(all_linkers)

        temp_path = Path('test_write_sdf.sdf')
        res = fl.write_sdf_linkers(sorted_linkers, temp_path)
        self.assertTrue(res)

        temp_path.unlink()

    def test_write_sqlite(self) -> None:
        smis = ['c1ccccc1CCc1cccnc1', 'c1ccccc1CCCCc2[nH]ccc2C=CCCc2ncncn2',
                'c1ccccc1C2CCC(OC2)c1cccnc1', 'c1ccccc1CCc1ccccc1']
        mols = [(smi, f'Str_{i + 1}') for i, smi in enumerate(smis)]
        all_linkers = fl.find_all_linkers(mols, 8, 5)
        sorted_linkers = fl.sort_linkers(all_linkers)

        temp_path = Path('test_write_sqlite.db')
        res = fl.write_sqlite_linkers(sorted_linkers, temp_path)
        self.assertTrue(res)

        temp_path.unlink()

    def test_linker_obj(self) -> None:
        l1 = fl.Linker('L1', 'O([*:1])[*:2]', 'c1ccncc1Oc1ccccc1',
                       'c1cncc([*:1])c1', 'c1ccc([*:2])cc1', 3, 3, 0, 1)
        self.assertTrue(l1.symmetrical)
        self.assertEqual(3, l1.path_length)
        self.assertEqual(3, l1.num_atoms)
        l2 = fl.Linker('L2', 'O([*:1])[*:2]', 'c1ccncc1Oc1ccccc1',
                       'c1cncc([*:1])c1', 'c1ccc([*:2])cc1', 3, 3, 0, 1)
        self.assertEqual(l1, l2)
        l3 = fl.Linker('L3', 'O([*:1])[*:2]', 'c1ccccc1Oc1cnccc1',
                       'c1ccc([*:1])cc1', 'c1cncc([*:2])c1', 3, 3, 0, 1)
        self.assertEqual(l1, l3)

    def test_find_linkers(self) -> None:
        smi = 'c1ccccc1Oc1cnccc1'
        mol_name, linkers = fl.find_linkers((smi, 'Str_0'))
        self.assertEqual(mol_name, 'Str_0')
        exp_linker = fl.Linker(name='Str_0', mol_smiles=smi,
                               linker_smiles='O([*:1])[*:2]',
                               num_atoms=3,
                               left_smiles='c1ccc([*:1])cc1',
                               right_smiles='c1cncc([*:2])c1',
                               linker_length=3,
                               num_donors=0,
                               num_acceptors=1)
        self.assertEqual(exp_linker, linkers[0])

    def test_hbond_counts(self) -> None:
        test_cases = [('O([*:1])[*:2]', 0, 1), ('[*:1]C(O)C[*:2]', 1, 1),
                      ('[*:1]C[*:2]', 0, 0), ('[*:1]NC(=O)[*:2]', 1, 2),
                      ('[*:1]c1c[nH]cc1[*:2]', 1, 1), ('O=C(N[*:1])N[*:2]', 2, 3)]
        for tc in test_cases:
            num_don, num_acc = fl.count_donors_and_acceptors(Chem.MolFromSmiles(tc[0]))
            self.assertEqual(tc[1], num_don, f'{tc[0]} : donor')
            self.assertEqual(tc[2], num_acc, f'{tc[0]} : acceptor')


if __name__ == '__main__':
    unittest.main()
