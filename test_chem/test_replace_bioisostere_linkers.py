import unittest

from os import environ
from pathlib import Path
from typing import Optional

from rdkit import Chem, rdBase

import df.replace_bioisostere_linkers as rbl


def find_replacements_file() -> Optional[Path]:
    ppath_dirs = environ['PYTHONPATH'].split(';')
    for ppd in ppath_dirs:
        ppd_par = Path(ppd).parent
        data_dir = ppd_par / 'Data'
        if data_dir.exists():
            repl_file = data_dir / 'chembl_32_bioisostere_linkers.db'
            return repl_file
    return None


def count_linker_atoms(mol: Chem.Mol) -> int:
    num_lats = 0
    for atom in mol.GetAtoms():
        try:
            _ = atom.GetIntProp('Linker')
            num_lats += 1
        except KeyError:
            pass
    return num_lats


class TestReplaceBioisostereLinkers(unittest.TestCase):

    def setUp(self) -> None:
        print(f'Using RDKit version {rdBase.rdkitVersion}')

    def test_db_files(self):
        db_file = str(Path(__file__).parent / 'resources' / 'no_such_database.db')
        self.assertRaises(FileNotFoundError, rbl.check_db_file, db_file)
        dbs = ['test_db_no_bios.db', 'test_db_no_linkers.db']
        for bad_db in dbs:
            bad_db_file = str(Path(__file__).parent / 'resources' / bad_db)
            self.assertRaises(ValueError, rbl.check_db_file, bad_db_file)

    def test_fetch_bioisosteres1(self) -> None:
        db_file = str(Path(__file__).parent / 'resources' / 'chembl_31_bios.db')
        bios = rbl.fetch_bioisosteres('C(C[*:4])[*:3]', db_file, -1, -1,
                                      False, False, False)
        self.assertEqual(36, len(bios))
        self.assertEqual('C(#C[*:4])[*:3]', bios[0])
        self.assertEqual('c1cc([*:4])ccc1[*:3]', bios[-1])
        bios = rbl.fetch_bioisosteres('C(C[*:4])[*:3]', db_file, -1, -1,
                                      False, False, True)
        self.assertEqual(35, len(bios))
        self.assertEqual('C(#C[*:4])[*:3]', bios[0])
        self.assertEqual('O[C@H](C[*:4])[*:3]', bios[-1])

    def test_replace0(self) -> None:
        db_file = str(Path(__file__).parent / 'resources' / 'chembl_31_bios.db')
        mol = Chem.MolFromSmiles('c12ccccc2cc[nH]1')
        new_mols, _, _ = rbl.replace_linkers(mol, db_file, 8, 5, -1, -1,
                                             False, False, False, -1)
        # should send back an empty list because molecule has no
        # linkers
        self.assertEqual(0, len(new_mols))

    def test_replace1(self) -> None:
        db_file = str(Path(__file__).parent / 'resources' / 'chembl_31_bios.db')
        # This bit of database has 22 replacements for *CC*.
        mol = Chem.MolFromSmiles('c1ccccc1CCc1cnccc1')
        new_mols, _, linker_smis = \
            rbl.replace_linkers(mol, db_file, 8, 5, -1, -1, False, False,
                                False, 10000)
        self.assertEqual(36, len(new_mols))
        self.assertEqual('C(#Cc1cccnc1)c1ccccc1', Chem.MolToSmiles(new_mols[0]))
        self.assertEqual(36, len(linker_smis))
        # Currently the linker atoms in the query copy aren't being tagged in
        # the output molecule.
        # num_lats = count_linker_atoms(query_cp)
        # self.assertEqual(2, num_lats)

    def test_replace1a(self) -> None:
        db_file = str(Path(__file__).parent / 'resources' / 'chembl_31_bios.db')
        mol = Chem.MolFromSmiles('c1cccc1CCC')
        new_mols, _, _ = rbl.replace_linkers(mol, db_file, 8, 5, -1, -1,
                                             False, False, False, -1)
        self.assertEqual(0, len(new_mols))

    def test_replace2(self) -> None:
        db_file = str(Path(__file__).parent / 'resources' / 'chembl_31_bios.db')
        mol = Chem.MolFromSmiles('c1ccccc1CCc1cnccc1OCOc1ccccc1')
        new_mols, _, linker_smis = \
            rbl.replace_linkers(mol, db_file, 8, 5, -1, -1, False, False,
                                False, 10000)
        # It's the same number as test_replace1 because this snippet of
        # db doesn't have a replacement for *OCO*.
        self.assertEqual(36, len(new_mols))
        self.assertEqual('C(#Cc1cnccc1OCOc1ccccc1)c1ccccc1',
                         Chem.MolToSmiles(new_mols[0]))
        self.assertEqual(36, len(linker_smis))
        # num_lats = count_linker_atoms(query_cp)
        # self.assertEqual(5, num_lats)

    def test_replace3(self) -> None:
        db_file = str(Path(__file__).parent / 'resources' / 'chembl_31_bios.db')
        mol = Chem.MolFromSmiles('c1ccccc1CCc1cnccc1OCOc1cc(SC2CCNCC2)ccc1')
        new_mols, _, _ = rbl.replace_linkers(mol, db_file, 8, 5, -1, -1,
                                             False, False, False, -1)
        # Expect 332 mols, 37 * 1 * 9 - 8 replacements of *S* plus
        # itself and minus the original input molecule.
        self.assertEqual(332, len(new_mols))
        self.assertEqual('c1ccc(CCc2cnccc2OCOc2cccc(OCC3CCNCC3)c2)cc1',
                         Chem.MolToSmiles(new_mols[0]))
        # num_lats = count_linker_atoms(query_cp)
        # self.assertEqual(6, num_lats)

    def test_bulk_replace1(self) -> None:
        db_file = str(Path(__file__).parent / 'resources' / 'chembl_31_bios.db')
        infile = str(Path(__file__).parent / 'resources' / 'test_bulk_replacement.smi')
        new_mols, query_cps, _ = rbl.bulk_replace_linkers(infile, db_file, 8, 5, -1, -1,
                                                          False, False, False, -1, 3)
        self.assertEqual(84, len(new_mols))
        self.assertEqual(84, len(query_cps))
        # make sure they're in the input order
        last_smi = None
        smis = []
        for qc in query_cps:
            smi = Chem.MolToSmiles(qc)
            if smi != last_smi:
                smis.append(smi)
                last_smi = smi
        # the 6th mol has no linkers, so produces no results.
        self.assertEqual(6, len(smis))
        self.assertEqual(smis, ['CCc1cc(-c2ccccc2)nnc1NCCN1CCOCC1',
                                'CCc1cc(-c2ccccc2)nnc1CCCN1CCOCC1',
                                'CCN1CCCC1Nc1cc(C)c(-c2ccccc2O)nn1',
                                'CCN1CCCC1Cc1cc(C)c(-c2ccccc2O)nn1',
                                'Cc1cc(-c2cccs2)nnc1NCCN1CCOCC1',
                                'CCc1cc(-c2ccccc2)nnc1C(=O)N1CCOCC1'])

    def test_replace_with_limits(self) -> None:
        db_file = str(Path(__file__).parent / 'resources' / 'chembl_31_bios.db')
        mol = Chem.MolFromSmiles('c1c[nH]cc1CC1CCNCC1')
        new_mols, _, _ = rbl.replace_linkers(mol, db_file, 8, 5,
                                             -1, -1, False, False, False, -1)
        self.assertEqual(25, len(new_mols))
        self.assertEqual('C(=C/C1CCNCC1)\c1cc[nH]c1', Chem.MolToSmiles(new_mols[1]))
        mol = Chem.MolFromSmiles('c1c[nH]cc1CC1CCNCC1')
        new_mols, _, _ = rbl.replace_linkers(mol, db_file, 8, 5,
                                             1, 1, False, False, False, -1)
        self.assertEqual(18, len(new_mols))
        self.assertEqual('C(=C/C1CCNCC1)\c1cc[nH]c1', Chem.MolToSmiles(new_mols[1]))

    def test_replace_with_hbond_matches(self) -> None:
        db_file = str(Path(__file__).parent / 'resources' / 'chembl_31_bios.db')
        mol = Chem.MolFromSmiles('c1c[nH]cc1OC1CCNCC1')
        new_mols, _, _ = rbl.replace_linkers(mol, db_file, 8, 5,
                                             -1, -1, False, False, False, -1)
        self.assertEqual(20, len(new_mols))
        self.assertEqual('C(=C/C1CCNCC1)\c1cc[nH]c1', Chem.MolToSmiles(new_mols[1]))

        mol = Chem.MolFromSmiles('c1c[nH]cc1OC1CCNCC1')
        new_mols, _, _ = rbl.replace_linkers(mol, db_file, 8, 5,
                                             -1, -1, True, True, False, -1)
        self.assertEqual(12, len(new_mols))
        self.assertEqual('c1cc(CCCOC2CCNCC2)c[nH]1', Chem.MolToSmiles(new_mols[1]))

        mol = Chem.MolFromSmiles('c1c[nH]cc1C(=O)NC1CCNCC1')
        new_mols, _, _ = rbl.replace_linkers(mol, db_file, 8, 5,
                                             -1, -1, True, True, False, -1)
        self.assertEqual(28, len(new_mols))
        self.assertEqual('c1cc(NCCC2CCNCC2)c[nH]1', Chem.MolToSmiles(new_mols[1]))

        mol = Chem.MolFromSmiles('c1c[nH]cc1C(O)C1CCNCC1')
        new_mols, _, _ = rbl.replace_linkers(mol, db_file, 8, 5,
                                             -1, -1, True, False, False, -1)
        self.assertEqual(1, len(new_mols))
        self.assertEqual('CC(O)(c1cc[nH]c1)C1CCNCC1', Chem.MolToSmiles(new_mols[0]))

    def test_no_ring_linkers(self) -> None:
        repl_file = find_replacements_file()
        self.assertIsNotNone(repl_file, "Couldn't find replacements file")
        mol = Chem.MolFromSmiles('c1c(F)cccc1c1ccccc1c1ncccc1')
        new_mols, _, _ = rbl.replace_linkers(mol, repl_file, 8, 5,
                                             -1, -1, True, False, False, -1)
        self.assertEqual(14, len(new_mols))

        new_mols, _, _ = rbl.replace_linkers(mol, repl_file, 8, 5,
                                             -1, -1, True, False, True, -1)
        self.assertEqual(0, len(new_mols))

        mol = Chem.MolFromSmiles('c1c(F)cccc1CCCc1ncccc1')
        new_mols, _, _ = rbl.replace_linkers(mol, repl_file, 8, 5,
                                             -1, -1, True, False, False, -1)
        self.assertEqual(63, len(new_mols))

        new_mols, _, _ = rbl.replace_linkers(mol, repl_file, 8, 5,
                                             -1, -1, True, False, True, -1)
        self.assertEqual(61, len(new_mols))

    def test_max_output_mols(self) -> None:
        repl_file = find_replacements_file()
        self.assertIsNotNone(repl_file, "Couldn't find replacements file")
        smi = 'Cc1ccc(C(=O)N2CCCC(c3cccc(F)c3)C2)cc1'
        mol = Chem.MolFromSmiles(smi)
        new_mols, _, _ = rbl.replace_linkers(mol, repl_file, 8, 5,
                                             1, 1, False, False, False, 100)
        self.assertEqual(50, len(new_mols))

        # The more complicated original.
        smi = 'Cc1ccc(C(=O)N2CCCC(c3n[nH]cc3-c3cccc(F)c3)C2)cc1NCc1cccnc1'
        mol = Chem.MolFromSmiles(smi)
        self.assertIsNotNone(repl_file, "Couldn't find replacements file")
        new_mols, _, _ = rbl.replace_linkers(mol, repl_file, 8, 5,
                                             1, 1, False, False, False, 100)
        self.assertEqual(100, len(new_mols))

    def test_bad_mol_1(self) -> None:
        # This one failed before because one of the possible linkers,
        # C1C(*)CCCN1C(=O)*, can also produce the linker *C(=O)*
        # because it's between the piperidine and a phenyl.  The fix
        # was not to allow the larger linker in these circumstances.
        repl_file = find_replacements_file()
        self.assertIsNotNone(repl_file, "Couldn't find replacements file")
        mol = Chem.MolFromSmiles('Cc1ccc(C(=O)N2CCCC(c3cccc(F)c3)C2)cc1')
        new_mols, _, _ = rbl.replace_linkers(mol, repl_file, 8, 5,
                                             1, 1, False, False, False, 100)
        self.assertEqual(50, len(new_mols))

        # The more complicated original.
        smi = 'Fc1cccc(c1)c1c[nH]nc1C1CN(C(=O)c2cc(C)c(NCc3cnccc3)cc2)CCC1'
        mol = Chem.MolFromSmiles(smi)
        self.assertIsNotNone(repl_file, "Couldn't find replacements file")
        new_mols, _, l1s = rbl.replace_linkers(mol, repl_file, 8, 5,
                                               1, 1, False, False, False, -1)
        self.assertEqual(47735, len(new_mols))

        # A simple positional isomer that once gave a different number
        # of analogues.
        smi = 'Cc1ccc(C(=O)N2CCCC(c3n[nH]cc3-c3cccc(F)c3)C2)cc1NCc1cccnc1'
        mol = Chem.MolFromSmiles(smi)
        new_mols, _, l2s = rbl.replace_linkers(mol, repl_file, 8, 5,
                                               1, 1, False, False, False, -1)
        self.assertEqual(47735, len(new_mols))

    def test_bad_mol_2(self) -> None:
        # These have 3-way linkers, which isn't allowed.
        # They come out in the wash, though, because each cut to make
        # the linker produces a substituent that is too big.
        repl_file = 'resources/chembl_31_bios.db'
        smi = 'c1ccncc1c1cc(c2c[nH]cc2)cc(c2cocc2)c1'
        mol = Chem.MolFromSmiles(smi)
        new_mols, _, _ = rbl.replace_linkers(mol, repl_file, 8, 5,
                                             1, 1, False, False, False, 100)
        self.assertEqual(0, len(new_mols))

        smi = 'c1[nH]ccc1C(c1cocc1)c1cccnc1'
        mol = Chem.MolFromSmiles(smi)
        new_mols, _, _ = rbl.replace_linkers(mol, repl_file, 8, 5,
                                             1, 1, False, False, False, 100)
        self.assertEqual(0, len(new_mols))


if __name__ == '__main__':
    unittest.main()
