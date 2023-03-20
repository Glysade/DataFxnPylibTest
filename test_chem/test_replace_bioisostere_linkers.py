import unittest

from rdkit import Chem, rdBase

import df.replace_bioisostere_linkers as rbl


class TestReplaceBioisostereLinkers(unittest.TestCase):

    def setUp(self) -> None:
        print(f'Using RDKit version {rdBase.rdkitVersion}')

    def test_db_files(self):
        db_file = 'resources/no_such_database.db'
        self.assertRaises(FileNotFoundError, rbl.check_db_file, db_file)
        dbs = ['resources/test_db_no_bios.db', 'resources/test_db_no_linkers.db']
        for bad_db in dbs:
            self.assertRaises(ValueError, rbl.check_db_file, bad_db)

    def test_fetch_bioisosteres1(self) -> None:
        db_file = 'resources/chembl_31_bios.db'
        bios = rbl.fetch_bioisosteres('C(C[*:2])[*:1]', db_file, -1, -1,
                                      False, False)
        self.assertEqual(22, len(bios))
        self.assertEqual('C([*:1])[*:2]', bios[0])
        self.assertEqual('CC(C[*:1])[*:2]', bios[-1])

    def test_replace0(self) -> None:
        db_file = 'resources/chembl_31_bios.db'
        new_mols = rbl.replace_linkers('c12ccccc2cc[nH]1', db_file, 8, 5, -1, -1,
                                       False, False)
        self.assertEqual(1, len(new_mols))
        # Should come back unaltered.
        self.assertEqual('c1ccc2[nH]ccc2c1', Chem.MolToSmiles(new_mols[0]))

    def test_replace1(self) -> None:
        db_file = 'resources/chembl_31_bios.db'
        # This bit of database has 22 replacements for *CC*.
        new_mols = rbl.replace_linkers('c1ccccc1CCc1cnccc1', db_file, 8, 5, -1, -1,
                                       False, False)
        self.assertEqual(22, len(new_mols))
        self.assertEqual('c1ccc(Cc2cccnc2)cc1', Chem.MolToSmiles(new_mols[0]))

    def test_replace1a(self) -> None:
        db_file = 'resources/chembl_31_bios.db'
        # This bit of database has 22 replacements for *CC*.
        new_mols = rbl.replace_linkers('c1cccc1CCC', db_file, 8, 5, -1, -1,
                                       False, False)
        self.assertEqual(0, len(new_mols))

    def test_replace2(self) -> None:
        db_file = 'resources/chembl_31_bios.db'
        new_mols = rbl.replace_linkers('c1ccccc1CCc1cnccc1OCOc1ccccc1',
                                       db_file, 8, 5, -1, -1, False, False)
        # It's the same number as test_replace1 because this snippet of
        # db doesn't have a replacement for *OCO*.
        self.assertEqual(22, len(new_mols))
        self.assertEqual('c1ccc(Cc2cnccc2OCOc2ccccc2)cc1',
                         Chem.MolToSmiles(new_mols[0]))

    def test_replace3(self) -> None:
        db_file = 'resources/chembl_31_bios.db'
        new_mols = rbl.replace_linkers('c1ccccc1CCc1cnccc1OCOc1cc(SC2CCNCC2)ccc1',
                                       db_file, 8, 5, -1, -1, False, False)
        # Expect 132 mols, 22 * 1 * 6 - 6 replacements of *S*
        self.assertEqual(132, len(new_mols))
        self.assertEqual('c1ccc(Cc2cnccc2OCOc2cccc(CC3CCNCC3)c2)cc1',
                         Chem.MolToSmiles(new_mols[0]))

    def test_bulk_replace1(self) -> None:
        db_file = 'resources/chembl_31_bios.db'
        infile = 'resources/test_bulk_replacement.smi'
        new_mols = rbl.bulk_replace_linkers(infile, db_file, 8, 5, -1, -1,
                                            False, False)
        self.assertEqual(37, len(new_mols))

    def test_replace_with_limits(self) -> None:
        db_file = 'resources/chembl_31_bios.db'
        new_mols = rbl.replace_linkers('c1c[nH]cc1CC1CCNCC1', db_file, 8, 5,
                                       -1, -1, False, False)
        self.assertEqual(19, len(new_mols))
        self.assertEqual('c1cc(CCCC2CCNCC2)c[nH]1', Chem.MolToSmiles(new_mols[1]))
        new_mols = rbl.replace_linkers('c1c[nH]cc1CC1CCNCC1', db_file, 8, 5,
                                       1, 1, False, False)
        self.assertEqual(14, len(new_mols))
        self.assertEqual('c1cc(OCC2CCNCC2)c[nH]1', Chem.MolToSmiles(new_mols[1]))

    def test_replace_with_hbond_matches(self) -> None:
        db_file = 'resources/chembl_31_bios.db'
        new_mols = rbl.replace_linkers('c1c[nH]cc1OC1CCNCC1', db_file, 8, 5,
                                       -1, -1, False, False)
        self.assertEqual(14, len(new_mols))
        self.assertEqual('c1cc(CC2CCNCC2)c[nH]1', Chem.MolToSmiles(new_mols[1]))

        new_mols = rbl.replace_linkers('c1c[nH]cc1OC1CCNCC1', db_file, 8, 5,
                                       -1, -1, True, True)
        self.assertEqual(9, len(new_mols))
        self.assertEqual('O=C(NC1CCNCC1)c1cc[nH]c1', Chem.MolToSmiles(new_mols[1]))

        new_mols = rbl.replace_linkers('c1c[nH]cc1C(=O)NC1CCNCC1', db_file, 8, 5,
                                       -1, -1, True, True)
        self.assertEqual(6, len(new_mols))
        self.assertEqual('O=C(Nc1cc[nH]c1)C1CCNCC1', Chem.MolToSmiles(new_mols[1]))

        new_mols = rbl.replace_linkers('c1c[nH]cc1C(O)C1CCNCC1', db_file, 8, 5,
                                       -1, -1, True, False)
        self.assertEqual(1, len(new_mols))
        self.assertEqual('CC(O)(c1cc[nH]c1)C1CCNCC1', Chem.MolToSmiles(new_mols[0]))


if __name__ == '__main__':
    unittest.main()
