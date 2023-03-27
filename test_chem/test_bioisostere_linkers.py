#!/usr/bin/env python

import sqlite3
import unittest

from os import close
from pathlib import Path
from tempfile import mkstemp
from typing import Union

from rdkit import rdBase, Chem

import df.find_bioisotere_linkers as fbl


def create_mols(smis: list[str]) -> list[Chem.Mol]:
    molecules = []
    for smi in smis:
        mol = Chem.MolFromSmiles(smi)
        if mol and mol.GetNumAtoms():
            molecules.append(mol)

    return molecules


def create_series(smis: list[str]) -> dict[str, Union[str, list[str]]]:
    names = [f'CHEMBL{i}' for i in range(len(smis))]
    series = {"doc": "test_doc",
              "assay": "test_assay",
              "SMILES": smis,
              "ChEMBLID": names}
    return series


class TestBioisostereLinkers(unittest.TestCase):

    def setUp(self) -> None:
        print(f'Using RDKit version {rdBase.rdkitVersion}.')

    def test_series1(self) -> None:
        # Straightforward case.
        smis = ['c1ncccc1Cc1ccccc1', 'c1ncccc1Oc1ccccc1', 'c1ncccc1CCc1ccccc1']
        series = create_series(smis)
        bioisosteres = fbl.find_bioisosteres_in_series(series, 8, 5)
        self.assertEqual(3, len(bioisosteres))
        self.assertEqual('C([*:1])[*:2]', bioisosteres[0]._linker1)
        self.assertEqual('O([*:1])[*:2]', bioisosteres[0]._linker2)
        self.assertEqual('C(C[*:2])[*:1]', bioisosteres[1]._linker1)
        self.assertEqual('C([*:1])[*:2]', bioisosteres[1]._linker2)
        self.assertEqual('C(C[*:2])[*:1]', bioisosteres[2]._linker1)
        self.assertEqual('O([*:1])[*:2]', bioisosteres[2]._linker2)

    def test_series2(self) -> None:
        # Symmetrical case, molecules in different order.
        smis = ['c1ncccc1Oc1ccccc1', 'c1ccccc1Cc1cnccc1']
        series = create_series(smis)
        bioisosteres = fbl.find_bioisosteres_in_series(series, 8, 5)
        self.assertEqual(1, len(bioisosteres))
        self.assertEqual('C([*:1])[*:2]', bioisosteres[0]._linker1)
        self.assertEqual('O([*:1])[*:2]', bioisosteres[0]._linker2)

    def test_series3(self) -> None:
        # once upon a time this gave 2 bioisosteres, because it's
        # an asymmetric linker, but that doesn't make sense. There
        # should be 1 bioisostere with 1 example.
        smis = ["c1ccc2c(c1)CCCc1cc(NCCCN3CCOCC3)nnc1-2",
                "c1ccc2c(c1)CCCc1cc(NCCN3CCOCC3)nnc1-2"]
        series = create_series(smis)
        bioisosteres = fbl.find_bioisosteres_in_series(series, 8, 5)
        self.assertEqual(1, len(bioisosteres), 'wrong number of bioisosteres')
        self.assertEqual(1, len(bioisosteres[0]._examples), 'wrong number of examples')

    def test_series4(self) -> None:
        smis = [
            "CCc1cc(-c2ccccc2)nnc1NCCN1CCOCC1",
            "CCc1cc(-c2ccccc2)nnc1CCCN1CCOCC1",
            "CCN1CCCC1Nc1cc(C)c(-c2ccccc2O)nn1",
            "CCN1CCCC1Cc1cc(C)c(-c2ccccc2O)nn1",
            "Cc1cc(-c2cccs2)nnc1NCCN1CCOCC1",
            "c1ccc2c(c1)CCCc1cc(NCCCN3CCOCC3)nnc1-2",
            "Cc1cc(NCCN2CCOCC2)nnc1-c1ccccc1",
            "CCCc1cc(NC2CCCN2CC)nnc1-c1ccccc1",
            "CCN(CC)C(C)(C)CNc1cc(C2CC2)c(-c2ccccc2)nn1",
            "CCN(CC)C(C)CNc1cc(C2CC2)c(-c2ccccc2)nn1",
            "CCN1CCCC1CNc1nnc(-c2ccccc2)cc1C",
            "Cc1cc(NCCN(C)C)ccc1-c1ccccc1",
            "CCN(CC)CCNc1nnc(-c2ccccc2)cc1C",
            "CCN(CC)CCNc1cc2c(nn1)-c1ccccc1CCC2",
            "CCN(CC)CCCNc1cc2c(nn1)-c1ccccc1CCC2",
            "Cc1cc(-c2ccc3ccccc3c2)nnc1NCCN1CCOCC1",
            "CCCc1nnc(NCC(C)(C)N(CC)CC)cc1-c1ccccc1",
            "CCCc1cc(NC2CCCN2CC)nnc1-c1ccccc1",
            "CCN(CC)C(C)(C)CNc1cc(-c2ccccc2)c(CC(C)C)nn1",
            "c1ccc2c(c1)CCCc1cc(NCCN3CCOCC3)nnc1-2",
            "Cc1cc(-c2ccccc2)nnc1NCCN(C)C",
            "CCCc1cc(NCCN(CC)CC)nnc1-c1ccccc1",
            "CCN1CCCC1Nc1cc(-c2ccccc2)c(C)nn1",
            "c1ccc(Cc2cc(-c3ccccc3)nnc2NCCN2CCOCC2)cc1",
            "CCN(CC)C(C)(C)CNc1cc(-c2ccccc2)c(C)nn1",
            "CCCc1cc(NCCN2CCOCC2)nnc1-c1ccccc1",
            "CCN(CC)C(C)(C)CNc1cc(C)c(-c2ccccc2)nn1",
            "CCCc1cc(NCC(C)(C)N(CC)CC)nnc1-c1ccccc1O",
            "CCN(CC)C(C)(C)CNc1cc(-c2ccccc2)c(-c2ccccc2)nn1",
            "Cc1c(NCCN2CCOCC2)nnc(-c2ccccc2)c1C",
            "c1ccc(-c2cc(NCCN3CCOCC3)nnc2-c2ccccc2)cc1",
            "CCN(CC)C(C)(C)CNc1cc(C)c(-c2ccccc2O)nn1",
            "Cc1cc(NCCCN(C)C)ccc1-c1ccccc1",
            "CCCc1cc(NC2CCCN2CC)nnc1-c1ccccc1",
            "Cc1cc(C2CCCCC2)nnc1NCCN1CCOCC1",
            "CCN1CCCC1Nc1cc(C)c(-c2ccccc2)nn1",
            "Cc1cc(-c2cccs2)nnc1NCCN1CCOCC1",
            "c1ccc(-c2cc(-c3ccccc3)c(NCCN3CCOCC3)nn2)cc1",
            "CCCc1cc(NCC(C)(C)N(CC)CC)nnc1-c1ccccc1",
            "CCc1cc(-c2ccccc2)nnc1NCCN1CCNCC1",
            "CCc1cc(-c2ccccc2)nnc1CCCN1CCNCC1"
        ]
        series = create_series(smis)
        bioisosteres = fbl.find_bioisosteres_in_series(series, 8, 5)
        exp_res = [('C(C[*:1])C[*:2]', 'C(C[*:2])N[*:1]', 2),
                   ('C(CN[*:1])C[*:2]', 'C(C[*:2])N[*:1]', 1),
                   ('C([*:1])[*:2]', 'N([*:1])[*:2]', 1)]
        self.assertEqual(3, len(bioisosteres))
        for bios, er in zip(bioisosteres, exp_res):
            self.assertEqual(er[0], bios.linker1_smiles)
            self.assertEqual(er[1], bios.linker2_smiles)
            self.assertEqual(er[2], bios.num_examples)

    def test_full1(self) -> None:
        cli_args = ['--input-series-file',
                    'resources/find_bioisosteres_test1.json',
                    '--min-examples', '3',
                    '--output-db-file']
        fh, db_file = mkstemp(suffix='.db', dir=Path.cwd())
        close(fh)
        cli_args.append(db_file)

        res = fbl.main(cli_args)
        self.assertTrue(res)
        db_path = Path(db_file)
        self.assertTrue(db_path.exists())
        conn = sqlite3.connect(db_file)
        sql = 'SELECT * FROM bioisosteres'
        bioisostere = conn.execute(sql).fetchall()
        self.assertEqual(1, len(bioisostere))
        self.assertEqual(bioisostere[0][1], 'C(C[*:2])[*:1]')
        self.assertEqual(bioisostere[0][2], 'C([*:1])[*:2]')
        self.assertEqual(bioisostere[0][3], 3)
        self.assertEqual(bioisostere[0][4], 4)
        conn.close()
        db_path.unlink()


if __name__ == '__main__':
    unittest.main()
