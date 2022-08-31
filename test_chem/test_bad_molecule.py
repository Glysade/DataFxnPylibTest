import base64
import gzip
import os
from typing import List
from unittest import TestCase
from unittest import main

from ruse.chem.chem_data_table_helper import data_table_column_to_mols
from ruse.rdkit.rdkit_utils import RDKitFormat, string_to_mol
from ruse.util.data_table import DataTable


class TestBadMolecule(TestCase):

    def test_mol_to_string(self) -> None:
        entries = self.load_sample_sdf()
        self.assertEqual(2, len(entries))
        mol1 = string_to_mol(RDKitFormat.sdf, entries[0])
        mol2 = string_to_mol(RDKitFormat.sdf, entries[1])
        self.assertIsNotNone(mol1)
        self.assertIsNone(mol2)

    def test_data_table_column_to_mols(self):
        data_table = self.to_data_table();
        mols1 = list(data_table_column_to_mols(data_table, 0));
        self.assertEqual(1, len(mols1))
        self.assertIsNotNone(mols1[0])
        mols2 = list(data_table_column_to_mols(data_table, 0, True))
        self.assertEqual(2, len(mols2))
        self.assertIsNotNone(mols2[0])
        self.assertIsNone(mols2[1])

    def test_data_table_column_to_mols2(self):
        data_table = self.to_data_table(True)
        mols1 = list(data_table_column_to_mols(data_table, 0))
        self.assertEqual(1, len(mols1))
        self.assertIsNotNone(mols1[0])
        mols2 = list(data_table_column_to_mols(data_table, 0, True))
        self.assertEqual(2, len(mols2))
        self.assertIsNone(mols2[0])
        self.assertIsNotNone(mols2[1])

    def test_data_table_with_5_rows(self):
        data_table = self.to_large_data_table()
        mols1 = list(data_table_column_to_mols(data_table, 0))
        self.assertEqual(4, len(mols1))
        for mol in mols1:
            self.assertIsNotNone(mol);
        mols2 = list(data_table_column_to_mols(data_table, 0, True))
        self.assertEqual(5, len(mols2))
        for i, mol in enumerate(mols2):
            if i == 2:
                self.assertIsNone(mol)
            else:
                self.assertIsNotNone(mol)

    def load_sample_sdf(self) -> List[str]:
        file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources", 'ruse-test-sample.sdf'))
        with open(file, 'r') as fh:
            text = fh.read()
        entries = text.split('$$$$\n')
        return [e + '$$$$\n' for e in entries if e != '']

    def encode_sample_sdf(self, reverse=False, entries=None):
        if not entries:
            entries = self.load_sample_sdf()
        if reverse:
            entries.reverse();
        return [self.encode(e) for e in entries]

    def encode(self, sdf):
        mol_bytes = bytes(sdf, 'utf-8')
        mol_gzip = gzip.compress(mol_bytes)
        mol_b64 = base64.b64encode(mol_gzip).decode('utf-8')
        return mol_b64

    def to_data_table(self, reverse=False):
        data_table = DataTable()
        data_table.data = [[e] for e in self.encode_sample_sdf(reverse)]
        data_table.columns = [
            DataTable.column_definition("Structure", "binary", 'chemical/x-sdf')]
        return data_table

    def to_large_data_table(self):
        data_table = DataTable()
        a, b = self.load_sample_sdf()
        data_table.data = [[self.encode(e)] for e in (a, a, b, a, a)]
        data_table.columns = [
            DataTable.column_definition("Structure", "binary", 'chemical/x-sdf')]
        return data_table


if __name__ == "__main__":
    main()
