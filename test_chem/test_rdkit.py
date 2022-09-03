"""
Copyright (C) 2017 Anodyne Informatics, LLC
"""

import gzip
import json
import os
import rdkit
from unittest import TestCase
from unittest import main

import psutil
from rdkit import Chem
from rdkit.Chem import AllChem

from ruse.rdkit.rdkit_utils import is_three_dimensional, remove_explicit_hydrogens
from ruse.chem.chem_data_table_helper import data_table_column_to_mols
from ruse.util.data_table import DataTable


class TestRDKit(TestCase):
    def test_supplier(self) -> None:
        """Can't use ForwardSDMolSupplier on a regular file handle"""
        print(f'RDKit version {rdkit.__version__}')
        file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources", 'test.sdf'))
        suppl = Chem.SDMolSupplier(file)
        for mol in suppl:
            self.assertGreater(mol.GetNumAtoms(), 5)
        # SDMolSupplier on a regular filed failed to close!  As of 2022.09.1pre it closes in the command line, but stays
        # open in PyCharm debug
        self.assertTrue(self.file_open(file))

        file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources", 'test2.sdf.gz'))
        fh = gzip.open(file, 'r')
        suppl = Chem.ForwardSDMolSupplier(fh)
        for mol in suppl:
            self.assertGreater(mol.GetNumAtoms(), 5)
        fh.close()
        self.assertFalse(self.file_open(file))

        file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources", 'test.sdf'))
        # opening forward supplier on a text file bombs
        fh = open(file, 'r')
        with self.assertRaises(ValueError):
            suppl = Chem.ForwardSDMolSupplier(fh)
        fh.close()
        self.assertFalse(self.file_open(file))

        # works ok if we have binary input
        file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources", 'test.sdf'))
        fh = open(file, 'rb')
        suppl = Chem.ForwardSDMolSupplier(fh)
        for mol in suppl:
            mol.GetNumAtoms()
        fh.close()
        self.assertFalse(self.file_open(file))

    def test_writer(self) -> None:
        file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources", 'test2.sdf.gz'))
        fh = gzip.open(file, 'r')
        suppl = Chem.ForwardSDMolSupplier(fh)
        out_file = 'test_out.sdf'
        writer = Chem.SDWriter(out_file)
        for mol in suppl:
            self.assertGreater(mol.GetNumAtoms(), 5)
            writer.write(mol)
        fh.close()

        writer.flush()
        writer.close()

        self.assertFalse(self.file_open(file))
        self.assertFalse(self.file_open(out_file))

        os.remove(out_file)

    @classmethod
    def file_open(cls, file) -> bool:
        process = psutil.Process(os.getpid())
        #open_files = [f.path for f in process.open_files()]
        #print("looking for file {} process has open paths {}".format(file, ', '.join(open_files)))
        test = [f.path for f in process.open_files() if f.path.endswith(file)]
        return len(test) > 0

    def test_dimension(self) -> None:
        file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources", 'test2.sdf.gz'));
        self.three_dimensional(file)

    def test_dimension2(self) -> None:
        file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources", '100conformers.sdf.gz'));
        self.three_dimensional(file)

    def testDimension3(self):
        mol = Chem.MolFromSmiles('c1ccccc1')
        self.assertTrue(mol is not None)
        self.assertFalse(is_three_dimensional(mol))

    def three_dimensional(self, file):
        fh = gzip.open(file, 'r')
        suppl = Chem.ForwardSDMolSupplier(fh)
        mol = next(suppl)
        self.assertTrue(mol is not None)
        self.assertGreater(mol.GetNumAtoms(), 5)
        self.assertTrue(mol.HasProp('_MolFileInfo'))
        dimension = mol.GetProp('_MolFileInfo')[20:22]
        self.assertEqual(dimension, '3D')
        self.assertTrue(is_three_dimensional(mol))
        fh.close()

    def test_smiles_supplier(self):
        file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources", '100smiles.smi'));
        supplier = Chem.SmilesMolSupplier(file, titleLine=False)
        mols = [m for m in supplier]
        self.assertEqual(len(mols), 100)

    def test_remove_h(self):
        smi1 = 'Nc1nc(NC(CO)Cc2ccc(O)cc2)nc2c1ncn2C1OC(c2nn[nH]n2)C(O)C1O'
        smi2 = 'Nc1nc(NC(CO)Cc2ccc(O)cc2)nc2c1ncn2C1OC(c2nnnn2)C(O)C1O'

        # AllChem.RemoveHs does not remove explicit H in smi1
        mol = Chem.MolFromSmiles(smi1)
        smi3 = Chem.MolToSmiles(AllChem.RemoveHs(mol, updateExplicitCount=True), True)
        self.assertNotEqual(smi2, smi3)

        # but remove_explicit_hydrogens does
        smi4 = Chem.MolToSmiles(remove_explicit_hydrogens(mol), True)
        self.assertEqual(smi2, smi4)


if __name__ == "__main__":
    main()
