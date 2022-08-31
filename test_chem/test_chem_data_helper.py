import glob
import os
from unittest import TestCase, main

from ruse.chem.chem_data_table_helper import create_data_table_from_mols, data_table_column_to_local_files
from ruse.rdkit.rdkit_utils import mol_supplier


class ChemDataHelperTestCase(TestCase):

    def test_data_table_column_to_local_files(self):

        input_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources", "single_conformers.sdf.gz"))
        with mol_supplier(input_file) as supplier:
            mols = [mol for mol in supplier if mol]

        table = create_data_table_from_mols(mols, content_type='chemical/x-sdf')
        self.assertEqual(table.columns[0]['properties']['ContentType'], 'chemical/x-sdf')
        files = data_table_column_to_local_files(table, 0, 'export_mols', 1)

        self.assertEqual(len(files), len(mols))
        globbed_files = glob.glob('export_mols_*')
        self.assertSetEqual(set(files), set(globbed_files))

        for f in globbed_files:
            os.remove(f)


if __name__ == "__main__":
    main()