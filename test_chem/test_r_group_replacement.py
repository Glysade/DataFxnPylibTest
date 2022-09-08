from pathlib import Path
from unittest import TestCase, main

from df.chem_helper import column_to_molecules
from df.data_transfer import DataFunction, DataFunctionRequest, DataFunctionResponse
from rdkit import Chem

from df.RGroupReplacement import RGroupReplacement, replace_rgroups


def run_script(in_file: str, test_class: DataFunction) -> DataFunctionResponse:
    with open(in_file, 'r') as fh:
        request_json = fh.read()

    request = DataFunctionRequest.parse_raw(request_json)
    response = test_class.execute(request)
    return response


class ScriptTest(TestCase):

    def test_script_smarts_core(self) -> None:
        file_in = Path(__file__).parent / 'resources' / 'test_r_group_replacement.json'
        rgr = RGroupReplacement()
        response = run_script(file_in, rgr)
        self.assertTrue(response)

    def test_script_smarts_core2(self) -> None:
        file_in = Path(__file__).parent / 'resources' / 'test_r_group_replacement1.json'
        rgr = RGroupReplacement()
        response = run_script(file_in, rgr)
        self.assertTrue(response)
        print(f'number of columns : {len(response.outputTables[0].columns)}')
        print(f'number of prods : {len(response.outputTables[0].columns[1].values)}')
        self.assertEqual(len(response.outputTables[0].columns), 2)
        self.assertEqual(len(response.outputTables[0].columns[1].values), 42)
        mols = column_to_molecules(response.outputTables[0].columns[1])
        self.assertEqual(Chem.MolToSmiles(mols[0]), 'N#Cc1cc(S(=O)(=O)Nc2ccco2)ccc1Oc1ccccc1-c1ccccc1')

    def test_script_molfile_core(self) -> None:
        file_in = Path(__file__).parent / 'resources' / 'test_r_group_replacement2.json'
        rgr = RGroupReplacement()
        response = run_script(file_in, rgr)
        self.assertTrue(response)
        print(f'number of columns : {len(response.outputTables[0].columns)}')
        print(f'number of prods : {len(response.outputTables[0].columns[1].values)}')
        self.assertEqual(len(response.outputTables[0].columns), 2)
        mols = column_to_molecules(response.outputTables[0].columns[1])
        self.assertEqual(len(response.outputTables[0].columns[1].values), 42)
        self.assertEqual(Chem.MolToSmiles(mols[0]), 'N#Cc1cc(S(=O)(=O)Nc2ccco2)ccc1Oc1ccccc1-c1ccccc1')


    def test_2_subs_on_atom(self) -> None:
        mols = [Chem.MolFromSmiles('Cl[C@H](OC)c1nc(O)c[nH]1')]
        core_query = Chem.MolFromSmarts('Cc1nccn1')
        analogue_table = replace_rgroups(mols, core_query, True, False,
                                         'Test')
        analogues = column_to_molecules(analogue_table.columns[1])
        self.assertEqual(len(analogues), 125)
        analogue_table = replace_rgroups(mols, core_query, False, True,
                                         'Test')
        analogues = column_to_molecules(analogue_table.columns[1])
        self.assertEqual(len(analogues), 504)

    def test_duplicate_removal(self) -> None:
        mols = [Chem.MolFromSmiles('Clc1ccccc1Cl')]
        core_query = Chem.MolFromSmarts('c1ccccc1')
        analogue_table = replace_rgroups(mols, core_query, True, False,
                                         'Test')
        analogues = column_to_molecules(analogue_table.columns[1])
        self.assertEqual(len(analogues), 15)

    def test_multiple_matches(self) -> None:
        mols = [Chem.MolFromSmiles('OC(=O)c1ccccc1Oc1ccc(N)cc1')]
        core_query = Chem.MolFromSmarts('c1ccccc1')
        analogue_table = replace_rgroups(mols, core_query, True, False,
                                         'Test')
        analogues = column_to_molecules(analogue_table.columns[1])
        self.assertEqual(len(analogues), 10)

    def test_dismiss_bidentate(self) -> None:
        mols = [Chem.MolFromSmiles('c1ccc2c(c1)COC2')]
        core_query = Chem.MolFromSmarts('c1ccccc1')
        analogue_table = replace_rgroups(mols, core_query, True, False,
                                         'Test')
        # we should just get the same molecule back again
        analogues = column_to_molecules(analogue_table.columns[1])
        self.assertEqual(len(analogues), 1)
        self.assertEqual(Chem.MolToSmiles(analogues[0]), 'c1ccc2c(c1)COC2')

    def test_do_nothing(self) -> None:
        smis = ['[nH]1c(C)nc(F)c1', '[nH]1c(C)nc(F)c1Cl', '[nH]1c(C)nc(O)c1',
                '[nH]1c(CO)nc(F)c1c1ccc(Cl)cc1']
        mols = [Chem.MolFromSmiles(s) for s in smis]
        core_query = Chem.MolFromSmarts('c1nccn1')
        analogue_table = replace_rgroups(mols, core_query, False, False,
                                         'Test')
        analogues = column_to_molecules(analogue_table.columns[1])
        self.assertEqual(len(analogues), 0)

if __name__ == '__main__':
    main()
