from pathlib import Path
from unittest import TestCase, main

from df.chem_helper import column_to_molecules
from df.data_transfer import DataFunction, DataFunctionRequest, DataFunctionResponse,\
    DataType
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
        self.assertEqual(len(response.outputTables[0].columns), 3)
        self.assertEqual(len(response.outputTables[0].columns[1].values), 710)
        self.assertEqual(len(response.outputTables[0].columns[2].values), 710)
        self.assertEqual(len(response.outputColumns[0].values), 5)
        self.assertEqual(response.outputColumns[0].values[0], 7)
        self.assertEqual(response.outputColumns[0].values[1], 98)
        mols = column_to_molecules(response.outputTables[0].columns[2])

    def test_script_smarts_core2(self) -> None:
        file_in = Path(__file__).parent / 'resources' / 'test_r_group_replacement1.json'
        rgr = RGroupReplacement()
        response = run_script(file_in, rgr)
        self.assertTrue(response)
        self.assertEqual(len(response.outputTables[0].columns), 3)
        self.assertEqual(len(response.outputTables[0].columns[1].values), 42)
        self.assertEqual(len(response.outputTables[0].columns[2].values), 42)
        mols = column_to_molecules(response.outputTables[0].columns[2])
        self.assertEqual(Chem.MolToSmiles(mols[0]), 'N#Cc1cc(S(=O)(=O)Nc2ccco2)ccc1Oc1ccccc1-c1ccccc1')
        self.assertEqual(len(response.outputColumns[0].values), 100)
        self.assertEqual(response.outputColumns[0].values[0], 3)
        self.assertEqual(response.outputColumns[0].values[1], 0)

    def test_script_molfile_core(self) -> None:
        file_in = Path(__file__).parent / 'resources' / 'test_r_group_replacement2.json'
        rgr = RGroupReplacement()
        response = run_script(file_in, rgr)
        self.assertTrue(response)
        self.assertEqual(len(response.outputTables[0].columns), 3)
        mols = column_to_molecules(response.outputTables[0].columns[2])
        self.assertEqual(len(response.outputTables[0].columns[1].values), 42)
        self.assertEqual(len(response.outputTables[0].columns[2].values), 42)
        self.assertEqual(Chem.MolToSmiles(mols[0]), 'N#Cc1cc(S(=O)(=O)Nc2ccco2)ccc1Oc1ccccc1-c1ccccc1')
        self.assertEqual(len(response.outputColumns[0].values), 100)
        self.assertEqual(response.outputColumns[0].values[0], 3)
        self.assertEqual(response.outputColumns[0].values[1], 0)

    def test_2_subs_on_atom(self) -> None:
        mols = [Chem.MolFromSmiles('Cl[C@H](OC)c1nc(O)c[nH]1')]
        ids = ['mol1']
        core_query = Chem.MolFromSmarts('Cc1nccn1')
        analogue_table, _ = replace_rgroups(mols, ids, DataType.STRING, core_query,
                                            True, False, 'Test')
        analogues = column_to_molecules(analogue_table.columns[2])
        self.assertEqual(len(analogues), 125)
        analogue_table, _ = replace_rgroups(mols, ids, DataType.STRING, core_query,
                                            False, True, 'Test')
        analogues = column_to_molecules(analogue_table.columns[2])
        self.assertEqual(len(analogues), 504)

    def test_duplicate_removal(self) -> None:
        mols = [Chem.MolFromSmiles('Clc1ccccc1Cl')]
        ids = ['mol1']
        core_query = Chem.MolFromSmarts('c1ccccc1')
        analogue_table, _ = replace_rgroups(mols, ids, DataType.STRING, core_query,
                                            True, False, 'Test')
        analogues = column_to_molecules(analogue_table.columns[2])
        self.assertEqual(len(analogues), 15)

    def test_multiple_matches(self) -> None:
        mols = [Chem.MolFromSmiles('OC(=O)c1ccccc1Oc1ccc(N)cc1')]
        ids = ['mol1']
        core_query = Chem.MolFromSmarts('c1ccccc1')
        analogue_table, _ = replace_rgroups(mols, ids, DataType.STRING, core_query,
                                            True, False, 'Test')
        analogues = column_to_molecules(analogue_table.columns[2])
        self.assertEqual(len(analogues), 10)

    def test_dismiss_bidentate(self) -> None:
        mols = [Chem.MolFromSmiles('c1ccc2c(c1)COC2')]
        ids = ['mol1']
        core_query = Chem.MolFromSmarts('c1ccccc1')
        analogue_table, analogue_counts = replace_rgroups(mols, ids, DataType.STRING, core_query,
                                                          True, False, 'Test')
        # we should just get an empty table
        analogues = column_to_molecules(analogue_table.columns[2])
        self.assertEqual(len(analogues), 0)
        self.assertEqual(len(analogue_counts), 1)
        self.assertEqual(analogue_counts[0], 0)

    def test_do_nothing(self) -> None:
        smis = ['[nH]1c(C)nc(F)c1', '[nH]1c(C)nc(F)c1Cl', '[nH]1c(C)nc(O)c1',
                '[nH]1c(CO)nc(F)c1c1ccc(Cl)cc1']
        mols = [Chem.MolFromSmiles(s) for s in smis]
        ids = ['mol1', 'mol2', 'mol3', 'mol4']
        core_query = Chem.MolFromSmarts('c1nccn1')
        analogue_table, analogue_counts = replace_rgroups(mols, ids, DataType.STRING, core_query,
                                                          False, False, 'Test')
        analogues = column_to_molecules(analogue_table.columns[2])
        self.assertEqual(len(analogues), 0)
        self.assertEqual(len(analogue_counts), 4)
        self.assertEqual(len([ac for ac in analogue_counts if ac == 0]), 4)

    def test_script_bad_smarts_core(self) -> None:
        # a poor core in this case, which results in r groups with
        # partial aromatic rings that used to throw an exception.
        file_in = Path(__file__).parent / 'resources' / 'test_r_group_replacement3.json'
        rgr = RGroupReplacement()
        response = run_script(file_in, rgr)
        self.assertTrue(response)
        self.assertEqual(len(response.outputTables[0].columns), 3)
        self.assertEqual(len(response.outputTables[0].columns[1].values), 176)
        self.assertEqual(len(response.outputTables[0].columns[2].values), 176)
        self.assertEqual(len(response.outputColumns[0].values), 15)
        self.assertEqual(response.outputColumns[0].values[0], 0)
        self.assertEqual(response.outputColumns[0].values[3], 75)

    def test_analogue_in_input_set(self) -> None:
        smis = ['COc1ccccc1N', 'CCOc1ccccc1N']
        mols = [Chem.MolFromSmiles(s) for s in smis]
        ids = ['mol1', 'mol2']
        core_query = Chem.MolFromSmarts('Nc1ccccc1')
        analogue_table, _ = replace_rgroups(mols, ids, DataType.STRING, core_query,
                                            True, False, 'Test')
        analogues = column_to_molecules(analogue_table.columns[2])
        an_smi = [Chem.MolToSmiles(an) for an in analogues]
        self.assertNotIn(smis[1], an_smi)

    def test_dont_colour_no_change_r_groups(self) -> None:
        # because there's no replacement for FCCCCCCCCCCCCO, it should
        # not be coloured.
        # Also check layer 1 and layer 2 colours are correct.
        smis = ['FCCCCCCCCCCCCOc1ccc(OCC)cc1N']
        mols = [Chem.MolFromSmiles(s) for s in smis]
        ids = ['mol1']
        core_query = Chem.MolFromSmarts('Nc1ccccc1')
        analogue_table, _ = replace_rgroups(mols, ids, DataType.STRING, core_query,
                                            True, True, 'Test')
        analogues = column_to_molecules(analogue_table.columns[2])
        render_string1 = 'COLOR #0000ff\nATOMS 1 2 3 4 5 6 7\nBONDS 1 2 3 4 7 5 6\nCOLOR #ff0000\nATOMS 22 23\nBONDS 21 22'
        self.assertEqual(analogues[0].GetProp('Renderer_Highlight'), render_string1)
        render_string2 = 'COLOR #0000ff\nATOMS 1 2 3 4 5 6 7\nBONDS 1 2 3 4 7 5 6\nCOLOR #ffa500\nATOMS 22\nBONDS 21'
        self.assertEqual(analogues[-1].GetProp('Renderer_Highlight'), render_string2)


if __name__ == '__main__':
    main()
