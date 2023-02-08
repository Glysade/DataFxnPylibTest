from pathlib import Path
from typing import Callable
from unittest import TestCase, main

from df.chem_helper import column_to_molecules
from df.data_transfer import DataFunctionRequest, DataFunctionResponse
from rdkit import Chem


def run_script(in_file: str, execute: Callable[[DataFunctionRequest], DataFunctionResponse]) -> DataFunctionResponse:
    with open(in_file, 'r') as fh:
        request_json = fh.read()

    return run_script_json(request_json, execute)


def run_script_json(in_json: dict,
                    execute: Callable[[DataFunctionRequest], DataFunctionResponse]) -> DataFunctionResponse:
    request = DataFunctionRequest.parse_raw(in_json)
    response = execute(request)
    return response


def read_mols(mol_file: str) -> list[Chem.Mol]:
    suppl = None
    if mol_file.endswith('.smi'):
        suppl = Chem.SmilesMolSupplier(mol_file, titleLine=False)
    elif mol_file.endswith('.sdf'):
        suppl = Chem.SDMolSupplier(mol_file)
    mols = []
    for mol in suppl:
        if mol:
            mols.append(mol)
    return mols


class ScriptTest(TestCase):

    def test_script(self) -> None:
        from MaximumCommonSubstructure_script import execute
        file_in = Path(__file__).parent / 'resources' / 'test_maximum_common_substructure.json'
        response = run_script(file_in, execute)
        self.assertTrue(response)
        self.assertEqual(len(response.outputColumns), 0)
        self.assertEqual(len(response.outputTables), 1)
        self.assertEqual(len(response.outputTables[0].columns), 8)
        self.assertEqual(len(response.outputTables[0].columns[0].values), 99)
        num_mcs = len(set(response.outputTables[0].columns[0].values))
        self.assertEqual(num_mcs, 19)
        self.assertEqual(response.outputTables[0].columns[0].values.count(1), 15)
        self.assertEqual(response.outputTables[0].columns[0].values.count(19), 2)
        self.assertEqual(response.outputTables[0].columns[5].values[0],
                         '[#6]1:&@[#6](:&@[#7]:&@[#6](:&@[#7]:&@[#6]:&@1)-&!@[#7])-&!@[#7]')
        self.assertEqual(response.outputTables[0].columns[5].values[-1],
                         '[#6]1:&@[#6](-&!@[#8]-&!@[#6]-&!@[#6]):&@[#6]:&@[#6]:&@[#6]:&@[#6]:&@1')
        self.assertEqual(response.outputTables[0].columns[4].values[0], 15)
        self.assertEqual(response.outputTables[0].columns[4].values[-1], 2)
        mols = column_to_molecules(response.outputTables[0].columns[6])
        self.assertEqual(mols[0].GetProp('Renderer_Highlight'),
                         'COLOR #ff0000\nATOMS\nBONDS 1 2 4 5 6 7 8 9')
        self.assertEqual(mols[-1].GetProp('Renderer_Highlight'),
                         'COLOR #ff0000\nATOMS\nBONDS 24 23 37 31 35 33 32 39 27')
        ids = response.outputTables[0].columns[7].values
        self.assertEqual(ids[0], '1kmv_lig_LII')
        self.assertEqual(ids[15], '1kmv_lig_LII')
        self.assertEqual(ids[-1], '3s7a_lig_684')

    def test_script1(self) -> None:
        from MaximumCommonSubstructure_script import execute
        file_in = Path(__file__).parent / 'resources' / 'test_maximum_common_substructure1.json'
        response = run_script(file_in, execute)
        self.assertTrue(response)
        self.assertEqual(response.outputTables[0].tableName,
                         'Structure <MOLFILE> MCSs EXHAUSTIVE')
        self.assertEqual(len(response.outputTables[0].columns[0].values), 33)
        self.assertEqual(response.outputTables[0].columns[4].values[0], 7)
        self.assertEqual(response.outputTables[0].columns[4].values[-1], 6)

    def test_script1_timeout(self) -> None:
        from MaximumCommonSubstructure_script import execute
        file_in = Path(__file__).parent / 'resources' / 'test_maximum_common_substructure1.json'
        with open(file_in, 'r') as fh:
            request_json = fh.read()
        import json
        json_dict = json.loads(request_json)
        json_dict['inputFields']['maxTime']['data'] = 3
        request_json = json.dumps(json_dict)
        response = run_script_json(request_json, execute)
        self.assertTrue(response)
        self.assertEqual(response.outputTables[0].tableName,
                         'Structure <MOLFILE> PARTIAL MCSs EXHAUSTIVE')
        self.assertEqual(len(response.outputTables[0].columns[0].values), 13)
        self.assertEqual(response.outputTables[0].columns[4].values[0], 7)
        self.assertEqual(response.outputTables[0].columns[4].values[-1], 6)

    def test_script1_greedy(self) -> None:
        from MaximumCommonSubstructure_script import execute
        file_in = Path(__file__).parent / 'resources' / 'test_maximum_common_substructure1.json'
        with open(file_in, 'r') as fh:
            request_json = fh.read()
        import json
        json_dict = json.loads(request_json)
        json_dict['inputFields']['method']['data'] = 'GREEDY'
        request_json = json.dumps(json_dict)
        response = run_script_json(request_json, execute)
        self.assertEqual(response.outputTables[0].tableName,
                         'Structure <MOLFILE> MCSs GREEDY')
        self.assertEqual(len(response.outputTables[0].columns[0].values), 12)
        self.assertEqual(response.outputTables[0].columns[4].values[0], 6)
        self.assertEqual(response.outputTables[0].columns[4].values[-1], 6)

    def test_timeout(self) -> None:
        from MaximumCommonSubstructure_script import findMCSs
        in_file = Path('__file__').parent / 'resources' / 'chembl_30_50.smi'
        mols = read_mols(str(in_file))
        mcss, timed_out = findMCSs(mols, 3, 6, 6, True, 60, 7)
        self.assertTrue(timed_out)
        self.assertEqual(len(mcss), 17)
        self.assertEqual(mcss[0]['numMols'], 14)
        self.assertEqual(mcss[-1]['numMols'], 3)

    def test_greedy(self) -> None:
        from MaximumCommonSubstructure_script import findMCSs
        in_file = Path('__file__').parent / 'resources' / 'P00374_2d.sdf'
        mols = read_mols(str(in_file))
        self.assertEqual(len(mols), 15)
        mcss, timed_out = findMCSs(mols, 3, 6, 6, True, 60, -1, 'GREEDY')
        self.assertFalse(timed_out)
        self.assertEqual(len(mcss), 5)
        self.assertEqual(mcss[0]['numMols'], 15)
        self.assertEqual(mcss[-1]['numMols'], 3)

    def test_exhaustive(self) -> None:
        from MaximumCommonSubstructure_script import findMCSs
        in_file = Path('__file__').parent / 'resources' / 'P00374_2d.sdf'
        mols = read_mols(str(in_file))
        mcss, timed_out = findMCSs(mols, 3, 6, 6, True, 60, -1, 'EXHAUSTIVE')
        self.assertFalse(timed_out)
        self.assertEqual(len(mcss), 16)
        self.assertEqual(mcss[0]['numMols'], 15)
        self.assertEqual(mcss[-1]['numMols'], 3)


if __name__ == '__main__':
    main()
