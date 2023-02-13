import re
from pathlib import Path
from typing import Callable
from unittest import TestCase, main

from df.chem_helper import column_to_molecules
from df.data_transfer import DataFunctionRequest, DataFunctionResponse
from yaml_to_script import extract_script

from rdkit import Chem

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


def run_script(in_file: str, execute: Callable[[DataFunctionRequest], DataFunctionResponse]) -> DataFunctionResponse:
    with open(in_file, 'r') as fh:
        request_json = fh.read()

    return run_script_json(request_json, execute)


def run_script_json(in_json: dict,
                    execute: Callable[[DataFunctionRequest], DataFunctionResponse]) -> DataFunctionResponse:
    request = DataFunctionRequest.parse_raw(in_json)
    response = execute(request)
    return response


def read_mols(mol_file: str) -> tuple[list[Chem.Mol], list[str]]:
    suppl = None
    if mol_file.endswith('.smi'):
        suppl = Chem.SmilesMolSupplier(mol_file, titleLine=False)
    elif mol_file.endswith('.sdf'):
        suppl = Chem.SDMolSupplier(mol_file)
    mols = []
    ids = []
    for i, mol in enumerate(suppl):
        if mol:
            mols.append(mol)
            try:
                ids.append(mol.GetProp('_Name'))
            except KeyError:
                ids.append(f'Str_{i}')
    return mols, ids


class ScriptTest(TestCase):

    def setUp(self) -> None:
        """
        It's important when doing production testing to use the script
        that's in the YAML, which is in the DataFxn repo, as that's the
        one that Spotfire will be running.  When developing, it's
        convenient to edit the script in this directory, transferring
        it to the YAML file once it's ready for production.
        If there's a script MaximumCommonSubstructure_script_dev.py,
        use that, otherwise make one from the YAML file and use that.
        MaximumCommonSubstructure_script_dev.py should not be in the
        repo, only ever in a development branch so testing a production
        repo will use the script derived from the YAML.
        Doing it this way allows for people editing the YAML directly,
        which is handing for minor tweaks, and those tweaks always
        being tested.
        Assume the DataFxn repo has been cloned alongside this one.
        """
        this_dir = Path(__file__).parent
        script_file = this_dir / 'MaximumCommonSubstructure_script_dev.py'
        if Path(script_file).exists():
            print(f'Using development script')
            from MaximumCommonSubstructure_script_dev import execute, findMCSs
            self._script_file = None
        else:
            data_fxn_dir = this_dir.parent.parent / 'DataFxns'
            print(f'Using production script in {data_fxn_dir}')
            yaml_file = data_fxn_dir / 'python' / 'local' / 'MaximumCommonSubstructure.yaml'
            self._script_file = this_dir / 'MaximumCommonSubstructure_script.py'
            if yaml_file.exists():
                script_lines = extract_script(yaml_file)
                with open(self._script_file, 'w') as f:
                    f.write(''.join(script_lines))

            from MaximumCommonSubstructure_script import execute, findMCSs
        self._execute = execute
        self._findMCSs = findMCSs

    def tearDown(self) -> None:
        # if a script file was made, tidy it up
        if self._script_file is not None and self._script_file.exists():
            self._script_file.unlink()

    def test_script(self) -> None:
        file_in = Path(__file__).parent / 'resources' / 'test_maximum_common_substructure.json'
        response = run_script(file_in, self._execute)
        self.assertTrue(response)
        self.assertEqual(response.outputTables[0].tableName, 'MCSs EXHAUSTIVE')
        self.assertEqual(len(response.outputColumns), 0)
        self.assertEqual(len(response.outputTables), 1)
        self.assertEqual(len(response.outputTables[0].columns), 8)
        self.assertEqual(len(response.outputTables[0].columns[0].values), 99)
        num_mcs = len(set(response.outputTables[0].columns[0].values))
        self.assertEqual(num_mcs, 19)
        self.assertEqual(response.outputTables[0].columns[0].values.count(1), 15)
        self.assertEqual(response.outputTables[0].columns[0].values.count(19), 2)
        self.assertEqual(response.outputTables[0].columns[6].values[0],
                         '[#6]1:&@[#6](:&@[#7]:&@[#6](:&@[#7]:&@[#6]:&@1)-&!@[#7])-&!@[#7]')
        self.assertEqual(response.outputTables[0].columns[6].values[-1],
                         '[#6]1:&@[#6](-&!@[#8]-&!@[#6]-&!@[#6]):&@[#6]:&@[#6]:&@[#6]:&@[#6]:&@1')
        self.assertEqual(response.outputTables[0].columns[3].values[0], 15)
        self.assertEqual(response.outputTables[0].columns[3].values[-1], 2)
        mols = column_to_molecules(response.outputTables[0].columns[2])
        self.assertEqual(mols[0].GetProp('Renderer_Highlight'),
                         'COLOR #14aadb\nATOMS\nBONDS 1 2 4 5 6 7 8 9')
        self.assertEqual(mols[-1].GetProp('Renderer_Highlight'),
                         'COLOR #14aadb\nATOMS\nBONDS 24 23 37 31 35 33 32 39 27')
        ids = response.outputTables[0].columns[7].values
        self.assertEqual(ids[0], '1kmv_lig_LII')
        self.assertEqual(ids[15], '1kmv_lig_LII')
        self.assertEqual(ids[-1], '3s7a_lig_684')

    def test_script1(self) -> None:
        file_in = Path(__file__).parent / 'resources' / 'test_maximum_common_substructure1.json'
        response = run_script(file_in, self._execute)
        self.assertTrue(response)
        self.assertEqual(response.outputTables[0].tableName,
                         'MCSs EXHAUSTIVE')
        self.assertEqual(len(response.outputTables[0].columns[0].values), 33)
        self.assertEqual(response.outputTables[0].columns[3].values[0], 7)
        self.assertEqual(response.outputTables[0].columns[3].values[-1], 6)

    def test_script1_timeout(self) -> None:
        file_in = Path(__file__).parent / 'resources' / 'test_maximum_common_substructure1.json'
        with open(file_in, 'r') as fh:
            request_json = fh.read()
        import json
        json_dict = json.loads(request_json)
        json_dict['inputFields']['maxTime']['data'] = 3
        request_json = json.dumps(json_dict)
        response = run_script_json(request_json, self._execute)
        self.assertTrue(response)
        self.assertEqual(response.outputTables[0].tableName,
                         'PARTIAL MCSs EXHAUSTIVE')
        self.assertEqual(len(response.outputTables[0].columns[0].values), 13)
        self.assertEqual(response.outputTables[0].columns[3].values[0], 7)
        self.assertEqual(response.outputTables[0].columns[3].values[-1], 6)

    def test_script1_greedy(self) -> None:
        file_in = Path(__file__).parent / 'resources' / 'test_maximum_common_substructure1.json'
        with open(file_in, 'r') as fh:
            request_json = fh.read()
        import json
        json_dict = json.loads(request_json)
        json_dict['inputFields']['method']['data'] = 'GREEDY'
        request_json = json.dumps(json_dict)
        response = run_script_json(request_json, self._execute)
        self.assertEqual(response.outputTables[0].tableName,
                         'MCSs GREEDY')
        self.assertEqual(len(response.outputTables[0].columns[0].values), 12)
        self.assertEqual(response.outputTables[0].columns[3].values[0], 6)
        self.assertEqual(response.outputTables[0].columns[3].values[-1], 6)

    def test_timeout(self) -> None:
        in_file = Path('__file__').parent / 'resources' / 'chembl_30_50.smi'
        mols, ids = read_mols(str(in_file))
        mcss, timed_out = self._findMCSs(mols, 3, 6, 6, ids, True, 60, 7)
        self.assertTrue(timed_out)
        self.assertEqual(len(mcss), 17)
        self.assertEqual(mcss[0]['numMols'], 14)
        self.assertEqual(mcss[-1]['numMols'], 3)

    def test_greedy(self) -> None:
        in_file = Path('__file__').parent / 'resources' / 'P00374_2d.sdf'
        mols, ids = read_mols(str(in_file))
        self.assertEqual(len(mols), 15)
        mcss, timed_out = self._findMCSs(mols, 3, 6, 6, ids, True, 60, -1,
                                         'GREEDY')
        self.assertFalse(timed_out)
        self.assertEqual(len(mcss), 5)
        self.assertEqual(mcss[0]['numMols'], 15)
        self.assertEqual(mcss[-1]['numMols'], 3)

    def test_exhaustive(self) -> None:
        in_file = Path('__file__').parent / 'resources' / 'P00374_2d.sdf'
        mols, ids = read_mols(str(in_file))
        mcss, timed_out = self._findMCSs(mols, 3, 6, 6, ids, True, 60, -1,
                                         'EXHAUSTIVE')
        self.assertFalse(timed_out)
        self.assertEqual(len(mcss), 16)
        self.assertEqual(mcss[0]['numMols'], 15)
        self.assertEqual(mcss[-1]['numMols'], 3)


if __name__ == '__main__':
    main()
