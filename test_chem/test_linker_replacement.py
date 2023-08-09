# Without adding the data function library to the project content roots PyCharm will not run this
# test unless the data function library is added to sys.paths
# import sys
# sys.path.append(r'C:\Users\david\Desktop\Spotfire12\Modules\Glysade Python_4.2.0.6\Python\Lib\site-packages')
# sys.path.append(r'C:\Users\david\AppData\Local\Temp\ScriptSync\DataFxnPylib\pylib_3.9.7')

# To enable the test to run within PyCharm without explicitly modifying sys.path
# Add C:\Users\david\AppData\Local\Temp\ScriptSync\DataFxnPylib\pylib_3.9.7 as a Project content root and comment the lines above

import json

from collections import defaultdict
from pathlib import Path
from unittest import TestCase, main

from df.LinkerReplacements import LinkerReplacements
from df.chem_helper import column_to_molecules
from df.data_transfer import (DataFunction, DataFunctionRequest,
                              DataFunctionResponse)

from rdkit import Chem


def run_script(in_file: str, test_class: DataFunction) -> DataFunctionResponse:
    with open(in_file, 'r') as fh:
        request_json = fh.read()

    request = DataFunctionRequest.parse_raw(request_json)
    response = test_class.execute(request)
    return response


def run_json(request_json: str, test_class: DataFunction) -> DataFunctionResponse:
    request = DataFunctionRequest.parse_raw(request_json)
    response = test_class.execute(request)
    return response


class ScriptTest(TestCase):

    def test_default_replacement(self) -> None:
        file_in = Path(__file__).parent / 'resources' / 'test_linker_replacement1.json'
        lr = LinkerReplacements()
        response = run_script(file_in, lr)
        self.assertTrue(response)
        self.assertEqual(1, len(response.outputTables))
        self.assertEqual(4, len(response.outputTables[0].columns))
        new_mols = column_to_molecules(response.outputTables[0].columns[2])
        self.assertEqual(667, len(new_mols))
        parent_mols = column_to_molecules(response.outputTables[0].columns[0])
        parent_ids = response.outputTables[0].columns[1].values
        linker_smis = response.outputTables[0].columns[3].values
        self.assertEqual(667, len(parent_mols))
        self.assertEqual(667, len(parent_ids))
        self.assertEqual(667, len(linker_smis))
        exp_par_counts = [91, 87, 90, 153, 91, 0, 112, 43]
        exp_par_smiles = ['CCc1cc(-c2ccccc2)nnc1NCCN1CCOCC1',
                          'CCc1cc(-c2ccccc2)nnc1CCCN1CCOCC1',
                          'CCN1CCCC1Nc1cc(C)c(-c2ccccc2O)nn1',
                          'CCN1CCCC1Cc1cc(C)c(-c2ccccc2O)nn1',
                          'Cc1cc(-c2cccs2)nnc1NCCN1CCOCC1',
                          'CCCc1ccccc1',
                          'CCc1cc(-c2ccccc2)nnc1C(=O)N1CCOCC1',
                          'Cc1cccc(-c2cccc(-c3cccc(F)c3)c2)c1']
        par_counts = defaultdict(int)
        for pm in parent_mols:
            par_counts[Chem.MolToSmiles(pm)] += 1
        for smi, count in zip(exp_par_smiles, exp_par_counts):
            self.assertEqual(count, par_counts[smi])
        self.assertEqual('Mol1', parent_ids[0])
        self.assertEqual('Mol 8', parent_ids[-1])

    def test_match_hbonds(self) -> None:
        file_in = Path(__file__).parent / 'resources' / 'test_linker_replacement1.json'
        with open(file_in, 'r') as fh:
            request_json = fh.read()
        request_dict = json.loads(request_json)
        request_dict['inputFields']['matchHbonds']['data'] = True
        request_json = json.dumps(request_dict)
        lr = LinkerReplacements()
        response = run_json(request_json, lr)

        self.assertTrue(response)
        self.assertEqual(1, len(response.outputTables))
        self.assertEqual(4, len(response.outputTables[0].columns))
        new_mols = column_to_molecules(response.outputTables[0].columns[2])
        self.assertEqual(305, len(new_mols))
        self.assertEqual(305, len(response.outputTables[0].columns[1].values))
        parent_mols = column_to_molecules(response.outputTables[0].columns[0])
        self.assertEqual(305, len(parent_mols))
        exp_par_counts = [57, 30, 56, 33, 57, 0, 46, 26]
        exp_par_smiles = ['CCc1cc(-c2ccccc2)nnc1NCCN1CCOCC1',
                          'CCc1cc(-c2ccccc2)nnc1CCCN1CCOCC1',
                          'CCN1CCCC1Nc1cc(C)c(-c2ccccc2O)nn1',
                          'CCN1CCCC1Cc1cc(C)c(-c2ccccc2O)nn1',
                          'Cc1cc(-c2cccs2)nnc1NCCN1CCOCC1',
                          'CCCc1ccccc1',
                          'CCc1cc(-c2ccccc2)nnc1C(=O)N1CCOCC1',
                          'Cc1cccc(-c2cccc(-c3cccc(F)c3)c2)c1']
        par_counts = defaultdict(int)
        for pm in parent_mols:
            par_counts[Chem.MolToSmiles(pm)] += 1
        for smi, count in zip(exp_par_smiles, exp_par_counts):
            self.assertEqual(count, par_counts[smi])

    def test_length_constraints(self) -> None:
        file_in = Path(__file__).parent / 'resources' / 'test_linker_replacement1.json'
        with open(file_in, 'r') as fh:
            request_json = fh.read()
        request_dict = json.loads(request_json)
        request_dict['inputFields']['plusDeltaBonds']['data'] = 0
        request_dict['inputFields']['minusDeltaBonds']['data'] = 0
        request_json = json.dumps(request_dict)
        lr = LinkerReplacements()
        response = run_json(request_json, lr)

        self.assertTrue(response)
        self.assertEqual(1, len(response.outputTables))
        self.assertEqual(4, len(response.outputTables[0].columns))
        new_mols = column_to_molecules(response.outputTables[0].columns[2])
        self.assertEqual(190, len(new_mols))
        parent_mols = column_to_molecules(response.outputTables[0].columns[0])
        self.assertEqual(190, len(parent_mols))
        self.assertEqual(190, len(response.outputTables[0].columns[1].values))
        exp_par_counts = [39, 36, 9, 24, 39, 0, 18, 25]
        exp_par_smiles = ['CCc1cc(-c2ccccc2)nnc1NCCN1CCOCC1',
                          'CCc1cc(-c2ccccc2)nnc1CCCN1CCOCC1',
                          'CCN1CCCC1Nc1cc(C)c(-c2ccccc2O)nn1',
                          'CCN1CCCC1Cc1cc(C)c(-c2ccccc2O)nn1',
                          'Cc1cc(-c2cccs2)nnc1NCCN1CCOCC1',
                          'CCCc1ccccc1',
                          'CCc1cc(-c2ccccc2)nnc1C(=O)N1CCOCC1',
                          'Cc1cccc(-c2cccc(-c3cccc(F)c3)c2)c1']
        par_counts = defaultdict(int)
        for pm in parent_mols:
            par_counts[Chem.MolToSmiles(pm)] += 1
        for smi, count in zip(exp_par_smiles, exp_par_counts):
            self.assertEqual(count, par_counts[smi])

    def test_all_constraints(self) -> None:
        file_in = Path(__file__).parent / 'resources' / 'test_linker_replacement1.json'
        with open(file_in, 'r') as fh:
            request_json = fh.read()
        request_dict = json.loads(request_json)
        request_dict['inputFields']['matchHbonds']['data'] = True
        request_dict['inputFields']['plusDeltaBonds']['data'] = 0
        request_dict['inputFields']['minusDeltaBonds']['data'] = 0
        request_json = json.dumps(request_dict)
        lr = LinkerReplacements()
        response = run_json(request_json, lr)

        self.assertTrue(response)
        self.assertEqual(1, len(response.outputTables))
        self.assertEqual(4, len(response.outputTables[0].columns))
        new_mols = column_to_molecules(response.outputTables[0].columns[2])
        self.assertEqual(93, len(new_mols))
        self.assertEqual(93, len(response.outputTables[0].columns[1].values))
        parent_mols = column_to_molecules(response.outputTables[0].columns[0])
        self.assertEqual(93, len(parent_mols))
        exp_par_counts = [25, 10, 2, 13, 25, 0, 4, 14]
        exp_par_smiles = ['CCc1cc(-c2ccccc2)nnc1NCCN1CCOCC1',
                          'CCc1cc(-c2ccccc2)nnc1CCCN1CCOCC1',
                          'CCN1CCCC1Nc1cc(C)c(-c2ccccc2O)nn1',
                          'CCN1CCCC1Cc1cc(C)c(-c2ccccc2O)nn1',
                          'Cc1cc(-c2cccs2)nnc1NCCN1CCOCC1',
                          'CCCc1ccccc1',
                          'CCc1cc(-c2ccccc2)nnc1C(=O)N1CCOCC1',
                          'Cc1cccc(-c2cccc(-c3cccc(F)c3)c2)c1']
        par_counts = defaultdict(int)
        for pm in parent_mols:
            par_counts[Chem.MolToSmiles(pm)] += 1
        for smi, count in zip(exp_par_smiles, exp_par_counts):
            self.assertEqual(count, par_counts[smi])

    def test_no_ring_linkers(self) -> None:
        file_in = Path(__file__).parent / 'resources' / 'test_linker_replacement1.json'
        with open(file_in, 'r') as fh:
            request_json = fh.read()
        request_dict = json.loads(request_json)
        request_dict['inputFields']['noRingLinkers']['data'] = True
        request_json = json.dumps(request_dict)
        lr = LinkerReplacements()
        response = run_json(request_json, lr)

        self.assertTrue(response)
        self.assertEqual(1, len(response.outputTables))
        self.assertEqual(4, len(response.outputTables[0].columns))
        new_mols = column_to_molecules(response.outputTables[0].columns[2])
        self.assertEqual(573, len(new_mols))
        self.assertEqual(573, len(response.outputTables[0].columns[1].values))
        parent_mols = column_to_molecules(response.outputTables[0].columns[0])
        self.assertEqual(573, len(parent_mols))
        exp_par_counts = [81, 85, 85, 137, 81, 0, 104, 0]
        exp_par_smiles = ['CCc1cc(-c2ccccc2)nnc1NCCN1CCOCC1',
                          'CCc1cc(-c2ccccc2)nnc1CCCN1CCOCC1',
                          'CCN1CCCC1Nc1cc(C)c(-c2ccccc2O)nn1',
                          'CCN1CCCC1Cc1cc(C)c(-c2ccccc2O)nn1',
                          'Cc1cc(-c2cccs2)nnc1NCCN1CCOCC1',
                          'CCCc1ccccc1',
                          'CCc1cc(-c2ccccc2)nnc1C(=O)N1CCOCC1',
                          'Cc1cccc(-c2cccc(-c3cccc(F)c3)c2)c1']
        par_counts = defaultdict(int)
        for pm in parent_mols:
            par_counts[Chem.MolToSmiles(pm)] += 1
        for smi, count in zip(exp_par_smiles, exp_par_counts):
            self.assertEqual(count, par_counts[smi])

    def test_different_number_linkers(self) -> None:
        file_in = Path(__file__).parent / 'resources' / 'test_linker_replacement2.json'
        lr = LinkerReplacements()
        response = run_script(file_in, lr)
        self.assertTrue(response)
        self.assertEqual(1, len(response.outputTables))
        self.assertEqual(6, len(response.outputTables[0].columns))
        new_mols = column_to_molecules(response.outputTables[0].columns[2])
        self.assertEqual(219, len(new_mols))
        parent_mols = column_to_molecules(response.outputTables[0].columns[0])
        parent_ids = response.outputTables[0].columns[1].values
        linker_smis1 = response.outputTables[0].columns[3].values
        linker_smis2 = response.outputTables[0].columns[4].values
        linker_smis3 = response.outputTables[0].columns[5].values
        self.assertEqual(219, len(parent_mols))
        self.assertEqual(219, len(parent_ids))
        self.assertEqual(219, len(linker_smis1))
        self.assertEqual(219, sum([(ls is not None and len(ls) > 0) for ls in linker_smis1]))
        self.assertEqual(219, len(linker_smis2))
        self.assertEqual(200, sum([(ls is not None and len(ls) > 0) for ls in linker_smis2]))
        self.assertEqual(219, len(linker_smis3))
        self.assertEqual(100, sum([(ls is not None and len(ls) > 0) for ls in linker_smis3]))
        exp_par_counts = [19, 100, 100]
        exp_par_smiles = ['c1ccc(Sc2ccccc2)cc1',
                          'O=C(Nc1cccc(Sc2ccccc2)c1)c1cccc(OC(=O)c2ccccc2)c1',
                          'O=C(c1ccccc1)c1cccc(Cc2ccccc2)c1']
        par_counts = defaultdict(int)
        for pm in parent_mols:
            par_counts[Chem.MolToSmiles(pm)] += 1
        for smi, count in zip(exp_par_smiles, exp_par_counts):
            self.assertEqual(count, par_counts[smi])
        self.assertEqual('Mol1', parent_ids[0])
        self.assertEqual('Mol3', parent_ids[-1])
        # It's not possible to check any other output values due to the
        # random sampling of the truncated sets.


if __name__ == '__main__':
    main()
