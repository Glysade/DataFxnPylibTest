
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
        self.assertEqual(2, len(response.outputTables[0].columns))
        new_mols = column_to_molecules(response.outputTables[0].columns[1])
        self.assertEqual(374, len(new_mols))
        parent_mols = column_to_molecules(response.outputTables[0].columns[0])
        self.assertEqual(374, len(parent_mols))
        exp_par_counts = [47, 52, 60, 98, 47, 0, 70]
        exp_par_smiles = ['CCc1cc(-c2ccccc2)nnc1NCCN1CCOCC1',
                          'CCc1cc(-c2ccccc2)nnc1CCCN1CCOCC1',
                          'CCN1CCCC1Nc1cc(C)c(-c2ccccc2O)nn1',
                          'CCN1CCCC1Cc1cc(C)c(-c2ccccc2O)nn1',
                          'Cc1cc(-c2cccs2)nnc1NCCN1CCOCC1',
                          'CCCc1ccccc1',
                          'CCc1cc(-c2ccccc2)nnc1C(=O)N1CCOCC1']
        par_counts = defaultdict(int)
        for pm in parent_mols:
            par_counts[Chem.MolToSmiles(pm)] += 1
        for smi, count in zip(exp_par_smiles, exp_par_counts):
            self.assertEqual(count, par_counts[smi])

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
        self.assertEqual(2, len(response.outputTables[0].columns))
        new_mols = column_to_molecules(response.outputTables[0].columns[1])
        self.assertEqual(171, len(new_mols))
        parent_mols = column_to_molecules(response.outputTables[0].columns[0])
        self.assertEqual(171, len(parent_mols))
        exp_par_counts = [30, 20, 35, 28, 30, 0, 28]
        exp_par_smiles = ['CCc1cc(-c2ccccc2)nnc1NCCN1CCOCC1',
                          'CCc1cc(-c2ccccc2)nnc1CCCN1CCOCC1',
                          'CCN1CCCC1Nc1cc(C)c(-c2ccccc2O)nn1',
                          'CCN1CCCC1Cc1cc(C)c(-c2ccccc2O)nn1',
                          'Cc1cc(-c2cccs2)nnc1NCCN1CCOCC1',
                          'c1cccc1CCC']
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
        self.assertEqual(2, len(response.outputTables[0].columns))
        new_mols = column_to_molecules(response.outputTables[0].columns[1])
        self.assertEqual(104, len(new_mols))
        parent_mols = column_to_molecules(response.outputTables[0].columns[0])
        self.assertEqual(104, len(parent_mols))
        exp_par_counts = [18, 21, 9, 21, 18, 0, 17]
        exp_par_smiles = ['CCc1cc(-c2ccccc2)nnc1NCCN1CCOCC1',
                          'CCc1cc(-c2ccccc2)nnc1CCCN1CCOCC1',
                          'CCN1CCCC1Nc1cc(C)c(-c2ccccc2O)nn1',
                          'CCN1CCCC1Cc1cc(C)c(-c2ccccc2O)nn1',
                          'Cc1cc(-c2cccs2)nnc1NCCN1CCOCC1',
                          'c1cccc1CCC']
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
        self.assertEqual(2, len(response.outputTables[0].columns))
        new_mols = column_to_molecules(response.outputTables[0].columns[1])
        self.assertEqual(104, len(new_mols))
        parent_mols = column_to_molecules(response.outputTables[0].columns[0])
        self.assertEqual(104, len(parent_mols))
        exp_par_counts = [11, 6, 2, 11, 11, 0, 4]
        exp_par_smiles = ['CCc1cc(-c2ccccc2)nnc1NCCN1CCOCC1',
                          'CCc1cc(-c2ccccc2)nnc1CCCN1CCOCC1',
                          'CCN1CCCC1Nc1cc(C)c(-c2ccccc2O)nn1',
                          'CCN1CCCC1Cc1cc(C)c(-c2ccccc2O)nn1',
                          'Cc1cc(-c2cccs2)nnc1NCCN1CCOCC1',
                          'c1cccc1CCC']
        par_counts = defaultdict(int)
        for pm in parent_mols:
            par_counts[Chem.MolToSmiles(pm)] += 1
        for smi, count in zip(exp_par_smiles, exp_par_counts):
            self.assertEqual(count, par_counts[smi])


if __name__ == '__main__':
    main()
