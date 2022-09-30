import json
import random
from collections import defaultdict
from pathlib import Path
from unittest import TestCase, main

from df.chem_helper import column_to_molecules
from df.data_transfer import DataFunction, DataFunctionRequest, DataFunctionResponse
from rdkit import Chem

from df.RGroupReplacement import RGroupReplacement


def run_script(in_file: str, test_class: DataFunction) -> DataFunctionResponse:
    with open(in_file, 'r') as fh:
        request_json = fh.read()

    request = DataFunctionRequest.parse_raw(request_json)
    response = test_class.execute(request)
    return response


def get_core_atom_cds(mol):
    # the atom properties have been lost by now, so dig the core atoms
    # out of the atom colouring
    lines = mol.GetProp('Renderer_Highlight').split('\n')
    core_ats = lines[1][6:].split()
    at_cds = []
    for ca in core_ats:
        ica = int(ca) - 1
        at_cds.append(mol.GetConformer().GetAtomPosition(ica))
    return at_cds


def calc_core_rmses(mols: list[Chem.Mol]) -> float:
    # check the cores are all lined up correctly.  Coordgen doesn't
    # make them exactly the same, but very close
    core_cds = get_core_atom_cds(mols[0])

    rmses = []
    for i, mol in enumerate(mols[1:], 1):
        this_core_cds = get_core_atom_cds(mol)
        rms = 0.0
        for cd, tcd in zip(core_cds, this_core_cds):
            rms += (cd.x - tcd.x) * (cd.x - tcd.x) + (cd.y - tcd.y) * (cd.y - tcd.y)
        rmses.append(rms)
    return max(rmses)

def analogues_in_parents(parents, parent_ids, analogues) -> list[str]:
    same_smis = []
    for par_id, par, an in zip(parent_ids, parents, analogues):
        par_smi = Chem.MolToSmiles(par)
        mol_smi = Chem.MolToSmiles(an)
        if par_smi == mol_smi:
            same_smis.append(par_id)
    return same_smis

class ScriptTest(TestCase):

    def test_layer1(self) -> None:
        file_in = Path(__file__).parent / 'resources' / 'test_r_group_replacement1.json'
        with open(file_in, 'r') as fh:
            request_json = fh.read()

        request_dict = json.loads(request_json)
        request_json = json.dumps(request_dict)
        rgr = RGroupReplacement()
        request = DataFunctionRequest.parse_raw(request_json)
        response = rgr.execute(request)

        self.assertTrue(response)
        parents = column_to_molecules(response.outputTables[0].columns[0])
        parent_ids = response.outputTables[0].columns[1].values
        mols = column_to_molecules(response.outputTables[0].columns[3])
        self.assertEqual(len(parents), 105)
        self.assertEqual(len(parent_ids), 105)
        self.assertEqual(len(response.outputTables[0].columns[2].values), 105)
        self.assertEqual(len(mols), 105)
        self.assertListEqual(response.outputColumns[0].values, [0, 5, 5, 4, 35, 5, 5, 35, 11])
        parent_counts = defaultdict(int)
        for par_id in response.outputTables[0].columns[1].values:
            parent_counts[par_id] += 1
        for i in range(9):
            self.assertEqual(parent_counts[f'Mol{i+1}'], response.outputColumns[0].values[i])
        self.assertEqual(Chem.MolToSmiles(mols[0]), 'Cc1ccccc1Cl')
        self.assertEqual(Chem.MolToSmiles(mols[-1]), 'COc1ccc(C)c(OC)n1')
        self.assertLess(calc_core_rmses(mols), 0.005)
        same_smis = analogues_in_parents(parents, parent_ids, mols)
        self.assertFalse(same_smis, f'Molecule(s) had same SMILES as parent : {" ".join(same_smis)}')

    def test_layer1_plus_layer2(self) -> None:
        file_in = Path(__file__).parent / 'resources' / 'test_r_group_replacement1.json'
        with open(file_in, 'r') as fh:
            request_json = fh.read()

        request_dict = json.loads(request_json)
        request_dict['inputFields']['useLayer2']['data'] = True
        request_json = json.dumps(request_dict)
        rgr = RGroupReplacement()
        request = DataFunctionRequest.parse_raw(request_json)
        response = rgr.execute(request)

        self.assertTrue(response)
        parents = column_to_molecules(response.outputTables[0].columns[0])
        parent_ids = response.outputTables[0].columns[1].values
        mols = column_to_molecules(response.outputTables[0].columns[3])

        self.assertEqual(len(parents), 431)
        self.assertEqual(len(parent_ids), 431)
        self.assertEqual(len(response.outputTables[0].columns[2].values), 431)
        self.assertEqual(len(mols), 431)
        self.assertListEqual(response.outputColumns[0].values, [0, 13, 12, 4, 209, 13, 14, 111, 55])
        parent_counts = defaultdict(int)
        for par_id in response.outputTables[0].columns[1].values:
            parent_counts[par_id] += 1
        for i in range(9):
            self.assertEqual(parent_counts[f'Mol{i+1}'], response.outputColumns[0].values[i])
        self.assertEqual(Chem.MolToSmiles(mols[0]), 'Cc1ccccc1Cl')
        self.assertEqual(Chem.MolToSmiles(mols[-1]), 'Cc1ccc(O)nc1S(=O)(=O)N(C)C')
        self.assertLess(calc_core_rmses(mols), 0.005)
        same_smis = analogues_in_parents(parents, parent_ids, mols)
        self.assertFalse(same_smis, f'Molecule(s) had same SMILES as parent : {" ".join(same_smis)}')

    def test_gyrase_multicore_decomp(self) -> None:
        file_in = Path(__file__).parent / 'resources' / 'test_r_group_replacement2.json'
        rgr = RGroupReplacement()
        response = run_script(file_in, rgr)
        self.assertTrue(response)
        self.assertEqual(len(response.outputColumns[0].values), 981)
        # the sum of the output column values (the number of analogues
        # each parent produced) should be the same as the lengths of
        # the columns in the output table
        self.assertEqual(sum(response.outputColumns[0].values), 1781)
        parents = column_to_molecules(response.outputTables[0].columns[0])
        parent_ids = response.outputTables[0].columns[1].values
        mols = column_to_molecules(response.outputTables[0].columns[3])

        same_smis = analogues_in_parents(parents, parent_ids, mols)
        self.assertFalse(same_smis, f'Molecule(s) had same SMILES as parent : {" ".join(same_smis)}')
        self.assertEqual(len(parents), 1781)
        self.assertEqual(len(parent_ids), 1781)
        self.assertEqual(len(mols), 1781)
        self.assertEqual(Chem.MolToSmiles(mols[0]), 'Cc1[nH]c(C(=O)NC2CCN(c3ccccn3)CC2)cc1Br')
        self.assertEqual(Chem.MolToSmiles(mols[-1]),
                         'CNC(=O)c1cc2c(-n3ccc(C(F)(F)F)n3)c(-c3cc(-c4n[nH]c(=O)o4)cnc3OC)cnc2[nH]1')
        self.assertLess(calc_core_rmses(mols), 0.005)

    def test_gyrase_multicore_decomp_full(self) -> None:
        file_in = Path(__file__).parent / 'resources' / 'test_r_group_replacement2.json'
        with open(file_in, 'r') as fh:
            request_json = fh.read()

        request_dict = json.loads(request_json)
        request_dict['inputFields']['useLayer2']['data'] = True
        request_dict['inputFields']['incOrigRGroups']['data'] = True
        request_json = json.dumps(request_dict)
        rgr = RGroupReplacement()
        request = DataFunctionRequest.parse_raw(request_json)
        response = rgr.execute(request)

        self.assertTrue(response)
        self.assertEqual(len(response.outputColumns[0].values), 981)
        # the sum of the output column values (the number of analogues
        # each parent produced) should be the same as the lengths of
        # the columns in the output table
        self.assertEqual(sum(response.outputColumns[0].values), 4072)
        parents = column_to_molecules(response.outputTables[0].columns[0])
        parent_ids = response.outputTables[0].columns[1].values
        mols = column_to_molecules(response.outputTables[0].columns[3])
        self.assertEqual(len(parents), 4072)
        self.assertEqual(len(parent_ids), 4072)
        self.assertEqual(len(mols), 4072)
        self.assertEqual(Chem.MolToSmiles(mols[0]), 'Cc1[nH]c(C(=O)NC2CCN(c3ccccn3)CC2)cc1Br')
        self.assertEqual(Chem.MolToSmiles(mols[-1]),
                         'COc1ncc(-c2n[nH]c(=O)o2)cc1-c1cnc2[nH]c(C(=O)N(C)C)cc2c1-n1ccc(C(F)(F)F)n1')
        same_smis = analogues_in_parents(parents, parent_ids, mols)
        self.assertFalse(same_smis, f'Molecule(s) had same SMILES as parent : {" ".join(same_smis)}')
        self.assertLess(calc_core_rmses(mols), 0.005)


if __name__ == '__main__':
    main()
