import json
from collections import defaultdict
from pathlib import Path
from unittest import TestCase, main

from df.chem_helper import column_to_molecules
from df.data_transfer import DataFunction, DataFunctionRequest, DataFunctionResponse
from rdkit import Chem

from df.RGroupReplacement import RGroupReplacement

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


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


# columns used in the output table.  So that when they are moved around
# or added to in the DataFxn code, there's only one place they need
#changing in the tests.
PAR_COL = 0
PAR_IDS_COL = 1
MOLS_COL = 2
CORES_COL = 5
CORE_NUMS_COL = 4
R_CHG_COL = 3


class ScriptTest(TestCase):

    def test_layer1(self) -> None:
        file_in = Path(__file__).parent / 'resources' / 'test_r_group_replacement1.json'
        rgr = RGroupReplacement()
        response = run_script(file_in, rgr)
        self.assertTrue(response)

        parents = column_to_molecules(response.outputTables[0].columns[PAR_COL])
        parent_ids = response.outputTables[0].columns[PAR_IDS_COL].values
        mols = column_to_molecules(response.outputTables[0].columns[MOLS_COL])
        changed_rgroups = response.outputTables[0].columns[R_CHG_COL].values
        core_nums = response.outputTables[0].columns[CORE_NUMS_COL].values

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
        mol1_highs = 'COLOR #00bfff\nATOMS\nBONDS 8'
        self.assertEqual(mols[0].GetProp('Renderer_Highlight'), mol1_highs)
        molm1_highs = 'COLOR #00bfff\nATOMS\nBONDS 9 11'
        self.assertEqual(mols[-1].GetProp('Renderer_Highlight'), molm1_highs)
        self.assertEqual(changed_rgroups[0], 'R5')
        self.assertEqual(parent_ids[34], 'Mol5')
        self.assertEqual(changed_rgroups[34], 'R4:R5')
        self.assertEqual(changed_rgroups[-1], 'R2')
        self.assertEqual(core_nums.count('1'), 54)
        self.assertEqual(core_nums.count('2'), 51)

    def test_layer2(self) -> None:
        file_in = Path(__file__).parent / 'resources' / 'test_r_group_replacement1.json'
        with open(file_in, 'r') as fh:
            request_json = fh.read()

        request_dict = json.loads(request_json)
        request_dict['inputFields']['useLayer1']['data'] = False
        request_dict['inputFields']['useLayer2']['data'] = True
        request_json = json.dumps(request_dict)
        rgr = RGroupReplacement()
        request = DataFunctionRequest.parse_raw(request_json)
        response = rgr.execute(request)

        self.assertTrue(response)
        parents = column_to_molecules(response.outputTables[0].columns[PAR_COL])
        parent_ids = response.outputTables[0].columns[PAR_IDS_COL].values
        mols = column_to_molecules(response.outputTables[0].columns[MOLS_COL])
        changed_rgroups = response.outputTables[0].columns[R_CHG_COL].values
        core_nums = response.outputTables[0].columns[CORE_NUMS_COL].values

        self.assertEqual(len(parents), 173)
        self.assertEqual(len(parent_ids), 173)
        self.assertEqual(len(response.outputTables[0].columns[2].values), 173)
        self.assertEqual(len(mols), 173)
        self.assertListEqual(response.outputColumns[0].values, [0, 8, 7, 0, 89, 8, 9, 26, 26])
        parent_counts = defaultdict(int)
        for par_id in response.outputTables[0].columns[1].values:
            parent_counts[par_id] += 1
        for i in range(9):
            self.assertEqual(parent_counts[f'Mol{i+1}'], response.outputColumns[0].values[i])
        self.assertEqual(Chem.MolToSmiles(mols[0]), 'Cc1ccccc1F')
        self.assertEqual(Chem.MolToSmiles(mols[-1]), 'Cc1ccc(O)nc1S(=O)(=O)N(C)C')
        self.assertLess(calc_core_rmses(mols), 0.005)
        same_smis = analogues_in_parents(parents, parent_ids, mols)
        self.assertFalse(same_smis, f'Molecule(s) had same SMILES as parent : {" ".join(same_smis)}')
        # check for a level 2 highlight
        mol6highs = 'COLOR #ffbf00\nATOMS\nBONDS 8 9 10 11 12'
        self.assertEqual(mols[6].GetProp('Renderer_Highlight'), mol6highs)
        self.assertEqual(changed_rgroups[0], 'R5')
        self.assertEqual(parent_ids[121], 'Mol8')
        self.assertEqual(changed_rgroups[132], 'R1:R3')
        self.assertEqual(changed_rgroups[-1], 'R1:R2')
        self.assertEqual(core_nums.count('1'), 112)
        self.assertEqual(core_nums.count('2'), 61)

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
        parents = column_to_molecules(response.outputTables[0].columns[PAR_COL])
        parent_ids = response.outputTables[0].columns[PAR_IDS_COL].values
        mols = column_to_molecules(response.outputTables[0].columns[MOLS_COL])
        changed_rgroups = response.outputTables[0].columns[R_CHG_COL].values
        core_nums = response.outputTables[0].columns[CORE_NUMS_COL].values

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
        # check for a level 2 highlight
        mol6highs = 'COLOR #ffbf00\nATOMS\nBONDS 8'
        self.assertEqual(mols[6].GetProp('Renderer_Highlight'), mol6highs)
        self.assertEqual(changed_rgroups[0], 'R5')
        self.assertEqual(parent_ids[273], 'Mol8')
        self.assertEqual(changed_rgroups[273], 'R1:R3')
        self.assertEqual(changed_rgroups[-1], 'R1:R2')
        self.assertEqual(core_nums.count('1'), 251)
        self.assertEqual(core_nums.count('2'), 180)

    def test_gyrase_multicore_decomp(self) -> None:
        file_in = Path(__file__).parent / 'resources' / 'test_r_group_replacement2.json'
        rgr = RGroupReplacement()
        response = run_script(file_in, rgr)
        self.assertTrue(response)
        self.assertEqual(len(response.outputColumns[0].values), 981)
        # the sum of the output column values (the number of analogues
        # each parent produced) should be the same as the lengths of
        # the columns in the output table
        self.assertEqual(sum(response.outputColumns[0].values), 1804)
        parents = column_to_molecules(response.outputTables[0].columns[PAR_COL])
        parent_ids = response.outputTables[0].columns[PAR_IDS_COL].values
        mols = column_to_molecules(response.outputTables[0].columns[MOLS_COL])
        changed_rgroups = response.outputTables[0].columns[R_CHG_COL].values
        core_nums = response.outputTables[0].columns[CORE_NUMS_COL].values

        same_smis = analogues_in_parents(parents, parent_ids, mols)
        self.assertFalse(same_smis, f'Molecule(s) had same SMILES as parent : {" ".join(same_smis)}')
        self.assertEqual(len(parents), 1804)
        self.assertEqual(len(parent_ids), 1804)
        self.assertEqual(len(mols), 1804)
        self.assertEqual(Chem.MolToSmiles(mols[0]), 'Cc1[nH]c(C(=O)NC2CCN(c3ccccn3)CC2)cc1Br')
        self.assertEqual(Chem.MolToSmiles(mols[-1]),
                         'CNC(=O)c1cc2c(-n3ccc(C(F)(F)F)n3)c(-c3cc(-c4n[nH]c(=O)o4)cnc3OC)cnc2[nH]1')
        self.assertLess(calc_core_rmses(mols), 0.005)
        self.assertEqual(changed_rgroups[0], 'R3')
        self.assertEqual(parent_ids[1059], 'AZ1669')
        self.assertEqual(changed_rgroups[1059], 'R2:R8')
        self.assertEqual(changed_rgroups[-1], 'R1')
        self.assertListEqual(sorted(list(set(core_nums))), ['1', '11', '5', '6', '7', '8'])
        self.assertEqual(core_nums.count('1'), 1751)
        self.assertEqual(core_nums.count('2'), 0)
        self.assertEqual(core_nums.count('5'), 17)
        self.assertEqual(core_nums.count('6'), 1)
        self.assertEqual(core_nums.count('7'), 22)
        self.assertEqual(core_nums.count('8'), 1)
        self.assertEqual(core_nums.count('11'), 12)

    def test_gyrase_multicore_decomp_full(self) -> None:
        file_in = Path(__file__).parent / 'resources' / 'test_r_group_replacement2.json'
        with open(file_in, 'r') as fh:
            request_json = fh.read()

        request_dict = json.loads(request_json)
        request_dict['inputFields']['useLayer2']['data'] = True
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
        parents = column_to_molecules(response.outputTables[0].columns[PAR_COL])
        parent_ids = response.outputTables[0].columns[PAR_IDS_COL].values
        mols = column_to_molecules(response.outputTables[0].columns[MOLS_COL])
        changed_rgroups = response.outputTables[0].columns[R_CHG_COL].values
        core_nums = response.outputTables[0].columns[CORE_NUMS_COL].values

        self.assertEqual(len(parents), 4072)
        self.assertEqual(len(parent_ids), 4072)
        self.assertEqual(len(mols), 4072)
        self.assertEqual(Chem.MolToSmiles(mols[0]), 'Cc1[nH]c(C(=O)NC2CCN(c3ccccn3)CC2)cc1Br')
        self.assertEqual(Chem.MolToSmiles(mols[-1]),
                         'COc1ncc(-c2n[nH]c(=O)o2)cc1-c1cnc2[nH]c(C(=O)N(C)C)cc2c1-n1ccc(C(F)(F)F)n1')
        same_smis = analogues_in_parents(parents, parent_ids, mols)
        self.assertFalse(same_smis, f'Molecule(s) had same SMILES as parent : {" ".join(same_smis)}')
        self.assertLess(calc_core_rmses(mols), 0.005)
        self.assertEqual(changed_rgroups[0], 'R3')
        self.assertEqual(parent_ids[1059], 'AZ1521')
        self.assertEqual(changed_rgroups[1059], 'R2')
        self.assertEqual(changed_rgroups[-1], 'R1')
        self.assertListEqual(sorted(list(set(core_nums))), ['1', '11', '5', '6', '7', '8'])
        self.assertEqual(core_nums.count('1'), 3939)
        self.assertEqual(core_nums.count('2'), 0)
        self.assertEqual(core_nums.count('5'), 38)
        self.assertEqual(core_nums.count('6'), 1)
        self.assertEqual(core_nums.count('7'), 63)
        self.assertEqual(core_nums.count('8'), 3)
        self.assertEqual(core_nums.count('11'), 28)

    def test_pauls_aromatic_rgroup_bug(self) -> None:
        file_in = Path(__file__).parent / 'resources' / 'test_r_group_replacement3.json'
        rgr = RGroupReplacement()
        response = run_script(file_in, rgr)
        self.assertTrue(response)


if __name__ == '__main__':
    main()
