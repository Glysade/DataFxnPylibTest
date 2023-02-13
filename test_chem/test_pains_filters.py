import inspect
import os

from pathlib import Path
from typing import Callable, Union
from unittest import TestCase, main

from df.chem_helper import column_to_molecules
from df.data_transfer import DataFunctionRequest, DataFunctionResponse

from rdkit import Chem

from yaml_to_script import extract_script

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


def run_script(in_file: str, execute: Callable[[DataFunctionRequest], DataFunctionResponse]) -> DataFunctionResponse:
    with open(in_file, 'r') as fh:
        request_json = fh.read()
    
    request = DataFunctionRequest.parse_raw(request_json)  
    response = execute(request)
    return response


def read_mols(mol_file: Union[str, Path]) -> list[Chem.Mol]:
    mol_file = str(mol_file)
    suppl = None
    if mol_file.endswith('.smi'):
        suppl = Chem.SmilesMolSupplier(mol_file, titleLine=False)
    elif mol_file.endswith(('.csv')):
        suppl = Chem.SmilesMolSupplier(mol_file, titleLine=True,
                                       delimiter=',', smilesColumn=1,
                                       nameColumn=0)
    elif mol_file.endswith('.sdf'):
        suppl = Chem.SDMolSupplier(mol_file)
    mols = []
    for mol in suppl:
        if mol:
            mols.append(mol)
    return mols


class ScriptTest(TestCase):

    def setUp(self) -> None:
        """
        It's important when doing production testing to use the script
        that's in the YAML, which is in the DataFxn repo, as that's the
        one that Spotfire will be running.  When developing, it's
        convenient to edit the script in this directory, transferring
        it to the YAML file once it's ready for production.
        If there's a script PAINSFilters_script_dev.py,
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
            from PAINSFilters_script_dev import (execute, highlight_molecules,
                                                 highlight_molecule,
                                                 run_pains)
            self._script_file = None
        else:
            data_fxn_dir = this_dir.parent.parent / 'DataFxns'
            print(f'Using production script in {data_fxn_dir}')
            yaml_file = data_fxn_dir / 'python' / 'local' / 'PAINSFilters.yaml'
            self._script_file = this_dir / 'PAINSFilters_script.py'
            if yaml_file.exists():
                script_lines = extract_script(yaml_file)
                with open(self._script_file, 'w') as f:
                    f.write(''.join(script_lines))

            from PAINSFilters_script import (execute, highlight_molecules,
                                             highlight_molecule,
                                             run_pains)
        self._execute = execute
        self._highlight_molecules = highlight_molecules
        self._highlight_molecule = highlight_molecule
        self._run_pains = run_pains

    def tearDown(self) -> None:
        # if a script file was made, tidy it up
        if self._script_file is not None and self._script_file.exists():
            self._script_file.unlink()
            pass  # so we can umcomment the unlink() easily

    def test_script(self) -> None:
        # print(f'Running test : {inspect.stack()[0].function}', flush=True)
        file_in = Path(__file__).parent / 'resources' / 'test_pains_filters.json'
        response = run_script(file_in, self._execute)
        self.assertTrue(response)
        self.assertEqual(len(response.outputColumns), 4)
        self.assertEqual(len(response.outputColumns[0].values), 51)
        num_true = len([p for p in response.outputColumns[0].values if p])
        self.assertEqual(num_true, 4)
        self.assertEqual(response.outputColumns[1].values[-1], 'azo_A(324)')
        self.assertEqual(response.outputColumns[1].name, 'PAINS NAMES SMILES')
        mols = column_to_molecules(response.outputColumns[3])
        self.assertIsNotNone(mols[-1])
        self.assertEqual(mols[-1].GetProp('Renderer_Highlight'),
                         'COLOR #ff0000\nATOMS 21 22\nBONDS 21')
        self.assertIsNone(mols[0])

    def test_multi_hits(self) -> None:
        # print(f'Running test : {inspect.stack()[0].function}', flush=True)
        mols = read_mols(Path(__file__).parent / 'resources' / 'chembl_multi_pains.csv')
        self.assertEqual(len(mols), 5)
        pains_hits = self._run_pains(mols)
        pains_mols = [p[2] for p in pains_hits]
        self._highlight_molecules(mols, pains_mols)
        # mol 0 has 2 overlapping hits
        self.assertEqual(len(pains_hits[0][1].split(',')), 2, pains_hits[0][1])
        self.assertEqual(mols[0].GetProp('Renderer_Highlight'),
                         'COLOR #ff0000\nATOMS 14 15 17 22 23 24 25 26\nBONDS 14 16 22 23 24 25 27 28')
        # mol 2 has 2 non-overlapping hits
        self.assertEqual(len(pains_hits[2][1].split(',')), 2)
        self.assertEqual(mols[2].GetProp('Renderer_Highlight'),
                         'COLOR #ff0000\nATOMS 1 2 3 4 5 6 7 8 9 19 20 21 22 23 24 25 26 27 28 29 30\nBONDS 1 2 4 5 6 7 8 18 19 20 21 23 24 25 26 27 28 29 30 31 33')
        # mol 3 has 1 hit
        self.assertEqual(len(pains_hits[3][1].split(',')), 1)
        self.assertEqual(mols[3].GetProp('Renderer_Highlight'),
                         'COLOR #ff0000\nATOMS 1 2 3 4 14 15 16 21\nBONDS 1 2 3 13 14 15 21 23')

    def test_results(self) -> None:
        # print(f'Running test : {inspect.stack()[0].function}', flush=True)
        try:
            rdbase = os.environ['RDBASE']
        except KeyError:
            print(f'ERROR : no RDBASE')
            return None

        mols = []
        pains_test_file = Path(rdbase) / 'Data' / 'Pains' / 'test_data' / 'test_set3.txt'
        with open(pains_test_file, 'r') as f:
            for mline in f:
                if mline[0] == '#':
                    continue
                mline_bits = mline.strip().split()
                mols.append(Chem.MolFromSmiles(mline_bits[2]))
                mols[-1].SetProp('_Name', mline_bits[3])
        pains_hits = self._run_pains(mols)
        num_true = len([p[0] for p in pains_hits if p[0]])
        # each of these test SMILES should match a PAINS search
        self.assertEqual(num_true, len(mols))
        num_false = len([p[0] for p in pains_hits if not p[0]])
        self.assertEqual(num_false, 0)
        for mol, hitset in zip(mols, pains_hits):
            # use hitset[2] to highlight mol
            self._highlight_molecule(mol, hitset[2])
            self.assertTrue(mol.HasProp('Renderer_Highlight'))


if __name__ == '__main__':
    main()
