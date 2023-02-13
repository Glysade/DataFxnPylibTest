from pathlib import Path
from typing import Callable
from unittest import TestCase, main

from df.chem_helper import column_to_molecules
from df.data_transfer import DataFunctionRequest, DataFunctionResponse
from yaml_to_script import extract_script

from rdkit import Chem
from rdkit.Chem import AllChem

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

PRINTED_GUFF = False

def run_script(in_file: str, execute: Callable[[DataFunctionRequest], DataFunctionResponse]) -> DataFunctionResponse:
    with open(in_file, 'r') as fh:
        request_json = fh.read()

    request = DataFunctionRequest.parse_raw(request_json)
    response = execute(request)
    return response


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
        global PRINTED_GUFF
        this_dir = Path(__file__).parent
        script_file = this_dir / 'MaximumCommonSubstructure_script_dev.py'
        if Path(script_file).exists():
            if not PRINTED_GUFF:
                print(f'Using development script')
                PRINTED_GUFF = True
            from ArbitraryReaction_script_dev import execute, run_reactions
            self._script_file = None
        else:
            data_fxn_dir = this_dir.parent.parent / 'DataFxns'
            yaml_file = data_fxn_dir / 'python' / 'local' / 'ArbitraryReaction.yaml'
            if not PRINTED_GUFF:
                print(f'Using production script in {yaml_file}')
                PRINTED_GUFF = True
            self._script_file = this_dir / 'ArbitraryReaction_script.py'
            if yaml_file.exists():
                script_lines = extract_script(yaml_file)
                with open(self._script_file, 'w') as f:
                    f.write(''.join(script_lines))

            from ArbitraryReaction_script import execute, run_reactions
        self._execute = execute
        self._run_reactions = run_reactions

    def tearDown(self) -> None:
        # if a script file was made, tidy it up
        if self._script_file is not None and self._script_file.exists():
            self._script_file.unlink()
            pass  # so we can umcomment the unlink() easily

    def test_script(self) -> None:
        file_in = Path(__file__).parent / 'resources' / 'test_arbitrary_reaction.json'
        response = run_script(file_in, self._execute)
        self.assertTrue(response)
        self.assertEqual(len(response.outputColumns), 1)
        prods = column_to_molecules(response.outputColumns[0])
        self.assertEqual(len(prods), 83)
        self.assertEqual(prods.count(None), 53)
        self.assertEqual(Chem.MolToSmiles(prods[2]), 'NC(=O)c1cccc(N)c1')

    def test_script1(self) -> None:
        # SMARTS query
        file_in = Path(__file__).parent / 'resources' / 'test_arbitrary_reaction1.json'
        response = run_script(file_in, self._execute)
        self.assertTrue(response)
        self.assertEqual(len(response.outputColumns), 1)
        prods = column_to_molecules(response.outputColumns[0])
        self.assertEqual(len(prods), 2)
        self.assertEqual(prods.count(None), 1)
        self.assertEqual(Chem.MolToSmiles(prods[0]), 'CCC=O.O[C@@H](F)Cl')
        self.assertIsNone(prods[1])

    def test_bad_phenols(self) -> None:
        # MolFile query
        file_in = Path(__file__).parent / 'resources' / 'test_arbitrary_reaction2.json'
        response = run_script(file_in, self._execute)
        self.assertTrue(response)
        self.assertEqual(len(response.outputColumns), 1)
        self.assertEqual(len(response.outputColumns[0].values), 1)
        self.assertIsNone(response.outputColumns[0].values[0])

    def test_rxn_file(self) -> None:
        rxn_file = Path(__file__).parent / 'resources' / 'amide_2_amine.rxn'
        rxn = AllChem.ReactionFromRxnFile(str(rxn_file))
        self.assertIsNotNone(rxn)
        mols = [Chem.MolFromSmiles('c1ccccc1C(=O)NC'),
                Chem.MolFromSmiles('C1CCCCC1C(=O)NC'),
                Chem.MolFromSmiles('CCNC(=O)c1ccc(C(=O)NC)cc1')]
        prods = self._run_reactions(mols, rxn, 'exhaustiveReaction')
        self.assertEqual(len(prods), 3)
        self.assertEqual(Chem.MolToSmiles(prods[0][1]), 'NC(=O)c1ccccc1')
        self.assertIsNone(prods[1][1])
        self.assertEqual(Chem.MolToSmiles(prods[2][1]), 'NC(=O)c1ccc(C(N)=O)cc1')

    def test_two_products(self) -> None:
        rxn_file = Path(__file__).parent / 'resources' / 'amide_split.rxn'
        rxn = AllChem.ReactionFromRxnFile(str(rxn_file))
        self.assertIsNotNone(rxn)
        smis = ['Nc1cccc(c1)C(=O)Nc1ccc(N)cc1',
                'c1ccccc1C(Cl)C(Cl)C(=O)Nc1ccccc1']
        mols = [Chem.MolFromSmiles(s) for s in smis]
        prods = self._run_reactions(mols, rxn, 'exhaustiveReaction')
        self.assertEqual(len(prods), 2)
        self.assertEqual(Chem.MolToSmiles(prods[0][1]),
                         'Nc1ccc(N)cc1.Nc1cccc(C=O)c1')
        self.assertIsNone(prods[1][1])

    def test_rxn_smarts(self) -> None:
        # this is the same as test_script1, but done directly
        rxn_smt = '[C:0][C:1](=[O:2])[O:3][C@:4]>>[C:0][C:1]=[O:2].[O:3][C@@:4]'
        rxn = AllChem.ReactionFromSmarts(rxn_smt)
        smis = ['CCC(=O)O[C@H](F)Cl',
                'c1ccccc1C(=O)O[C@@](F)Cl']
        mols = [Chem.MolFromSmiles(s) for s in smis]
        prods = self._run_reactions(mols, rxn, 'exhaustiveReaction')
        self.assertEqual(Chem.MolToSmiles(prods[0][1]), 'CCC=O.O[C@@H](F)Cl')
        self.assertIsNone(prods[1][1])

    def test_single_mode(self) -> None:
        file_in = Path(__file__).parent / 'resources' / 'test_arbitrary_reaction3.json'
        response = run_script(file_in, self._execute)
        self.assertTrue(response)
        prods = column_to_molecules(response.outputColumns[0])
        self.assertEqual(len(prods), 1)
        self.assertEqual(Chem.MolToSmiles(prods[0]), 'CN(C)c1cnc2cccc(O)c2c1')

    def test_exhaustive_mode(self) -> None:
        file_in = Path(__file__).parent / 'resources' / 'test_arbitrary_reaction4.json'
        response = run_script(file_in, self._execute)
        self.assertTrue(response)
        prods = column_to_molecules(response.outputColumns[0])
        self.assertEqual(len(prods), 1)
        self.assertEqual(Chem.MolToSmiles(prods[0]), 'CN(C)c1cnc2c(O)c(O)c(O)c(O)c2c1')

    def test_multiple_mode(self) -> None:
        file_in = Path(__file__).parent / 'resources' / 'test_arbitrary_reaction5.json'
        response = run_script(file_in, self._execute)
        self.assertTrue(response)
        prods = column_to_molecules(response.outputTables[0].columns[1])
        self.assertEqual(Chem.MolToSmiles(prods[0]), 'CN(C)c1cnc2cccc(O)c2c1')
        self.assertEqual(Chem.MolToSmiles(prods[10]), 'CN(C)c1cnc2cc(O)c(O)c(O)c2c1')
        self.assertEqual(Chem.MolToSmiles(prods[14]), 'CN(C)c1cnc2c(O)c(O)c(O)c(O)c2c1')

if __name__ == '__main__':
    main()
