from pathlib import Path
from typing import Callable
from unittest import TestCase, main

from df.chem_helper import column_to_molecules
from df.data_transfer import DataFunctionRequest, DataFunctionResponse
from rdkit import Chem
from rdkit.Chem import AllChem

from ArbitraryReaction_script import run_reactions


def run_script(in_file: str, execute: Callable[[DataFunctionRequest], DataFunctionResponse]) -> DataFunctionResponse:
    with open(in_file, 'r') as fh:
        request_json = fh.read()

    request = DataFunctionRequest.parse_raw(request_json)
    response = execute(request)
    return response


class ScriptTest(TestCase):

    def test_script(self) -> None:
        from ArbitraryReaction_script import execute
        file_in = Path(__file__).parent / 'resources' / 'test_arbitrary_reaction.json'
        response = run_script(file_in, execute)
        self.assertTrue(response)
        self.assertEqual(len(response.outputColumns), 1)
        prods = column_to_molecules(response.outputColumns[0])
        self.assertEqual(len(prods), 83)
        self.assertEqual(prods.count(None), 53)
        self.assertEqual(Chem.MolToSmiles(prods[2]), 'NC(=O)c1cccc(N)c1')

    def test_script1(self) -> None:
        from ArbitraryReaction_script import execute
        file_in = Path(__file__).parent / 'resources' / 'test_arbitrary_reaction1.json'
        response = run_script(file_in, execute)
        self.assertTrue(response)
        self.assertEqual(len(response.outputColumns), 1)
        prods = column_to_molecules(response.outputColumns[0])
        self.assertEqual(len(prods), 2)
        self.assertEqual(prods.count(None), 1)
        self.assertEqual(Chem.MolToSmiles(prods[0]), 'CCC=O.O[C@@H](F)Cl')
        self.assertIsNone(prods[1])

    def test_rxn_file(self) -> None:
        rxn_file = Path(__file__).parent / 'resources' / 'amide_2_amine.rxn'
        rxn = AllChem.ReactionFromRxnFile(str(rxn_file))
        self.assertIsNotNone(rxn)
        mols = [Chem.MolFromSmiles('c1ccccc1C(=O)NC'),
                Chem.MolFromSmiles('C1CCCCC1C(=O)NC'),
                Chem.MolFromSmiles('CCNC(=O)c1ccc(C(=O)NC)cc1')]
        prods = run_reactions(mols, rxn)
        self.assertEqual(len(prods), 3)
        self.assertEqual(Chem.MolToSmiles(prods[0]), 'NC(=O)c1ccccc1')
        self.assertIsNone(prods[1])
        self.assertEqual(Chem.MolToSmiles(prods[2]), 'NC(=O)c1ccc(C(N)=O)cc1')

    def test_two_products(self) -> None:
        rxn_file = Path(__file__).parent / 'resources' / 'amide_split.rxn'
        rxn = AllChem.ReactionFromRxnFile(str(rxn_file))
        self.assertIsNotNone(rxn)
        smis = ['Nc1cccc(c1)C(=O)Nc1ccc(N)cc1',
                'c1ccccc1C(Cl)C(Cl)C(=O)Nc1ccccc1']
        mols = [Chem.MolFromSmiles(s) for s in smis]
        prods = run_reactions(mols, rxn)
        self.assertEqual(len(prods), 2)
        self.assertEqual(Chem.MolToSmiles(prods[0]),
                         'Nc1ccc(N)cc1.Nc1cccc(C=O)c1')
        self.assertIsNone(prods[1])

    def test_rxn_smarts(self) -> None:
        # this is the same as test_script1, but done directly
        rxn_smt = '[C:0][C:1](=[O:2])[O:3][C@:4]>>[C:0][C:1]=[O:2].[O:3][C@@:4]'
        rxn = AllChem.ReactionFromSmarts(rxn_smt)
        smis = ['CCC(=O)O[C@H](F)Cl',
                'c1ccccc1C(=O)O[C@@](F)Cl']
        mols = [Chem.MolFromSmiles(s) for s in smis]
        prods = run_reactions(mols, rxn)
        self.assertEqual(Chem.MolToSmiles(prods[0]), 'CCC=O.O[C@@H](F)Cl')
        self.assertIsNone(prods[1])


if __name__ == '__main__':
    main()
