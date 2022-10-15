from json import dumps, loads
from pathlib import Path
from typing import Callable
from unittest import TestCase, main

from df.chem_helper import column_to_molecules
from df.data_transfer import DataFunctionRequest, DataFunctionResponse
from rdkit import Chem

from MolNormalize_script import execute


def run_script(in_file: str, execute: Callable[[DataFunctionRequest], DataFunctionResponse]) -> DataFunctionResponse:
    with open(in_file, 'r') as fh:
        request_json = fh.read()

    return run_script_json(request_json, execute)


def run_script_json(in_json: dict,
                    execute: Callable[[DataFunctionRequest], DataFunctionResponse]) -> DataFunctionResponse:
    request = DataFunctionRequest.parse_raw(in_json)
    response = execute(request)
    return response


class ScriptTest(TestCase):

    def test_full_norms(self) -> None:
        file_in = Path(__file__).parent / 'resources' / 'test_mol_normalize1.json'
        response = run_script(file_in, execute)
        self.assertTrue(response)
        mols = column_to_molecules(response.outputColumns[0])
        self.assertEqual(len(mols), 26)
        full_norms = ["CCC(=O)F", "O=c1c2ccccc2[nH]c2ccncc12", "CCS(C)(=O)=O",
                      "CC[S@@+](C)[O-]", "CCP(C)(=N)N", "CN=[N+]=[N-]",
                      "C=[N+]=[N-]", "CC#N", "CN(C)c1ccccn1", "CN(C)C=CC=N",
                      "O=C(O)c1ccccc1", "O=C([O-])c1ccccc1", "CCN(C)C",
                      "CC[N+](C)(C)C", "CN(C)CC(CO)C(=O)O",
                      "C[N+](C)(C)CC(CO)C(=O)[O-]", "O=C(O)c1ccccc1",
                      "O=C(O)c1cccc(O)c1", "O=C([O-])c1ccccc1",
                      "O=C1C(CC[C@H](O)c2ccc(F)cc2)[C@@H](c2ccc(O)cc2)N1c1ccc(F)cc1",
                      "CC[O-]", "CC[O-]", "C[Hg]C", "O=C([O-])CCc1ccccc1[Mg]Br",
                      "COc1ccc2[nH]c([S@+]([O-])Cc3ncc(C)c(OC)c3C)nc2c1",
                      "COc1ccc2[nH]c([S+]([O-])Cc3ncc(C)c(OC)c3C)nc2c1"]
        mol_smis = [Chem.MolToSmiles(m) for m in mols]
        self.assertListEqual(mol_smis, full_norms)

    def test_normalize(self) -> None:
        self.maxDiff = None
        file_in = Path(__file__).parent / 'resources' / 'test_mol_normalize1.json'
        with open(file_in, 'r') as fh:
            request_json = fh.read()
        json_dict = loads(request_json)
        json_dict['inputFields']['normalizeID']['data'] = True
        json_dict['inputFields']['neutralizeID']['data'] = False
        json_dict['inputFields']['reionizeID']['data'] = False
        json_dict['inputFields']['parentID']['data'] = False
        json_dict['inputFields']['metalID']['data'] = False
        json_dict['inputFields']['tautomerID']['data'] = False

        response = run_script_json(dumps(json_dict), execute)
        self.assertTrue(response)
        mols = column_to_molecules(response.outputColumns[0])
        self.assertEqual(len(mols), 26)
        full_norms = ["CC=C(O)F", "Oc1c2ccccc2nc2ccncc12", "CCS(C)(=O)=O",
                      "CC[S@@+](C)[O-]", "CCP(C)(N)=[NH2+]", "CN=[N+]=[N-]",
                      "C=[N+]=[N-]", "CC#N", "CN(C)c1cccc[nH+]1", "CN(C)C=CC=[NH2+]",
                      "O=C([O-])c1ccccc1", "O=C([O-])c1ccccc1.[Na+]", "CC[NH+](C)C",
                      "CC[N+](C)(C)C", "C[NH+](C)CC(C[O-])C(=O)[O-]",
                      "C[N+](C)(C)CC(C[O-])C(=O)[O-]", "O=C([O-])c1ccccc1",
                      "O=C(O)c1cccc([O-])c1", "O=C([O-])c1ccccc1.[Na+]",
                      "O=C1[C@H](CC[C@H](O)c2ccc(F)cc2)[C@@H](c2ccc(O)cc2)N1c1ccc(F)cc1.c1ccccc1",
                      "CCO[Fe]", "CCO[AlH2]", "C[Hg]C", "O=C(CCc1ccccc1[Mg]Br)O[Na]",
                      "COc1ccc2nc([S@+]([O-])Cc3ncc(C)c(OC)c3C)[nH]c2c1",
                      "COc1ccc2[n-]c([S+]([O-])Cc3ncc(C)c(OC)c3C)nc2c1.COc1ccc2[n-]c([S+]([O-])Cc3ncc(C)c(OC)c3C)nc2c1.[Mg+2]"]
        mol_smis = [Chem.MolToSmiles(m) for m in mols]
        self.assertListEqual(mol_smis, full_norms)

    def test_plus_neutralize(self) -> None:
        self.maxDiff = None
        file_in = Path(__file__).parent / 'resources' / 'test_mol_normalize1.json'
        with open(file_in, 'r') as fh:
            request_json = fh.read()
        json_dict = loads(request_json)
        json_dict['inputFields']['normalizeID']['data'] = True
        json_dict['inputFields']['neutralizeID']['data'] = True
        json_dict['inputFields']['reionizeID']['data'] = False
        json_dict['inputFields']['parentID']['data'] = False
        json_dict['inputFields']['metalID']['data'] = False
        json_dict['inputFields']['tautomerID']['data'] = False
        # json_dict['inputColumns']['43f9cc7b-82d2-426a-ba04-cfb9daa6664fsSMILES']['values'] = json_dict['inputColumns']['43f9cc7b-82d2-426a-ba04-cfb9daa6664fsSMILES']['values'][8:9]
        # print(json_dict['inputColumns']['43f9cc7b-82d2-426a-ba04-cfb9daa6664fsSMILES']['values'])

        response = run_script_json(dumps(json_dict), execute)
        self.assertTrue(response)
        mols = column_to_molecules(response.outputColumns[0])
        self.assertEqual(len(mols), 26)
        full_norms = ["CC=C(O)F", "Oc1c2ccccc2nc2ccncc12", "CCS(C)(=O)=O",
                      "CC[S@@+](C)[O-]", "CCP(C)(=N)N", "CN=[N+]=[N-]",
                      "C=[N+]=[N-]", "CC#N", "CN(C)c1ccccn1", "CN(C)C=CC=N",
                      "O=C(O)c1ccccc1", "O=C([O-])c1ccccc1.[Na+]", "CCN(C)C",
                      "CC[N+](C)(C)C", "CN(C)CC(CO)C(=O)O",
                      "C[N+](C)(C)CC(CO)C(=O)[O-]", "O=C(O)c1ccccc1",
                      "O=C(O)c1cccc(O)c1", "O=C([O-])c1ccccc1.[Na+]",
                      "O=C1[C@H](CC[C@H](O)c2ccc(F)cc2)[C@@H](c2ccc(O)cc2)N1c1ccc(F)cc1.c1ccccc1",
                      "CCO[Fe]", "CCO[AlH2]", "C[Hg]C", "O=C(CCc1ccccc1[Mg]Br)O[Na]",
                      "COc1ccc2nc([S@+]([O-])Cc3ncc(C)c(OC)c3C)[nH]c2c1",
                      "COc1ccc2[nH]c([S+]([O-])Cc3ncc(C)c(OC)c3C)nc2c1.COc1ccc2[nH]c([S+]([O-])Cc3ncc(C)c(OC)c3C)nc2c1.[Mg+2]"]
        mol_smis = [Chem.MolToSmiles(m) for m in mols]
        self.assertListEqual(mol_smis, full_norms)

    def test_plus_reionize(self) -> None:
        self.maxDiff = None
        file_in = Path(__file__).parent / 'resources' / 'test_mol_normalize1.json'
        with open(file_in, 'r') as fh:
            request_json = fh.read()
        json_dict = loads(request_json)
        json_dict['inputFields']['normalizeID']['data'] = True
        json_dict['inputFields']['neutralizeID']['data'] = True
        json_dict['inputFields']['reionizeID']['data'] = True
        json_dict['inputFields']['parentID']['data'] = False
        json_dict['inputFields']['metalID']['data'] = False
        json_dict['inputFields']['tautomerID']['data'] = False

        response = run_script_json(dumps(json_dict), execute)
        self.assertTrue(response)
        self.assertEqual(len(response.outputColumns), 1)
        mols = column_to_molecules(response.outputColumns[0])
        self.assertEqual(len(mols), 26)
        full_norms = ["CC=C(O)F", "Oc1c2ccccc2nc2ccncc12", "CCS(C)(=O)=O",
                      "CC[S@@+](C)[O-]", "CCP(C)(=N)N", "CN=[N+]=[N-]",
                      "C=[N+]=[N-]", "CC#N", "CN(C)c1ccccn1", "CN(C)C=CC=N",
                      "O=C(O)c1ccccc1", "O=C([O-])c1ccccc1.[Na+]", "CCN(C)C",
                      "CC[N+](C)(C)C", "CN(C)CC(CO)C(=O)O",
                      "C[N+](C)(C)CC(CO)C(=O)[O-]", "O=C(O)c1ccccc1",
                      "O=C(O)c1cccc(O)c1", "O=C([O-])c1ccccc1.[Na+]",
                      "O=C1[C@H](CC[C@H](O)c2ccc(F)cc2)[C@@H](c2ccc(O)cc2)N1c1ccc(F)cc1.c1ccccc1",
                      "CCO[Fe]", "CCO[AlH2]", "C[Hg]C", "O=C(CCc1ccccc1[Mg]Br)O[Na]",
                      "COc1ccc2nc([S@+]([O-])Cc3ncc(C)c(OC)c3C)[nH]c2c1",
                      "COc1ccc2[nH]c([S+]([O-])Cc3ncc(C)c(OC)c3C)nc2c1.COc1ccc2[nH]c([S+]([O-])Cc3ncc(C)c(OC)c3C)nc2c1.[Mg+2]"]
        mol_smis = [Chem.MolToSmiles(m) for m in mols]
        self.assertListEqual(mol_smis, full_norms)

    def test_plus_parent(self) -> None:
        self.maxDiff = None
        file_in = Path(__file__).parent / 'resources' / 'test_mol_normalize1.json'
        with open(file_in, 'r') as fh:
            request_json = fh.read()
        json_dict = loads(request_json)
        json_dict['inputFields']['normalizeID']['data'] = True
        json_dict['inputFields']['neutralizeID']['data'] = True
        json_dict['inputFields']['reionizeID']['data'] = True
        json_dict['inputFields']['parentID']['data'] = True
        json_dict['inputFields']['metalID']['data'] = False
        json_dict['inputFields']['tautomerID']['data'] = False

        response = run_script_json(dumps(json_dict), execute)
        self.assertTrue(response)
        mols = column_to_molecules(response.outputColumns[0])
        self.assertEqual(len(mols), 26)
        full_norms = ["CC=C(O)F", "Oc1c2ccccc2nc2ccncc12", "CCS(C)(=O)=O",
                      "CC[S@@+](C)[O-]", "CCP(C)(=N)N", "CN=[N+]=[N-]",
                      "C=[N+]=[N-]", "CC#N", "CN(C)c1ccccn1", "CN(C)C=CC=N",
                      "O=C(O)c1ccccc1", "O=C([O-])c1ccccc1", "CCN(C)C",
                      "CC[N+](C)(C)C", "CN(C)CC(CO)C(=O)O",
                      "C[N+](C)(C)CC(CO)C(=O)[O-]", "O=C(O)c1ccccc1",
                      "O=C(O)c1cccc(O)c1", "O=C([O-])c1ccccc1",
                      "O=C1[C@H](CC[C@H](O)c2ccc(F)cc2)[C@@H](c2ccc(O)cc2)N1c1ccc(F)cc1",
                      "CC[O-]", "CC[O-]", "C[Hg]C", "O=C([O-])CCc1ccccc1[Mg]Br",
                      "COc1ccc2nc([S@+]([O-])Cc3ncc(C)c(OC)c3C)[nH]c2c1",
                      "COc1ccc2[nH]c([S+]([O-])Cc3ncc(C)c(OC)c3C)nc2c1"]
        mol_smis = [Chem.MolToSmiles(m) for m in mols]
        self.assertListEqual(mol_smis, full_norms)


if __name__ == '__main__':
    main()
