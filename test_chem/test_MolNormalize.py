from json import dumps, loads
from pathlib import Path
from typing import Callable
from unittest import TestCase, main

from df.chem_helper import column_to_molecules
from df.data_transfer import DataFunctionRequest, DataFunctionResponse
from yaml_to_script import extract_script
from rdkit import Chem

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

PRINTED_GUFF = False

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

    def setUp(self) -> None:
        """
        It's important when doing production testing to use the script
        that's in the YAML, which is in the DataFxn repo, as that's the
        one that Spotfire will be running.  When developing, it's
        convenient to edit the script in this directory, transferring
        it to the YAML file once it's ready for production.
        If there's a script MolNormalize_script_dev.py,
        use that, otherwise make one from the YAML file and use that.
        Arbitrary_script_dev.py should not be in the
        repo, only ever in a development branch so testing a production
        repo will use the script derived from the YAML.
        Doing it this way allows for people editing the YAML directly,
        which is handing for minor tweaks, and those tweaks always
        being tested.
        Assume the DataFxn repo has been cloned alongside this one.
        """
        global PRINTED_GUFF
        this_dir = Path(__file__).parent
        script_file = this_dir / 'MolNormalize_script_dev.py'
        if Path(script_file).exists():
            if not PRINTED_GUFF:
                print(f'Using development script')
                PRINTED_GUFF = True
            from MolNormalize_script_dev import execute, run_reactions
            self._script_file = None
        else:
            data_fxn_dir = this_dir.parent.parent / 'DataFxns'
            yaml_file = data_fxn_dir / 'python' / 'local' / 'MolNormalize.yaml'
            if not PRINTED_GUFF:
                print(f'Using production script in {yaml_file}')
                PRINTED_GUFF = True
            self._script_file = this_dir / 'MolNormalize_script.py'
            if yaml_file.exists():
                script_lines = extract_script(yaml_file)
                with open(self._script_file, 'w') as f:
                    f.write(''.join(script_lines))

            from MolNormalize_script import execute
        self._execute = execute

    def tearDown(self) -> None:
        # if a script file was made, tidy it up
        if self._script_file is not None and self._script_file.exists():
            self._script_file.unlink()
            pass  # so we can comment the unlink() easily

    def test_full_norms(self) -> None:
        file_in = Path(__file__).parent / 'resources' / 'test_mol_normalize1.json'
        response = run_script(file_in, self._execute)
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

        response = run_script_json(dumps(json_dict), self._execute)
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

        response = run_script_json(dumps(json_dict), self._execute)
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

        response = run_script_json(dumps(json_dict), self._execute)
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

        response = run_script_json(dumps(json_dict), self._execute)
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
