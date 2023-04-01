from pathlib import Path
from typing import Callable
from unittest import TestCase, main

from df.chem_helper import column_to_molecules
from df.data_transfer import DataFunctionRequest, DataFunctionResponse
from yaml_to_script import extract_script

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
        If there's a script Hbonds_script_dev.py,
        use that, otherwise make one from the YAML file and use that.
        Hbonds_script_dev.py should not be in the
        repo, only ever in a development branch so testing a production
        repo will use the script derived from the YAML.
        Doing it this way allows for people editing the YAML directly,
        which is handing for minor tweaks, and those tweaks always
        being tested.
        Assume the DataFxn repo has been cloned alongside this one.
        """
        this_dir = Path(__file__).parent
        script_file = this_dir / 'Hbonds_script_dev.py'
        global PRINTED_GUFF
        if Path(script_file).exists():
            if not PRINTED_GUFF:
                print(f'Using development script')
                PRINTED_GUFF = True
            from Hbonds_script_dev import execute
            self._script_file = None
        else:
            data_fxn_dir = this_dir.parent.parent / 'DataFxns'
            yaml_file = data_fxn_dir / 'python' / 'local' / 'HBondSMARTS.yaml'
            if not PRINTED_GUFF:
                print(f'Using production script in {yaml_file}')
                PRINTED_GUFF = True
            self._script_file = this_dir / 'Hbonds_script.py'
            if yaml_file.exists():
                script_lines = extract_script(yaml_file)
                with open(self._script_file, 'w') as f:
                    f.write(''.join(script_lines))

            from Hbonds_script import execute
        self._execute = execute

    def tearDown(self) -> None:
        # if a script file was made, tidy it up
        if self._script_file is not None and self._script_file.exists():
            self._script_file.unlink()
            pass  # so we can comment the unlink() easily

    def test_script(self) -> None:
        file_in = Path(__file__).parent / 'resources' / 'test_hbonds.json'
        response = run_script(file_in, self._execute)
        self.assertTrue(response)
        self.assertEqual(len(response.outputColumns), 3)
        self.assertEqual(len(response.outputColumns[0].values), 53)
        don_counts = response.outputColumns[1].values
        acc_counts = response.outputColumns[2].values
        self.assertEqual(don_counts[0], 2)
        self.assertEqual(acc_counts[0], 2)
        self.assertEqual(don_counts[-2], 0)
        self.assertEqual(acc_counts[-2], 0)
        self.assertEqual(don_counts[-1], 3)
        self.assertEqual(acc_counts[-1], 3)
        mols = column_to_molecules(response.outputColumns[0])
        # the last molecule has 3 donors, 3 acceptors, but 2 of each
        # are both, so coloured orange.
        self.assertEqual(mols[-1].GetProp('Renderer_Highlight'),
                         'COLOR #b3b3ff\nATOMS 12\nBONDS\n'
                         'COLOR #ff9494\nATOMS 5\nBONDS\n'
                         'COLOR #d5b3ef\nATOMS 4 8\nBONDS')


if __name__ == '__main__':
    main()
