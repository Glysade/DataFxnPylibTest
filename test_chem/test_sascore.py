from pathlib import Path
from typing import Callable
from unittest import TestCase, main

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
        If there's a script SyntheticAccessibilityScore_script_dev.py,
        use that, otherwise make one from the YAML file and use that.
        Arbitrary_script_dev.py should not be in the
        repo, only ever in a development branch so testing a production
        repo will use the script derived from the YAML.
        Doing it this way allows for people editing the YAML directly,
        which is handy for minor tweaks, and those tweaks always
        being tested.
        Assume the DataFxn repo has been cloned alongside this one.
        """
        global PRINTED_GUFF
        this_dir = Path(__file__).parent
        script_file = this_dir / 'SyntheticAccessibilityScore_script_dev.py'
        if Path(script_file).exists():
            if not PRINTED_GUFF:
                print(f'Using development script')
                PRINTED_GUFF = True
            from SyntheticAccessibilityScore_script_dev import execute
            self._script_file = None
        else:
            data_fxn_dir = this_dir.parent.parent / 'DataFxns'
            yaml_file = data_fxn_dir / 'python' / 'local' / 'SyntheticAccessibilityScore.yaml'
            if not PRINTED_GUFF:
                print(f'Using production script in {yaml_file}')
                PRINTED_GUFF = True
            self._script_file = this_dir / 'QEDScore_script.py'
            if yaml_file.exists():
                script_lines = extract_script(yaml_file)
                with open(self._script_file, 'w') as f:
                    f.write(''.join(script_lines))

            from SyntheticAccessibilityScore_script import execute
        self._execute = execute

    def tearDown(self) -> None:
        # if a script file was made, tidy it up
        if self._script_file is not None and self._script_file.exists():
            self._script_file.unlink()
            pass  # so we can comment the unlink() easily

    def test_script(self) -> None:
        file_in = Path(__file__).parent / 'resources' / 'test_sascore1.json'
        response = run_script(file_in, self._execute)
        scores = response.outputColumns[0].values
        self.assertEqual(15, len(scores))
        self.assertAlmostEqual(3.4, scores[0], 1)
        self.assertAlmostEqual(3.75, scores[-1], 1)
        self.assertTrue(response)


if __name__ == '__main__':
    main()
