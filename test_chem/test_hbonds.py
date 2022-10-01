from pathlib import Path
from typing import Callable
from unittest import TestCase, main

from df.chem_helper import column_to_molecules
from df.data_transfer import DataFunctionRequest, DataFunctionResponse


def run_script(in_file: str, execute: Callable[[DataFunctionRequest], DataFunctionResponse]) -> DataFunctionResponse:
    with open(in_file, 'r') as fh:
        request_json = fh.read()
    
    request = DataFunctionRequest.parse_raw(request_json)  
    response = execute(request)
    return response


class ScriptTest(TestCase):
    
    def test_script(self) -> None:
        from Hbonds_script import execute
        file_in = Path(__file__).parent / 'resources' / 'test_hbonds.json'
        response = run_script(file_in, execute)
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
                         'COLOR #dc143c\nATOMS 12\nBONDS\n'
                         'COLOR #00bfff\nATOMS 5\nBONDS\n'
                         'COLOR #ffbf00\nATOMS 4 8\nBONDS')


if __name__ == '__main__':
    main()
