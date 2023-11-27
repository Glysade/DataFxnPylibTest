from collections import defaultdict
import glob
import importlib
import os
from typing import Callable, Tuple
from unittest import TestCase, main, skip

from df.data_transfer import DataFunctionRequest, DataFunctionResponse, DataFunction, DataType

def clean_output_files() -> None:
    patterns = ['*.xml', '*.fasta', 'ig*.out', '*_in.dnd', '*_out.aln', '*.tree', '*.pdb', '*.phr', '*.pin',
                '*.pot', '*.psq', '*.ptf', '*.pto']
    for patt in patterns:
        for f in glob.glob(patt):
            os.remove(f)


def request_from_file(in_file: str) -> DataFunctionRequest:
    with open(in_file, 'r') as fh:
        request_json = fh.read()

    # for testing on WSL convert windows paths to Linux
    if os.name == 'posix':
        request_json = request_json.replace(r'C:\\db\\', r'/mnt/c/db/')
        request_json = request_json.replace(r'mmp\\', r'mmp/')
    request = DataFunctionRequest.parse_raw(request_json)
    return request


def run_named_data_function(in_file: str) -> Tuple[DataFunctionRequest, DataFunctionResponse]:
    request = request_from_file(in_file)
    data_function_class_name = request.serviceName
    module = importlib.import_module(f'df.{data_function_class_name}')
    class_ = getattr(module, data_function_class_name)
    df: DataFunction = class_()
    response = df.execute(request)
    return request, response


def run_script(in_file: str, execute: Callable[[DataFunctionRequest], DataFunctionResponse]) -> Tuple[
    DataFunctionRequest, DataFunctionResponse]:
    request = request_from_file(in_file)
    response = execute(request)
    return request, response


class mAb_DataFunctionTest(TestCase):

    def test_antibody_structure_prediction(self) -> None:
        # test simplest case
        file_in = os.path.join(os.path.dirname(__file__), 'resources', 'AbStructPred_Basic_in.json')
        _, response = run_named_data_function(file_in)

        self.assertTrue(response)

        self.assertEqual(1, len(response.outputTables))
        self.assertEqual(8, len(response.outputTables[0].columns[0].values))

        expected_table = defaultdict(list)
        with open(os.path.join(os.path.dirname(__file__), 'resources', 'AbStructPred_Basic_out.txt')) as df_out:
            header = df_out.readline().strip().split('\t')
            for l in df_out:
                l_tokens = l.strip().split('\t')
                for idx, h in enumerate(header):
                    expected_table[h].append(l_tokens[idx])

            for column in response.outputTables[0].columns:
                if column.name in expected_table:
                    expected_values = expected_table[column.name]
                    actual_values = column.values
                    if column.dataType == DataType.INTEGER:
                        expected_values = [int(v) for v in expected_values]
                    self.assertTrue(all([exp == act for exp, act in zip(expected_values, actual_values)]))

        # test most complex case
        file_in = os.path.join(os.path.dirname(__file__), 'resources', 'AbStructPred_AllOptions_in.json')
        _, response = run_named_data_function(file_in)

        self.assertTrue(response)

        self.assertEqual(1, len(response.outputTables))
        self.assertEqual(32, len(response.outputTables[0].columns[0].values))

        expected_table = defaultdict(list)
        with open(os.path.join(os.path.dirname(__file__), 'resources', 'AbStructPred_AllOptions_out.txt')) as df_out:
            header = df_out.readline().strip().split('\t')
            for l in df_out:
                l_tokens = l.strip().split('\t')
                for idx, h in enumerate(header):
                    expected_table[h].append(l_tokens[idx])

            for column in response.outputTables[0].columns:
                if column.dataType != DataType.BINARY and column.name in expected_table:
                    expected_values = expected_table[column.name]
                    actual_values = column.values
                    if column.dataType == DataType.INTEGER:
                        expected_values = [int(v) for v in expected_values]
                    self.assertTrue(all([exp == act for exp, act in zip(expected_values, actual_values)]))

    @classmethod
    def tearDownClass(cls) -> None:
        clean_output_files()


if __name__ == '__main__':
    main()
