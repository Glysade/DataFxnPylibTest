import glob
import importlib
import os
import shutil
from typing import Callable, Tuple
from unittest import TestCase, main, skip

from df.data_transfer import DataFunctionRequest, DataFunctionResponse, DataFunction

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
        file_in = os.path.join(os.path.dirname(__file__), 'resources', 'AbStructPred_Basic_in.json')
        _, response = run_named_data_function(file_in)

    @classmethod
    def tearDownClass(cls) -> None:
        clean_output_files()


if __name__ == '__main__':
    main()
