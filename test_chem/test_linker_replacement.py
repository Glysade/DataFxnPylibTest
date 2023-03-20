
# Without adding the data function library to the project content roots PyCharm will not run this
# test unless the data function library is added to sys.paths
# import sys
# sys.path.append(r'C:\Users\david\Desktop\Spotfire12\Modules\Glysade Python_4.2.0.6\Python\Lib\site-packages')
# sys.path.append(r'C:\Users\david\AppData\Local\Temp\ScriptSync\DataFxnPylib\pylib_3.9.7')

# To enable the test to run within PyCharm without explicitly modifying sys.path
# Add C:\Users\david\AppData\Local\Temp\ScriptSync\DataFxnPylib\pylib_3.9.7 as a Project content root and comment the lines above

from pathlib import Path
from unittest import TestCase, main

from df.LinkerReplacements import LinkerReplacements
from df.data_transfer import (DataFunction, DataFunctionRequest,
                              DataFunctionResponse)


def run_script(in_file: str, test_class: DataFunction) -> DataFunctionResponse:
    with open(in_file, 'r') as fh:
        request_json = fh.read()

    request = DataFunctionRequest.parse_raw(request_json)
    response = test_class.execute(request)
    return response


class ScriptTest(TestCase):
    
    def test_script(self) -> None:
        file_in = Path(__file__).parent / 'resources' / 'test_linker_replacement1.json'
        lr = LinkerReplacements()
        response = run_script(file_in, lr)
        self.assertTrue(response)


if __name__ == '__main__':
    main()
