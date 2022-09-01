import importlib
import os
import uuid
from typing import List
from unittest import TestCase, main

from df.chem_helper import molecules_to_column
from df.data_transfer import DataType, DataFunctionResponse, DataFunctionRequest, ColumnData, InputField, DataFunction
from ruse.rdkit.rdkit_utils import file_to_mols


def run_data_function_module(request: DataFunctionRequest) -> DataFunctionResponse:
    class_name = request.serviceName
    module = importlib.import_module(f'df.{class_name}')
    class_ = getattr(module, class_name)
    df: DataFunction = class_()
    response = df.execute(request)
    return response


class ChemClientTest(TestCase):

    def _check_columns_present(self, number_values: int, number_columns: int, output_columns: List[ColumnData]) -> None:
        self.assertEqual(number_columns, len(output_columns))
        for column_number in range(number_columns):
            self.assertEqual(number_values, len(output_columns[column_number].values))

    def test_exact_mass(self):
        mol_file = os.path.join(os.path.dirname(__file__), '../test_chem/resources/Enamine_1000.smi')
        mols = file_to_mols(mol_file)
        column_input_field = InputField(id='structureColumn', data='columnId',
                                        dataType=DataType.STRING)
        input_column = molecules_to_column(mols, 'Structures', DataType.STRING)
        request = DataFunctionRequest(serviceName='ExactMass', inputColumns={'columnId': input_column},
                                      inputFields={'structureColumn': column_input_field},
                                      id=str(uuid.uuid4()))
        response = run_data_function_module(request)
        self._check_columns_present(len(mols), 1, response.outputColumns)


if __name__ == "__main__":
    main()
