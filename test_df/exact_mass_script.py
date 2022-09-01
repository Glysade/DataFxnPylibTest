
# The contents of the data function definition script appear below
# the code should define an execute method

# Once this is working copy the code below to the script tag in the data function definition

from rdkit.Chem.Descriptors import ExactMolWt
from df.chem_helper import column_to_molecules
from df.data_transfer import DataFunction, DataFunctionRequest, DataFunctionResponse, DataType, ColumnData, \
    string_input_field


def execute(request: DataFunctionRequest) -> DataFunctionResponse:
    column_id = string_input_field(request, 'structureColumn')
    input_column = request.inputColumns[column_id]
    mols = column_to_molecules(input_column)
    weights = [None if m is None else ExactMolWt(m) for m in mols]
    output_column = ColumnData(name=f'{input_column.name} Exact Mass', dataType=DataType.DOUBLE, values=weights)
    response = DataFunctionResponse(outputColumns=[output_column])
    return response

