
# The contents of the data function definition script appear below
# the code should define an execute method

# Once this is working copy the code below to the script tag in the data function definition

from typing import Optional

from df.chem_helper import column_to_molecules, molecules_to_column
from df.data_transfer import DataFunctionRequest, DataFunctionResponse, DataType, ColumnData, \
    string_input_field
from rdkit import Chem
from rdkit.Chem.rdDeprotect import Deprotect
from rdkit.Chem.rdchem import Mol


def compare_molecules(mol1: Optional[Mol], mol2: Optional[Mol]) -> bool:
    if not mol1 or not mol2:
        return False
    return Chem.MolToSmiles(mol1, True) != Chem.MolToSmiles(mol2, True)


def execute(request: DataFunctionRequest) -> DataFunctionResponse:
    column_id = string_input_field(request, 'structureColumn')
    input_column = request.inputColumns[column_id]
    input_molecules = column_to_molecules(input_column)
    deprotected_molecules = [None if m is None else Deprotect(m) for m in input_molecules]
    changed = [compare_molecules(mol1, mol2) for mol1, mol2 in zip(input_molecules, deprotected_molecules)]
    output_molecules_column = molecules_to_column(deprotected_molecules, f'Deprotected {input_column.name}', DataType.STRING)
    changed_column = ColumnData(name='Changed', dataType=DataType.BOOLEAN, values=changed)
    response = DataFunctionResponse(outputColumns=[output_molecules_column, changed_column])
    return response

