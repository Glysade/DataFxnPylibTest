
# The contents of the data function definition script appear below
# the code should define an execute method

# Once this is working copy the code below to the script tag in the data function definition

from df.bio_helper import column_to_sequences, sequences_to_column
from df.data_transfer import DataFunction, DataFunctionRequest, DataFunctionResponse, string_input_field


def execute(request: DataFunctionRequest) -> DataFunctionResponse:
    column_id = string_input_field(request, 'sequenceColumn')
    input_column = request.inputColumns[column_id]
    input_sequences = column_to_sequences(input_column)
    codon_table_name = string_input_field(request, 'codonTableName', 'Standard')
    output_sequences = [None if s is None else s.translate(codon_table_name) for s in input_sequences]
    output_column = sequences_to_column(output_sequences, f'Translated {input_column.name}', genbank_output=False)
    response = DataFunctionResponse(
        outputColumns=[output_column])
    return response

