
# The contents of the data function definition script appear below
# the code should define an execute method

# Once this is working copy the code below to the script tag in the data function definition

from df.bio_helper import column_to_sequences, sequences_to_column
from df.data_transfer import DataFunctionRequest, DataFunctionResponse, string_input_field

from Bio.Data import CodonTable

def execute(request: DataFunctionRequest) -> DataFunctionResponse:
    column_id = string_input_field(request, 'sequenceColumn')
    input_column = request.inputColumns[column_id]
    input_sequences = column_to_sequences(input_column)
    codon_table_name = string_input_field(request, 'codonTableName', 'Standard')
    init_site_method = string_input_field(request, 'initMethod')

    # reduce sequences based on initiation site
    init_sequences = input_sequences.copy()
    if init_site_method == 'ATG':
        for n, s in enumerate(input_sequences):
            idx = s.seq.upper().find('ATG')
            if idx < 0:
                init_sequences[n] = None
            else:
                init_sequences[n] = s[idx:]
    elif init_site_method == 'table':
        codon_table = CodonTable.unambiguous_dna_by_name[codon_table_name]
        init_codons = codon_table.start_codons
        for n, s in enumerate(input_sequences):
            idx = [v for v in [s.seq.upper().find(codon) for codon in init_codons] if v != -1]
            if len(idx) == 0:
                init_sequences[n] = None
            else:
                init_sequences[n] = s[min(idx):]

    output_sequences = [None if s is None else s.translate(codon_table_name) for s in init_sequences]
    output_column = sequences_to_column(output_sequences, f'Translated {input_column.name}', genbank_output=False)
    response = DataFunctionResponse(
        outputColumns=[output_column])
    return response

