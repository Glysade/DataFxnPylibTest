import os
import uuid
from typing import List
from unittest import TestCase, main, skip

from Bio import SeqIO

from df.data_transfer import ColumnData, DataFunctionRequest, InputField, DataType
from test_df.test_chem_client import run_data_function_module
from test_df.test_data_functions import clean_output_files


def blast_db() -> str:
    if os.name == 'nt':
        return 'c:\\db\\ncbi'
    elif os.name == 'posix':
        return '/mnt/c/db/ncbi'
    else:
        raise ValueError(f'Unknown os {os.name}')


class BioClientTest(TestCase):

    def _check_columns_present(self, number_values: int, number_columns: int, output_columns: List[ColumnData]) -> None:
        self.assertEqual(number_columns, len(output_columns))
        for column_number in range(number_columns):
            self.assertEqual(number_values, len(output_columns[column_number].values))

    def test_antibody_numbering(self) -> None:
        seq_file = os.path.join(os.path.dirname(__file__), '../test_bio/resources/antibody_sequences.fasta')
        antibody_sequences = list(SeqIO.parse(seq_file, 'fasta'))[10:15]
        values = [str(s.seq) for s in antibody_sequences]

        input_field = InputField(dataType=DataType.STRING, data='sequenceColumn', id=str(uuid.uuid4()))
        column = ColumnData(name='sequence', dataType=DataType.STRING,
                            contentType='chemical/x-sequence', values=values)
        request = DataFunctionRequest(serviceName='AntibodyNumbering', inputColumns={'sequenceColumn': column},
                                      inputFields={'sequenceColumn': input_field},
                                      id=str(uuid.uuid4()))
        response = run_data_function_module(request)
        self._check_columns_present(len(values), 1, response.outputColumns)

    def test_sequence_alignment(self) -> None:
        seq_file = os.path.join(os.path.dirname(__file__), '../test_bio/resources/phylodrome.fasta')
        sequences = list(SeqIO.parse(seq_file, 'fasta'))
        values = [str(s.seq) for s in sequences]
        column = ColumnData(name='sequence', dataType=DataType.STRING,
                            contentType='chemical/x-sequence', values=values)
        request = DataFunctionRequest(serviceName='MultipleSequenceAlignment', inputColumns={'sequenceColumn': column},
                                      id=str(uuid.uuid4()))
        response = run_data_function_module(request)
        self._check_columns_present(len(values), 3, response.outputColumns)

    def test_blast_table_search(self) -> None:
        seq_file = os.path.join(os.path.dirname(__file__), '../test_bio/resources/antibody_sequences.fasta')
        sequences = list(SeqIO.parse(seq_file, 'fasta'))
        values = [str(s.seq) for s in sequences]
        column = ColumnData(name='sequence', dataType=DataType.STRING,
                            contentType='chemical/x-sequence', values=values)
        query_seq = values[0]
        query_input_field = InputField(id='query', contentType='chemical/x-sequence', data=query_seq,
                                       dataType=DataType.STRING)
        column_input_field = InputField(id='sequenceColumn', data='sequenceColumn',
                                        dataType=DataType.STRING)
        request = DataFunctionRequest(serviceName='BlastTableSearch', inputColumns={'sequenceColumn': column},
                                      inputFields={'query': query_input_field, 'sequenceColumn': column_input_field},
                                      id=str(uuid.uuid4()))
        response = run_data_function_module(request)
        self._check_columns_present(len(values), 4, response.outputColumns)

    def test_pairwise_alignment(self) -> None:
        seq_file = os.path.join(os.path.dirname(__file__), '../test_bio/resources/antibody_sequences.fasta')
        sequences = list(SeqIO.parse(seq_file, 'fasta'))[10:15]
        values = [str(s.seq) for s in sequences]

        column = ColumnData(name='sequence', dataType=DataType.STRING,
                            contentType='chemical/x-sequence', values=values)

        query_seq = values[0]
        query_input_field = InputField(id='query', contentType='chemical/x-sequence', data=query_seq,
                                       dataType=DataType.STRING)
        request = DataFunctionRequest(serviceName='PairwiseSequenceAlignment', inputColumns={'sequenceColumn': column},
                                      inputFields={'query': query_input_field}, id=str(uuid.uuid4()))
        response = run_data_function_module(request)
        self._check_columns_present(len(values), 4, response.outputColumns)

    def test_local_blast_text_sequence_search(self) -> None:
        seq_file = os.path.join(os.path.dirname(__file__), '..', 'test_bio', 'resources', 'phylodrome.fasta')
        with open(seq_file) as fh:
            sequence = next(SeqIO.parse(fh, 'fasta'))
        query = str(sequence.seq)[100:150]
        query_input_field = InputField(id='query', contentType='chemical/x-sequence', data=str(query),
                                       dataType=DataType.STRING)
        max_hits_input_field = InputField(id='maxHits', data=20, dataType=DataType.INTEGER)
        method_input_field = InputField(id='method', data='BLASTP', dataType=DataType.STRING)
        database_input_field = InputField(id='databaseName', data='human_nr', dataType=DataType.STRING)
        blast_db_field = InputField(id='blastDbPath', data=blast_db(),
                                    dataType=DataType.STRING)
        request = DataFunctionRequest(serviceName='BlastLocalTextSearch',
                                      inputFields={'query': query_input_field, 'maxHits': max_hits_input_field,
                                                   'method': method_input_field, 'databaseName': database_input_field,
                                                   'blastDbPath': blast_db_field},
                                      id=str(uuid.uuid4()))
        response = run_data_function_module(request)
        self.assertEqual(0, len(response.outputColumns))
        self.assertEqual(1, len(response.outputTables))
        self._check_columns_present(20, 8, response.outputTables[0].columns)

    def test_local_blast_table_sequence_search(self) -> None:
        seq_file = os.path.join(os.path.dirname(__file__), '..', 'test_bio', 'resources', 'phylodrome.fasta')
        sequences = list(SeqIO.parse(seq_file, 'fasta'))[0:3]
        values = [str(s.seq)[100:150] for s in sequences]

        column = ColumnData(name='sequence', dataType=DataType.STRING,
                            contentType='chemical/x-sequence', values=values)
        max_hits_input_field = InputField(id='maxHits', data=10, dataType=DataType.INTEGER)
        method_input_field = InputField(id='method', data='BLASTP', dataType=DataType.STRING)
        database_input_field = InputField(id='databaseName', data='human_nr', dataType=DataType.STRING)
        column_input_field = InputField(id='sequenceColumn', data='sequenceColumn',
                                        dataType=DataType.STRING)
        blast_db_field = InputField(id='blastDbPath', data=blast_db(),
                                    dataType=DataType.STRING)
        request = DataFunctionRequest(serviceName='BlastLocalColumnSearch', inputColumns={'sequenceColumn': column},
                                      inputFields={'maxHits': max_hits_input_field,
                                                   'method': method_input_field, 'databaseName': database_input_field,
                                                   'blastDbPath': blast_db_field, 'sequenceColumn': column_input_field},
                                      id=str(uuid.uuid4()))
        response = run_data_function_module(request)
        self.assertEqual(0, len(response.outputColumns))
        self.assertEqual(1, len(response.outputTables))
        self._check_columns_present(30, 10, response.outputTables[0].columns)

    @skip("Web blast search can take too long")
    def test_web_blast_search(self) -> None:
        seq_file = os.path.join(os.path.dirname(__file__), '..', 'test_bio', 'resources', 'phylodrome.fasta')
        with open(seq_file) as fh:
            sequence = next(SeqIO.parse(fh, 'fasta'))
        query = str(sequence.seq)[100:150]
        query_input_field = InputField(id='query', contentType='chemical/x-sequence', data=str(query),
                                       dataType=DataType.STRING)
        max_hits_input_field = InputField(id='maxHits', data=20, dataType=DataType.INTEGER)
        database_input_field = InputField(id='databaseName', data='nr', dataType=DataType.STRING)
        request = DataFunctionRequest(serviceName='BlastWebSearch',
                                      inputFields={'query': query_input_field, 'maxHits': max_hits_input_field,
                                                   'databaseName': database_input_field},
                                      id=str(uuid.uuid4()))
        response = run_data_function_module(request)
        self.assertEqual(0, len(response.outputColumns))
        self.assertEqual(1, len(response.outputTables))
        self._check_columns_present(20, 8, response.outputTables[0].columns)

    def test_translate_sequence(self) -> None:
        seq_file = os.path.join(os.path.dirname(__file__), '..', 'test_bio', 'resources', 'example_nt.gb')
        sequences = list(SeqIO.parse(seq_file, 'gb'))
        values = [str(s.seq) for s in sequences]

        column = ColumnData(name='sequence', dataType=DataType.STRING,
                            contentType='chemical/x-sequence', values=values)
        codon_table_field = InputField(id='codonTableName', dataType=DataType.STRING, data='Standard')
        column_input_field = InputField(id='sequenceColumn', data='sequenceColumn',
                                        dataType=DataType.STRING)
        request = DataFunctionRequest(serviceName='TranslateSequences',
                                      inputFields={'codonTableName': codon_table_field,
                                                   'sequenceColumn': column_input_field},
                                      inputColumns={'sequenceColumn': column},
                                      id=str(uuid.uuid4()))
        response = run_data_function_module(request)
        self.assertEqual(1, len(response.outputColumns))
        self._check_columns_present(21, 1, response.outputColumns)

    def test_enzyme_restriction(self) -> None:
        seq_file = os.path.join(os.path.dirname(__file__), '..', 'test_bio', 'resources', 'example_nt.gb')
        sequences = list(SeqIO.parse(seq_file, 'gb'))
        values = [str(s.seq) for s in sequences]

        column = ColumnData(name='sequence', dataType=DataType.STRING,
                            contentType='chemical/x-sequence', values=values)
        enzyme_field = InputField(id='enzymes', data=['EcoRI', 'PsiI'], dataType=DataType.STRING_LIST)
        column_input_field = InputField(id='sequenceColumn', data='sequenceColumn',
                                        dataType=DataType.STRING)
        request = DataFunctionRequest(serviceName='EnzymeRestriction',
                                      inputFields={'enzymes': enzyme_field, 'sequenceColumn': column_input_field},
                                      inputColumns={'sequenceColumn': column},
                                      id=str(uuid.uuid4()))
        response = run_data_function_module(request)
        self.assertEqual(1, len(response.outputColumns))
        self._check_columns_present(21, 1, response.outputColumns)

    def test_translate_orf(self) -> None:
        seq_file = os.path.join(os.path.dirname(__file__), '..', 'test_bio', 'resources', 'example_nt.gb')
        sequences = list(SeqIO.parse(seq_file, 'gb'))
        values = [str(s.seq) for s in sequences]

        column = ColumnData(name='sequence', dataType=DataType.STRING,
                            contentType='chemical/x-sequence', values=values)
        protein_length_input_field = InputField(id='minimumProteinLength', dataType=DataType.INTEGER, data=50)
        codon_table_field = InputField(id='codonTableName', dataType=DataType.STRING, data='Standard')
        column_input_field = InputField(id='sequenceColumn', data='sequenceColumn',
                                        dataType=DataType.STRING)

        request = DataFunctionRequest(serviceName='TranslateOpenReadingFrames',
                                      inputFields={'minimumProteinLength': protein_length_input_field,
                                                   'codonTableName': codon_table_field,
                                                   'sequenceColumn': column_input_field},
                                      inputColumns={'sequenceColumn': column},
                                      id=str(uuid.uuid4()))
        response = run_data_function_module(request)
        self.assertEqual(0, len(response.outputColumns))
        self.assertEqual(1, len(response.outputTables))
        self._check_columns_present(274, 8, response.outputTables[0].columns)

    @classmethod
    def tearDownClass(cls) -> None:
        clean_output_files()


if __name__ == "__main__":
    main()
