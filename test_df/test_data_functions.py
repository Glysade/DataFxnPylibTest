import glob
import importlib
import os
import shutil
from typing import Callable, Tuple
from unittest import TestCase, main, skip

from df.MmpdbColumnSearch import generate_mmpdb_dir
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


class DataFunctionTest(TestCase):

    def test_named_data_function_antibody_numbering(self) -> None:
        file_in = os.path.join(os.path.dirname(__file__), 'resources', 'antibody_numbering.json')
        _, response = run_named_data_function(file_in)
        self.assertTrue(response)
        self.assertEqual(1, len(response.outputColumns))
        self.assertEqual(0, len(response.outputTables))
        self.assertEqual(18, len(response.outputColumns[0].values))

    def test_named_data_function_blast_local_column_search(self) -> None:
        file_in = os.path.join(os.path.dirname(__file__), 'resources', 'blast_local_column_search.json')
        _, response = run_named_data_function(file_in)
        self.assertTrue(response)
        self.assertEqual(0, len(response.outputColumns))
        self.assertEqual(1, len(response.outputTables))
        self.assertEqual(10, len(response.outputTables[0].columns))
        self.assertEqual(1800, len(response.outputTables[0].columns[0].values))

    def test_named_data_function_blast_local_text_search(self) -> None:
        file_in = os.path.join(os.path.dirname(__file__), 'resources', 'blast_local_text_search.json')
        _, response = run_named_data_function(file_in)
        self.assertTrue(response)
        self.assertEqual(0, len(response.outputColumns))
        self.assertEqual(1, len(response.outputTables))
        self.assertEqual(9, len(response.outputTables[0].columns))
        self.assertEqual(101, len(response.outputTables[0].columns[0].values))

    @skip('Blast web search takes too long')
    def test_named_data_function_blast_web_search(self) -> None:
        file_in = os.path.join(os.path.dirname(__file__), 'resources', 'blast_web_search.json')
        _, response = run_named_data_function(file_in)
        self.assertTrue(response)
        self.assertEqual(0, len(response.outputColumns))
        self.assertEqual(1, len(response.outputTables))
        self.assertEqual(10, len(response.outputTables[0].columns))
        self.assertEqual(101, len(response.outputTables[0].columns[0].values))

    def test_named_data_function_blast_tablr_search(self) -> None:
        file_in = os.path.join(os.path.dirname(__file__), 'resources', 'blast_table_search.json')
        _, response = run_named_data_function(file_in)
        self.assertTrue(response)
        self.assertEqual(4, len(response.outputColumns))
        self.assertEqual(0, len(response.outputTables))
        self.assertEqual(18, len(response.outputColumns[0].values))

    def test_script_deprotect(self) -> None:
        from test_df.deprotect import execute
        file_in = os.path.join(os.path.dirname(__file__), 'resources', 'deprotect.json')
        _, response = run_script(file_in, execute)
        self.assertTrue(response)
        self.assertEqual(2, len(response.outputColumns))
        self.assertEqual(0, len(response.outputTables))
        self.assertEqual(100, len(response.outputColumns[0].values))

    def test_named_data_function_exact_mass(self) -> None:
        file_in = os.path.join(os.path.dirname(__file__), 'resources', 'exact_mass_df.json')
        _, response = run_named_data_function(file_in)
        self.assertTrue(response)
        self.assertEqual(1, len(response.outputColumns))
        self.assertEqual(0, len(response.outputTables))
        self.assertEqual(100, len(response.outputColumns[0].values))

    def test_exact_mass_script(self) -> None:
        from test_df.exact_mass_script import execute
        file_in = os.path.join(os.path.dirname(__file__), 'resources', 'exact_mass_script.json')
        _, response = run_script(file_in, execute)
        self.assertTrue(response)

    def test_named_data_function_mmpdb_table_column_search(self) -> None:
        file_in = os.path.join(os.path.dirname(__file__), 'resources', 'mmpdb_table_column_search.json')
        request = request_from_file(file_in)
        # unfortunately we can't delete the mmpdb directory after the test as the current process will
        # be locking it- so remove any old test data now
        mmpdb_dir = generate_mmpdb_dir(request)
        if os.path.exists(mmpdb_dir):
            shutil.rmtree(mmpdb_dir)
        for _ in range(2):
            request, response = run_named_data_function(file_in)
            self.assertTrue(response)
            self.assertEqual(0, len(response.outputColumns))
            self.assertEqual(1, len(response.outputTables))
            self.assertEqual(10, len(response.outputTables[0].columns))
            self.assertEqual(81, len(response.outputTables[0].columns[0].values))

    def test_named_data_function_mmpdb_database_search(self) -> None:
        file_in = os.path.join(os.path.dirname(__file__), 'resources', 'mmpdb_database_search.json')
        _, response = run_named_data_function(file_in)
        self.assertTrue(response)

        self.assertEqual(0, len(response.outputColumns))
        self.assertEqual(1, len(response.outputTables))
        self.assertEqual(15, len(response.outputTables[0].columns))
        self.assertEqual(63, len(response.outputTables[0].columns[0].values))

    def test_named_data_function_igblastn_local_text_search(self) -> None:
        file_in = os.path.join(os.path.dirname(__file__), 'resources', 'igblastn_local_text_search.json')
        _, response = run_named_data_function(file_in)
        self.assertTrue(response)

        self.assertEqual(0, len(response.outputColumns))
        self.assertEqual(1, len(response.outputTables))
        self.assertEqual(6, len(response.outputTables[0].columns))
        self.assertEqual(10, len(response.outputTables[0].columns[0].values))

    def test_named_data_function_igblastp_local_text_search(self) -> None:
        file_in = os.path.join(os.path.dirname(__file__), 'resources', 'igblastp_local_text_search.json')
        _, response = run_named_data_function(file_in)
        self.assertTrue(response)

        self.assertEqual(0, len(response.outputColumns))
        self.assertEqual(1, len(response.outputTables))
        self.assertEqual(6, len(response.outputTables[0].columns))
        self.assertEqual(4, len(response.outputTables[0].columns[0].values))

    def test_named_data_function_multiple_sequence_alignment(self) -> None:
        file_in = os.path.join(os.path.dirname(__file__), 'resources', 'multiple_sequence_alignment.json')
        _, response = run_named_data_function(file_in)
        self.assertTrue(response)

        self.assertEqual(3, len(response.outputColumns))
        self.assertEqual(0, len(response.outputTables))
        self.assertEqual(18, len(response.outputColumns[0].values))

    def test_named_data_function_pairwise_sequence_alignment(self) -> None:
        file_in = os.path.join(os.path.dirname(__file__), 'resources', 'pairwise_sequence_alignment.json')
        _, response = run_named_data_function(file_in)
        self.assertTrue(response)

        self.assertEqual(4, len(response.outputColumns))
        self.assertEqual(0, len(response.outputTables))
        self.assertEqual(18, len(response.outputColumns[0].values))

    def test_named_data_function_identify_enzyme_restriction_sites(self) -> None:
        file_in = os.path.join(os.path.dirname(__file__), 'resources', 'identify_enzyme_restriction_sites.json')
        _, response = run_named_data_function(file_in)
        self.assertTrue(response)

        self.assertEqual(1, len(response.outputColumns))
        self.assertEqual(0, len(response.outputTables))
        self.assertEqual(21, len(response.outputColumns[0].values))

    def test_named_data_function_translate_sequences(self) -> None:
        file_in = os.path.join(os.path.dirname(__file__), 'resources', 'translate_sequences.json')
        _, response = run_named_data_function(file_in)
        self.assertTrue(response)

        self.assertEqual(1, len(response.outputColumns))
        self.assertEqual(0, len(response.outputTables))
        self.assertEqual(21, len(response.outputColumns[0].values))

    def test_translate_sequences_script(self) -> None:
        from test_df.translate_sequences_script import execute
        file_in = os.path.join(os.path.dirname(__file__), 'resources', 'translate_sequences_script.json')
        _, response = run_script(file_in, execute)
        self.assertTrue(response)
        self.assertEqual(1, len(response.outputColumns))
        self.assertEqual(0, len(response.outputTables))
        self.assertEqual(21, len(response.outputColumns[0].values))

    def test_named_data_function_translate_sequences_all_reading_frames(self) -> None:
        file_in = os.path.join(os.path.dirname(__file__), 'resources', 'translate_sequences_all_reading_frames.json')
        _, response = run_named_data_function(file_in)
        self.assertTrue(response)

        self.assertEqual(0, len(response.outputColumns))
        self.assertEqual(1, len(response.outputTables))
        self.assertEqual(8, len(response.outputTables[0].columns))
        self.assertEqual(68, len(response.outputTables[0].columns[0].values))

    def test_named_data_function_reaction_table_search(self) -> None:
        file_in = os.path.join(os.path.dirname(__file__), 'resources', 'reaction_table_search.json')
        _, response = run_named_data_function(file_in)
        self.assertTrue(response)

        self.assertEqual(0, len(response.outputColumns))
        self.assertEqual(1, len(response.outputTables))
        self.assertEqual(7, len(response.outputTables[0].columns))
        self.assertEqual(2, len(response.outputTables[0].columns[0].values))

    def test_extract_genbank_regions(self) -> None:
        file_in = os.path.join(os.path.dirname(__file__), 'resources', 'extract_genbank_regions.json')
        _, response = run_named_data_function(file_in)
        self.assertTrue(response)

        self.assertEqual(14, len(response.outputColumns))
        self.assertEqual(0, len(response.outputTables))
        self.assertEqual(18, len(response.outputColumns[0].values))

    def test_named_near_neighbors_category_transition_search(self) -> None:
        file_in = os.path.join(os.path.dirname(__file__), 'resources',
                               'near_neighbors_free_wilson_transition_search.json')
        _, response = run_named_data_function(file_in)
        self.assertTrue(response)
        self.assertEqual(1, len(response.outputColumns))
        self.assertEqual(0, len(response.outputTables))
        self.assertEqual(10628, len(response.outputColumns[0].values))

    def test_named_antibody_numbering_with_bad_sequence(self) -> None:
        file_in = os.path.join(os.path.dirname(__file__), 'resources',
                               'antibody_numbering_bad_sequence.json')
        _, response = run_named_data_function(file_in)
        self.assertTrue(response)
        self.assertEqual(1, len(response.outputColumns))
        self.assertEqual(0, len(response.outputTables))
        values = response.outputColumns[0].values
        self.assertEqual(6, len(values))
        self.assertTrue(values[2] is None)

    def test_named_ab_component_analysis(self) -> None:
        file_in = os.path.join(os.path.dirname(__file__), 'resources',
                               'ab_component_analysis.json')
        _, response = run_named_data_function(file_in)
        self.assertTrue(response)
        self.assertEqual(0, len(response.outputColumns))
        self.assertEqual(1, len(response.outputTables))
        self.assertEqual(9, len(response.outputTables[0].columns))
        for col in response.outputTables[0].columns:
            self.assertEqual(220, len(col.values))

    def test_duplicate_msa_pair(self) -> None:
        file_in = os.path.join(os.path.dirname(__file__), 'resources', 'msa_error_duplicate_pair.json')
        _, response = run_named_data_function(file_in)
        self.assertTrue(response)
        self.assertEqual(3, len(response.outputColumns))
        self.assertEqual(0, len(response.outputTables))
        self.assertEqual(2, len(response.outputColumns[0].values))

    @classmethod
    def tearDownClass(cls) -> None:
        clean_output_files()


if __name__ == '__main__':
    main()
