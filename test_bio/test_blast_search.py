"""
Copyright (C) 2017 Anodyne Informatics, LLC
"""

import os
from unittest import TestCase, main, skip

from Bio import SeqIO

from ruse.bio.bio_util import is_defined_sequence
from ruse.bio.blast_parse import MultipleBlastResults
from ruse.bio.blast_search import BlastSearch, BlastWebSearch, BlastDatabase, BlastSearchType, BlastCreateAndSearch
from ruse.bio.blast_utils import SequenceMatcher


class BlastDatabaseTestCase(TestCase):

    def test_db_build(self) -> None:
        blast_database = self.build_database("ls_orchid.fasta", "orchid")
        blast_database.clean_database()

    def test_db_build2(self) -> None:
        blast_database = self.build_database("phylodrome.fasta", "phylodrome")
        blast_database.clean_database()

    def build_database(self, input_file: str, name: str) -> BlastDatabase:
        input_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources", input_file))
        sequences = list(SeqIO.parse(input_file, 'fasta'))
        blast_database = BlastDatabase(name)
        blast_database.build_database(sequences)
        return blast_database

    def test_search(self) -> None:
        blast_database = self.build_database("ls_orchid.fasta", "orchid")
        search = self.search_database("ls_orchid.fasta", "orchid", "blastn")
        blast_database.clean_database()
        search.clean_search()

    @classmethod
    def search_database(cls, sequence_file: str, database_name: str, search_type: str) -> BlastSearch:
        sequence_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources", sequence_file))
        with open(sequence_file, 'r') as fh:
            query = next(SeqIO.parse(fh, 'fasta'))
        search = BlastSearch()
        search.search_blast_database(query, database_name, BlastSearchType.from_string(search_type))
        return search


class BlastCreateAndSearchTestCase(TestCase):

    def test_search(self) -> None:
        self.search_database("ls_orchid.fasta")

    @classmethod
    def search_database(cls, sequence_file: str) -> BlastCreateAndSearch:
        sequence_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources", sequence_file))
        sequences = list(SeqIO.parse(sequence_file, 'fasta'))
        blast = BlastCreateAndSearch()
        blast.search_blast_sequences(sequences[0], sequences)
        blast.clean_search()
        return blast


class BlastWebSearchTestCase(TestCase):

    @skip("Skipping web blast search due to possible timeouts")
    def test_search(self) -> None:
        self.search_database("ls_orchid.fasta", 'nucl')

    @skip("Skipping web blast search due to possible timeouts")
    def test_search_prot(self) -> None:
        self.search_database("phylodrome.fasta", 'prot')

    def search_database(self, sequence_file: str, query_type: str) -> str:
        sequence_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources", sequence_file))
        with open(sequence_file, 'r') as fh:
            query = next(SeqIO.parse(sequence_file, 'fasta'))

        search = BlastWebSearch()
        search.search_blast_database(query)

        if query_type == 'nucl':
            self.assertEqual(search.database_name, 'nt')
            self.assertEqual(search.search_type, BlastSearchType.BLASTN)
        else:
            self.assertEqual(search.database_name, 'nr')
            self.assertEqual(search.search_type, BlastSearchType.BLASTP)

        search.clean_search()
        return search.search_name


class MultipleBlastLocalSearchTestCase(TestCase):
    """Search multiple sequences records against a local blast database"""

    def test_dna_search(self) -> None:
        self.search_local_database('rat_nt', 'blastn', 'rat.fasta')

    def test_protein_search(self) -> None:
        results = self.search_local_database('human_nr', 'blastp', 'phylodrome.fasta')
        self.assertTrue('desert hedgehog' in results.query_hits[0].hits[0].target_def)

    def search_local_database(self, database: str, method: str, query_file: str):
        sequence_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources", query_file))
        sequences = list(SeqIO.parse(sequence_file, 'fasta'))
        search = BlastSearch()
        search.multiple_query_search_blast_database(sequences, database, BlastSearchType.from_string(method))
        self.assertTrue(search.error is None)
        results = MultipleBlastResults()
        results.parse(search.output_file())
        search.clean_search()
        self.assertEqual(len(sequences), len(results.query_hits))
        self.assertEqual(results.database, database)
        for query, query_result in zip(sequences, results.query_hits):
            self.assertEqual(query_result.database, database)
            self.assertEqual(query_result.query_def, query.description)
            for hit in query_result.hits:
                result_query_seq = hit.query.replace('-', '')
                self.assertTrue(result_query_seq in str(query.seq))

        for i in range(2):
            if i == 1:
                results.retrieve_targets(True)
            else:
                results.retrieve_local_targets()
            for query_result in results.query_hits:
                for hit in query_result.hits:
                    self.assertTrue(hit.target_record is not None);
            results.complement_target_sequence()
            for query_result in results.query_hits:
                for hit in query_result.hits:
                    result_target_seq = hit.target.replace('-', '')
                    if results.search_type == BlastSearchType.BLASTN:
                        self.assertTrue(SequenceMatcher.match_sequences(str(hit.target_record.seq), result_target_seq) != -1)
                    else:
                        self.assertTrue(result_target_seq in str(hit.target_record.seq))
                    hit.target_record = None

        return results


if __name__ == "__main__":
    main()
