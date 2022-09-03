"""
Copyright (C) 2017 Anodyne Informatics, LLC
"""

from typing import List, Tuple
from unittest import TestCase, main

from Bio.SeqRecord import SeqRecord

from ruse.bio.bio_util import is_dna, is_protein
from ruse.bio.blast_parse import BlastResults
from ruse.bio.blast_search import BlastCreateAndSearch, BlastSearchType
from test_pylib.test_bio.helper import get_dna_query, get_protein_query, get_dna_targets, get_protein_targets


class BlastDatabaseTestCase(TestCase):
    def test_dna_query_dan_targets(self) -> None:
        for search_type in None, BlastSearchType.BLASTN:
            blast, results = self.run_search(get_dna_query(), get_dna_targets(), search_type);
            self.assertTrue(blast.search_type, BlastSearchType.BLASTN)
            self.assertTrue(blast.search.search_type, BlastSearchType.BLASTN)
            self.assertTrue(blast.query_type, 'nucl')
            self.assertTrue(blast.target_type, 'nucl')
            self.assertTrue(is_dna(results.hits[0].query))
            self.assertTrue(is_dna(results.hits[0].target))

    def test_protein_query_protein_targets(self) -> None:
        for search_type in None, BlastSearchType.BLASTP:
            blast, results = self.run_search(get_protein_query(), get_protein_targets(), search_type)
            self.assertTrue(blast.search_type, BlastSearchType.BLASTP)
            self.assertTrue(blast.search.search_type, BlastSearchType.BLASTP)
            self.assertTrue(blast.query_type, 'prot')
            self.assertTrue(blast.target_type, 'prot')
            self.assertFalse(is_dna(results.hits[0].query))
            self.assertFalse(is_dna(results.hits[0].target))
            self.assertTrue(is_protein(results.hits[0].query))
            self.assertTrue(is_protein(results.hits[0].target))

    def test_dna_query_protein_targets(self) -> None:
        for search_type in None, BlastSearchType.BLASTX:
            blast, results = self.run_search(get_dna_query(), get_protein_targets(), search_type);
            self.assertTrue(blast.search_type, BlastSearchType.BLASTX)
            self.assertTrue(blast.search.search_type, BlastSearchType.BLASTX)
            self.assertTrue(blast.query_type, 'nucl')
            self.assertTrue(blast.target_type, 'prot')
            self.assertFalse(is_dna(results.hits[0].query))
            self.assertFalse(is_dna(results.hits[0].target))
            self.assertTrue(is_protein(results.hits[0].query))
            self.assertTrue(is_protein(results.hits[0].target))

    def test_protein_query_dna_targets(self) -> None:
        for search_type in None, BlastSearchType.TBLASTN:
            blast, results = self.run_search(get_protein_query(), get_dna_targets(), search_type);
            self.assertTrue(blast.search_type, BlastSearchType.TBLASTN)
            self.assertTrue(blast.search.search_type, BlastSearchType.TBLASTN)
            self.assertTrue(blast.query_type, 'prot')
            self.assertTrue(blast.target_type, 'nucl')
            self.assertFalse(is_dna(results.hits[0].query))
            self.assertFalse(is_dna(results.hits[0].target))
            self.assertTrue(is_protein(results.hits[0].query))
            self.assertTrue(is_protein(results.hits[0].target))

    def test_tblastx(self) -> None:
        blast, results = self.run_search(get_dna_query(), get_dna_targets(), BlastSearchType.TBLASTX);
        self.assertTrue(blast.search_type, BlastSearchType.TBLASTX)
        self.assertTrue(blast.search.search_type, BlastSearchType.TBLASTX)
        self.assertTrue(blast.query_type, 'nucl')
        self.assertTrue(blast.target_type, 'nucl')
        self.assertFalse(is_dna(results.hits[0].query))
        self.assertFalse(is_dna(results.hits[0].target))
        self.assertTrue(is_protein(results.hits[0].query))
        self.assertTrue(is_protein(results.hits[0].target))

    @classmethod
    def run_search(cls, query: SeqRecord, targets: List[SeqRecord], search_type: BlastSearchType = None) -> Tuple[
        BlastCreateAndSearch, BlastResults]:
        blast = BlastCreateAndSearch();
        blast.search_blast_sequences(query, targets, search_type, options={'evalue': 10000});
        results = BlastResults()
        results.parse(blast.search.output_file())
        blast.clean_search();
        return blast, results;


if __name__ == "__main__":
    main()
