"""
Copyright (C) 2017 Anodyne Informatics, LLC
"""

import os
from unittest import TestCase, main

from Bio import SeqIO

from ruse.bio import blast_parse
from ruse.bio.blast_parse import build_common_alignments


class BlastDatabaseTestCase(TestCase):
    """Test parse of blast search results"""

    def test_parse_orchid(self) -> None:
        """Test parse of a local search"""
        file = 'ls_orchid_blast_nt_search.xml'
        file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources", file))
        results = blast_parse.BlastResults()
        results.parse(file)

        sequence_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources", "ls_orchid.fasta"))
        with open(sequence_file, 'r') as fh:
            query = next(SeqIO.parse(fh, 'fasta'))
        query = str(query.seq)
        self.assertEqual(query, results.hits[0].query)
        self.assertEqual(query, results.hits[0].target)
        self.assertEqual(0, results.hits[0].evalue)

        data = results.to_data()
        self.assertEqual(query, data['hits'][0]['query'])
        self.assertEqual(query, data['hits'][0]['target'])
        self.assertEqual(0, data['hits'][0]['evalue'])

    def test_retrieve_protein_targets(self) -> None:
        """Test parse of web search results and retrieval of genbank records"""
        file = "example-blast-web-protein-results.xml"
        file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources", file))
        self.retrieve_targets(file)

    def test_retrieve_dna_targets(self) -> None:
        """Test parse of web search results and retrieval of genbank records"""
        file = "example-blast-web-dna-results.xml"
        file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources", file))
        self.retrieve_targets(file)

    def retrieve_targets(self, file) -> None:
        """Test that we can download full length sequences for the blast matches and that the sequences match"""
        results = blast_parse.BlastResults()
        results.parse(file)
        results.retrieve_targets()

        for hit in results.hits:
            self.assertTrue(hit.target is not None)
            if hit.target_record:
                self.assertTrue(hit.target.replace('-', '') in str(hit.target_record.seq))

    def test_parse_multiple_queries(self) -> None:
        xml_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources", "phylodrome_multi_search.xml"))
        results = blast_parse.MultipleBlastResults()
        results.parse(xml_file)
        self.assertEqual(len(results.query_hits), 5)
        self.assertEqual(results.database, 'human_nr')
        self.assertEqual(results.query_hits[3].query_def, 'CAEEL_hedgehog-like_GRounDhog_grd-11_NP_507923.1')
        self.assertEqual(results.query_hits[3].query_id, 'Query_4')
        result = results.query_hits[0]
        self.assertEqual(len(result.hits), 30)
        self.assertTrue(result.hits[0].target_def.startswith('desert hedgehog'))
        self.assertEqual(result.hits[0].accession, 'NP_066382')

    def test_parse_to_alignment(self) -> None:
        seq_file = os.path.join(os.path.dirname(__file__), '..', 'test_bio', 'resources', 'phylodrome.fasta')
        with open(seq_file) as fh:
            sequence = next(x for i, x in enumerate(SeqIO.parse(fh, 'fasta')) if i == 1)
        xml_file = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "resources", "phylodrome_sequence2_web_blast.xml"))
        results = blast_parse.BlastResults()
        results.parse(xml_file)

        alignments = build_common_alignments(sequence, results.hits)
        query_align = alignments[0]
        for align in alignments:
            self.assertEqual(len(query_align), len(align))
        self.assertEqual(len(results.hits)+1, len(alignments))


if __name__ == "__main__":
    main()
