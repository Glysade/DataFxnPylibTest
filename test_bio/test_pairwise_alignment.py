import os
from typing import List
from unittest import TestCase, main

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from ruse.bio.pairwise_alignment import needle_pairwise_alignment, water_pairwise_alignment


class PairwiseSequenceAlignmentTestCase(TestCase):

    def test_needle_simple_protein_alignment(self) -> None:
        """Test alignment of two protein sequences"""

        input_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources/Human_kinase_protein.fasta"))
        sequences = list(SeqIO.parse(input_file, 'fasta'))
        self.assertEqual(len(sequences), 516)

        query = sequences[10]
        targets = sequences[15:20]
        self.run_needle_alignments(query, targets)

    def test_needle_genbank_align(self) -> None:
        """Test alignment using genbank sequences feature mapping"""
        input_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources/example_nt.gb"))
        sequences = list(SeqIO.parse(input_file, 'gb'))[:10]

        query = sequences[3]
        targets = sequences[4:8]

        self.run_needle_alignments(query, targets)

    def test_water_simple_protein_alignment(self) -> None:
        """Test alignment of two protein sequences"""

        input_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources/Human_kinase_protein.fasta"))
        sequences = list(SeqIO.parse(input_file, 'fasta'))
        self.assertEqual(len(sequences), 516)

        query = sequences[10]
        targets = sequences[15:20]
        self.run_water_alignments(query, targets)

    def test_water_genbank_align(self) -> None:
        """Test alignment using genbank sequences feature mapping"""
        input_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources/example_nt.gb"))
        sequences = list(SeqIO.parse(input_file, 'gb'))[:10]

        query = sequences[3]
        targets = sequences[4:8]

        self.run_water_alignments(query, targets)

    def run_needle_alignments(self, query: SeqRecord, targets: List[SeqRecord]) -> None:

        results = needle_pairwise_alignment(query, targets)
        self.assertEqual(len(results), len(targets))

        for result, target in zip(results, targets):
            self.assertEqual(query.id, result.query.id)
            self.assertEqual(query.seq, result.query.seq.ungap('-'))
            self.assertEqual(target.id, result.target.id)

            self.assertEqual(target.seq, result.target.seq.ungap('-'))

            if query.features:
                for a, b in zip(query.features, result.query.features):
                    self.assertEqual(a.type, b.type)
            if target.features:
                for a, b in zip(target.features, result.target.features):
                    self.assertEqual(a.type, b.type)

    def run_water_alignments(self, query: SeqRecord, targets: List[SeqRecord]) -> None:

        results = water_pairwise_alignment(query, targets)
        self.assertEqual(len(results), len(targets))

        for result, target in zip(results, targets):
            self.assertEqual(query.id, result.query.id)
            self.assertTrue(result.query.seq.ungap('-') in query.seq)
            self.assertEqual(target.id, result.target.id)

            self.assertTrue(result.target.seq.ungap('-') in target.seq)

            start = str(query.seq).index(str(result.query.seq.ungap('-')))
            end = start + len(result.query.seq)
            matched_query = query[start:end]

            start = str(target.seq).index(str(result.target.seq.ungap('-')))
            end = start + len(result.target.seq)
            matched_target = target[start:end]

            if matched_query.features:
                for a, b in zip(matched_query.features, result.query.features):
                    self.assertEqual(a.type, b.type)
            if target.features:
                for a, b in zip(matched_target.features, result.target.features):
                    self.assertEqual(a.type, b.type)


if __name__ == "__main__":
    main()
