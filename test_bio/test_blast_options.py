"""
Copyright (C) 2017 Anodyne Informatics, LLC
"""

import os
from typing import List
from unittest import TestCase, main

from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord

from ruse.bio.blast_search import BlastCreateAndSearch


class BlastOptions(TestCase):
    # meaningful options:
    #
    # evalue	real	10.0	Expect value (E) for saving hits
    # num_descriptions	integer	500	Show one-line descriptions for this number of database sequences.
    # num_alignments	integer	250	Show alignments for this number of database sequences.
    # max_target_seqs	Integer	500	Number of aligned sequences to keep.
    #
    # max_target_seqs cannot be used with num_alignments or num_descriptions in XML output format

    def test_options(self) -> None:
        sequences = self.load_sequences()
        blast = BlastCreateAndSearch()
        options = {'evalue': 20000.0, 'num_descriptions': 5000, 'num_alignments': 5000}
        blast.search_blast_sequences(sequences[0], sequences, options=options);
        with open(blast.search.output_file(), 'r') as fh:
            blast_record = NCBIXML.read(fh)
        self.assertTrue(float(blast_record.expect), 20000.0)
        blast.clean_search()

    def test_options2(self) -> None:
        sequences = self.load_sequences()
        blast = BlastCreateAndSearch()
        options = {'evalue': 200000.0, 'max_target_seqs': 5000}
        blast.search_blast_sequences(sequences[0], sequences, options=options)
        with open(blast.search.output_file(), 'r') as fh:
            blast_record = NCBIXML.read(fh)
        self.assertTrue(float(blast_record.expect), 200000.0)
        blast.clean_search()

    def test_bad_options(self) -> None:
        sequences = self.load_sequences()
        blast = BlastCreateAndSearch()
        options = {'evalue': 20000.0, 'num_descriptions': 5000, 'max_target_seqs': 5000}
        self.assertRaises(ValueError, blast.search_blast_sequences, sequences[0], sequences, options=options)
        blast.clean_search()

    @classmethod
    def load_sequences(cls) -> List[SeqRecord]:
        file = 'ls_orchid.fasta'
        sequence_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources", file))
        return list(SeqIO.parse(sequence_file, 'fasta'))


if __name__ == "__main__":
    main()
