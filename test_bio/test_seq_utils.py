"""
Copyright (C) 2017 Anodyne Informatics, LLC
"""
import os
from unittest import TestCase, main

from Bio import SeqIO

from ruse.bio.bio_util import extract_feature
from ruse.bio.blast_utils import retrieve_entrez_records


class TestSeqUtils(TestCase):
    """A class for testing blast related utilities"""

    def test_download_pdb_entry(self) -> None:
        records = retrieve_entrez_records('Protein', ['4KMH_A'])
        self.assertEqual(len(records), 1)

    def test_extract_features(self) -> None:
        file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources", "example_nt.gb"))
        sequences_in = list(SeqIO.parse(file, 'gb'))
        sequences_out = [extract_feature(s, type='genE', qualifier_name='geNe', qualifier_value='hiSI') for s in sequences_in]
        self.assertEqual(len(sequences_in), 21)
        self.assertEqual(len([s for s in sequences_out if s]), 14)

if __name__ == "__main__":
    main()
