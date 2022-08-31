"""
Copyright (C) 2017 Anodyne Informatics, LLC
"""

import os
import unittest
from typing import List, Dict
from unittest import TestCase, main, skip

from Bio import AlignIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from ruse.bio.phylo_tree import PhyloTree
from ruse.bio.sequence_align import MultipleSequenceAlignment, SequenceAlignmentMethod, _build_new_features


class MultipleSequenceAlignmentTestCase(TestCase):
    """A class to test MSA"""

    def test_simple_alignment(self) -> None:
        """Test alignment of some DNA sequences"""

        input_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources/ls_orchid.fasta"))
        sequences = list(SeqIO.parse(input_file, 'fasta'))
        self.assertEqual(len(sequences), 94)
        self.do_simple_alignment(sequences)

    def test_simple_protein_alignment(self) -> None:
        """Test alignment of some protein sequences"""

        input_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources/Human_kinase_protein.fasta"))
        sequences = list(SeqIO.parse(input_file, 'fasta'))
        self.assertEqual(len(sequences), 516)
        sequences = sequences[:25]
        self.do_simple_alignment(sequences)

    def test_simple_muscle_alignment(self) -> None:
        """Test alignment of some DNA sequences using Muscle instead of ClustalO"""

        input_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources/ls_orchid.fasta"))
        sequences = list(SeqIO.parse(input_file, 'fasta'))
        self.assertEqual(len(sequences), 94)
        self.do_simple_alignment(sequences, alignment_method=SequenceAlignmentMethod.MUSCLE)

    @skip("ClustalW not supported (and test not working)")
    def test_simple_clustalw_alignment(self) -> None:
        """Test alignment of some DNA sequences using ClustalW instead of ClustalO"""

        input_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources/ls_orchid.fasta"))
        sequences = list(SeqIO.parse(input_file, 'fasta'))
        self.assertEqual(len(sequences), 94)
        self.do_simple_alignment(sequences, alignment_method=SequenceAlignmentMethod.CLUSTALW)

    @skip("Covid RNA alignment takes ~ 2 hours")
    def test_covid_alignment(self) -> None:
        """Test alignment of 23 coronavirus sequences"""

        input_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources/SARS-MERS-SARS2.fasta"))
        sequences = list(SeqIO.parse(input_file, 'fasta'))
        self.assertEqual(len(sequences), 23)
        self.do_simple_alignment(sequences, alignment_method=SequenceAlignmentMethod.MUSCLE)

    def do_simple_alignment(self, sequences: List[SeqRecord],
                            alignment_method: SequenceAlignmentMethod = SequenceAlignmentMethod.CLUSTALO) -> None:
        options = {}
        msa = MultipleSequenceAlignment()
        msa.align_sequences(options, sequences, alignment_method=alignment_method)

        tree = PhyloTree(msa.tree)
        self.assertEqual(tree.data_tree, tree.tree_to_data_recursive())
        self.assertEqual(set([seq.id for seq in sequences]), set(self.extract_ids_from_tree(tree.data_tree)))

        output_sequences = msa.aligned_sequences
        self.assertEqual(len(sequences), len(output_sequences))
        for s1, s2 in zip(sequences, output_sequences):
            self.assertEqual(s1.id, s2.id)
            self.assertEqual(s1.seq, s2.seq.ungap('-'))

        msa.add_distances()
        msa.clean_files()
        for path in ["{}_in.fasta", "{}_out.aln"]:
            self.assertFalse(os.path.exists(path.format(msa.basename)))

    def test_genbank_align(self) -> None:
        """Test alignment using genbank sequences feature mapping"""
        input_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources/example_nt.gb"))
        sequences = list(SeqIO.parse(input_file, 'gb'))[:5]

        msa = MultipleSequenceAlignment()
        msa.align_sequences({}, sequences)
        phyloTree = PhyloTree(msa.tree)
        self.assertEqual(phyloTree.data_tree, phyloTree.tree_to_data_recursive())
        self.assertEqual(set([seq.id for seq in sequences]), set(self.extract_ids_from_tree(phyloTree.data_tree)))

        msa.copy_features_to_aligned_sequences()
        output_sequences = msa.aligned_sequences
        for s1, s2 in zip(sequences, output_sequences):
            for a, b in zip(s1.features, s2.features):
                self.assertEqual(a.type, b.type)

        msa.clean_files()

    @classmethod
    def extract_ids_from_tree(cls, tree: Dict):
        ids = []
        cls.extract_ids_from_tree_data(tree, ids)
        return ids

    @classmethod
    def extract_ids_from_tree_data(cls, tree: Dict, ids: List[str]):
        if 'name' in tree:
            ids.append(tree['name'])
        if 'children' in tree:
            for child in tree['children']:
                cls.extract_ids_from_tree_data(child, ids)

    def test_feature_mapping(self) -> None:
        """Validate that feature mapping of genbank files is working correctly"""

        input_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources/example_nt.gb"))
        sequences = list(SeqIO.parse(input_file, 'gb'))
        align_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources/example_nt.aln"))
        align = AlignIO.read(align_file, "clustal")
        aligned_sequences = []
        for seq in align:
            aligned_sequences.append(seq)

        msa = MultipleSequenceAlignment()
        msa.input_sequences = sequences
        msa.aligned_sequences = aligned_sequences
        msa.alignment = align

        msa.copy_features_to_aligned_sequences()

        for input_seq, aligned_seq in zip(msa.input_sequences, msa.aligned_sequences):
            self.assertEqual(input_seq.seq, aligned_seq.seq.ungap('-'))
            self.assertEqual(len(input_seq.features), len(aligned_seq.features))
            for input_feature, aligned_feature in zip(input_seq.features, aligned_seq.features):
                in_feature_seq = input_feature.extract(input_seq.seq);
                out_feature_seq = aligned_feature.extract(aligned_seq.seq)
                self.assertEqual(in_feature_seq, out_feature_seq.ungap('-'))

    def test_feature_copying(self) -> None:
        """Validate that feature mapping of genbank files is working correctly using self-mapping"""

        input_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources/NT_079573.5.gb"))
        sequence = list(SeqIO.parse(input_file, 'gb'))[0]

        new_record = SeqRecord(sequence.seq, id='test')
        self.assertEqual(len(new_record.features), 0)
        positions = range(len(sequence) + 1);

        _build_new_features(sequence, new_record, positions)
        self.assertEqual(len(sequence.features), len(new_record.features))
        for feature, new_feature in zip(sequence.features, new_record.features):
            in_feature_seq = feature.extract(sequence.seq)
            out_feature_seq = new_feature.extract(new_record.seq)
            self.assertEqual(in_feature_seq, out_feature_seq)


if __name__ == "__main__":
    main()
