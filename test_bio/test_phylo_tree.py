"""
Copyright (C) 2017 Anodyne Informatics, LLC
"""

import os
from typing import List, Dict
from unittest import TestCase, main, skip

from Bio import AlignIO, Phylo
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from ruse.bio.phylo_tree import PhyloTree, PhyloTreeBuilder


class PhyloTreeTestCase(TestCase):

    @skip("skipping raxml phylo tree construction- too time consuming")
    def test_raxml_phylo_tree(self):
        phylo_tree_builder = PhyloTreeBuilder()
        align_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources/kinase_alignment.aln"))
        phylo_tree_builder.build_raxml_tree(align_file)
        phylo_tree_builder.cleanup()

    def test_fasttree_phylo_tree(self):
        phylo_tree_builder = PhyloTreeBuilder()
        align_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources/kinase_alignment.aln"))
        t = phylo_tree_builder.build_fasttree_tree(align_file)
        # Phylo.draw_ascii(t)
        phylo_tree_builder.cleanup()
        ids = [s.id for s in phylo_tree_builder.alignment]

        tree = PhyloTree(t)
        distances = tree.leaf_distances()
        for d in distances.values():
            self.assertTrue(d >= 0)
        for label in distances:
            self.assertTrue(label in ids)
        for id in ids:
            self.assertTrue(id in distances)


if __name__ == "__main__":
    main()
