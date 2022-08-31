"""
Copyright (C) 2017 Anodyne Informatics, LLC
"""

import os
from Bio import SeqIO


def get_protein_query():
    input_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources", "phylodrome.fasta"))
    with open(input_file, 'r') as fh:
        return next(SeqIO.parse(fh, 'fasta'))


def get_protein_targets():
    input_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources", "phylodrome.fasta"))
    with open(input_file, 'r') as fh:
        return list(SeqIO.parse(fh, 'fasta'))


def get_dna_query():
    input_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources", "ls_orchid.fasta"))
    with open(input_file, 'r') as fh:
        return next(SeqIO.parse(fh, 'fasta'))


def get_dna_targets():
    input_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources", "ls_orchid.fasta"))
    with open(input_file, 'r') as fh:
        return list(SeqIO.parse(fh, 'fasta'))
