import os
import uuid
from typing import List
from unittest import TestCase, main

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from ruse.bio.applications import IgblastpCommandLine
from ruse.bio.blast_parse import build_common_alignments
from ruse.bio.igblast_parse import parse_igblastp_results, IgResult, parse_igblastn_results


class TestIgBlast(TestCase):

    def test_igblastp_run(self):

        id = str(uuid.uuid4())
        v_database = os.path.join(os.environ['IGDATA'], 'database', 'imgt.Homo_sapiens.AA.V.f.orf')
        organism = 'human'
        domain_system = 'imgt'
        format = '"7 std qseq sseq btop"'

        query_file = os.path.join(os.path.dirname(__file__), '..', 'test_bio', 'resources', 'igblastp.fasta')
        out_file = f'igblastp_{id}.out'
        try:
            command = IgblastpCommandLine(germline_db_V=v_database, organism=organism, domain_system=domain_system,
                                          outfmt=format, query=query_file, out=out_file)
            stdout, stderr = command()
            self.assertFalse(stderr)
            self.assertFalse(stdout)
            self.assertTrue(os.path.exists(out_file))
        finally:
            if os.path.exists(out_file):
                os.remove(out_file)

    def test_igblastp_parse(self):
        query_file = os.path.join(os.path.dirname(__file__), '..', 'test_bio', 'resources', 'igblastp.fasta')
        with open(query_file, 'r') as fh:
            query = SeqIO.read(fh, 'fasta')
        outfile = os.path.join(os.path.dirname(__file__), '..', 'test_bio', 'resources', 'igblastp.out')
        igresult: IgResult = parse_igblastp_results(outfile)
        alignments: List[SeqRecord] = build_common_alignments(query, igresult.hits)
        self.assertEqual(6, len(igresult.regions))
        self.assertEqual(3, len(igresult.hits))
        self.assertEqual(4, len(alignments))

    def test_igblastn_parse(self):
        query_file = os.path.join(os.path.dirname(__file__), '..', 'test_bio', 'resources', 'Y14934.gb')
        with open(query_file, 'r') as fh:
            query = SeqIO.read(fh, 'gb')
        outfile = os.path.join(os.path.dirname(__file__), '..', 'test_bio', 'resources', 'igblastn.out')
        igresult: IgResult = parse_igblastn_results(outfile)
        alignments: List[SeqRecord] = build_common_alignments(query, igresult.hits)
        self.assertEqual(6, len(igresult.regions))
        self.assertEqual(9, len(igresult.hits))
        self.assertEqual(10, len(alignments))


if __name__ == '__main__':
    main()
