import os
from unittest import TestCase, main, skip

from Bio import SeqIO

from ruse.bio.enzyme_cleavage import enzyme_cleavage_sites


class EnzymeCleavageTestCase(TestCase):

    @skip('No longer using EMBOSS restrict')
    def test_enzyme_alignment(self):
        input_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "resources/example_nt.gb"))
        sequences = list(SeqIO.parse(input_file, 'gb'))[:10]

        out_sequences = enzyme_cleavage_sites(sequences)

        output_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "output/test_enzyme_cleavage_nt.gb"))
        SeqIO.write(out_sequences, output_file, 'gb')


if __name__ == "__main__":
    main()
