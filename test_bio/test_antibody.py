import os
import re
from typing import Optional, List
from unittest import TestCase, main, skip

from Bio import Seq, SeqIO
from Bio.SeqRecord import SeqRecord

from ruse.bio.antibody import label_antibody_sequences, align_antibody_sequences, AntibodyAlignmentResult


class TestAntibody(TestCase):

    def test_label_antibody_sequences(self) -> None:
        seq_file = os.path.join(os.path.dirname(__file__), 'resources/antibody_sequences.fasta')
        antibody_sequences = list(SeqIO.parse(seq_file, 'fasta'))[10:15]
        label_antibody_sequences(antibody_sequences, 'chothia')

    def test_process_multiple_domains(self):
        seq_file = os.path.join(os.path.dirname(__file__), 'resources/antibody_sequences.fasta')
        antibody_sequences = list(s for s in SeqIO.parse(seq_file, 'fasta') if s.id == '3h3b_C')
        self.align_antibody_sequences(antibody_sequences, 'chothia')

    def test_process_tricky(self) -> None:
        seq_file = os.path.join(os.path.dirname(__file__), 'resources/antibody_sequences.fasta')
        antibody_sequences = list(s for s in SeqIO.parse(seq_file, 'fasta') if s.id == '4hjj_L')
        self.align_antibody_sequences(antibody_sequences, 'chothia')

    def test_failed(self) -> None:
        # only the first of these guys is processed by anarci
        ids = {'1mhp_Y', '4od3_H', '4od1_H', '4odh_H', '4ocs_H'}
        seq_file = os.path.join(os.path.dirname(__file__), 'resources/antibody_sequences.fasta')
        antibody_sequences = list(s for s in SeqIO.parse(seq_file, 'fasta') if s.id in ids)
        results = self.align_antibody_sequences(antibody_sequences, 'chothia', True)
        for seq in results.aligned_sequences[1:]:
            self.assertTrue(seq is None)

    def test_not_antibodies(self) -> None:
        seq_file = os.path.join(os.path.dirname(__file__), 'resources/phylodrome.fasta')
        sequences = list(SeqIO.parse(seq_file, 'fasta'))
        align_information = align_antibody_sequences(sequences, 'chothia')
        self.assertEqual(align_information.numbering, [])
        self.assertEqual(align_information.regions, [])

    def test_process_tricky2(self) -> None:
        def find_position(pos):
            return next((n for n in results.numbering if n.query_position == pos), None)

        ids = {'1cly_H', '1cly_L', '1clz_H', '1clz_L', '1clo_H', '1clo_L', '1cl7_H', '1cl7_L', '1a7o_H', '1a7o_L',
               '1a7n_H', '1a7n_L', '1a7r_H', '1a7r_L', '1a7q_H', '1a7q_L', '1a7p_H', '1a7p_L', '12e8_H',
               '12e8_M', '1bvk_A', '1bvk_B', '1bvl_A', '1bvl_C', '1a4k_A', '1a4k_H', '1bgx_H', '1bgx_L',
               '1a5f_H', '1a5f_    L', '1cr9_H', '1cr9_L', '1bfv_H', '1bfv_L', '1bfo_A', '1bfo_C', '1dba_H',
               '1dba_L', '1dbb_H', '1dbb_L', '1dbj_H', '1a2y_A', '1a2y_B', '1bey_H'}
        seq_file = os.path.join(os.path.dirname(__file__), 'resources/antibody_sequences.fasta')
        antibody_sequences = list(s for s in SeqIO.parse(seq_file, 'fasta') if s.id in ids)
        results = self.align_antibody_sequences(antibody_sequences, 'kabat')
        heavy_region = results.regions[1]

        cdr1_s = find_position(heavy_region.cdr1_start.query_position)
        cdr1_e = find_position(heavy_region.cdr1_end.query_position)
        cdr2_s = find_position(heavy_region.cdr2_start.query_position)
        cdr2_e = find_position(heavy_region.cdr2_end.query_position)
        cdr3_s = find_position(heavy_region.cdr3_start.query_position)
        cdr3_e = find_position(heavy_region.cdr3_end.query_position)

        self.assertEqual(cdr1_s.label(), "H31")
        self.assertEqual(cdr1_e.label(), "H35A")
        self.assertEqual(cdr2_s.label(), "H50")
        self.assertEqual(cdr2_e.label(), "H65")
        self.assertEqual(cdr3_s.label(), "H95")
        self.assertEqual(cdr3_e.label(), "H102")

    def test_some_antibody_sequences_kabat(self):
        self.align_example_antibody_sequences(100, 'kabat')

    def test_some_antibody_sequences_chothia(self):
        self.align_example_antibody_sequences(100, 'chothia')

    @skip("Skipping time consuming test of ~1000 antibodies")
    def test_all_antibody_sequences(self):
        self.align_example_antibody_sequences(None, 'chothia')

    def align_example_antibody_sequences(self, count: Optional[int], scheme: Optional[str]) -> None:
        seq_file = os.path.join(os.path.dirname(__file__), 'resources/antibody_sequences.fasta')
        antibody_sequences = list(SeqIO.parse(seq_file, 'fasta'))
        if count:
            antibody_sequences = antibody_sequences[0:count]
        if not scheme:
            scheme = 'chothia'
        self.align_antibody_sequences(antibody_sequences, scheme)

    def align_antibody_sequences(self, antibody_sequences: List[SeqRecord], scheme: str,
                                 bad_input=False) -> AntibodyAlignmentResult:
        align_information = align_antibody_sequences(antibody_sequences, scheme)
        numbering = align_information.numbering
        aligned_sequences = align_information.aligned_sequences

        # thses guys all fail under anarci
        skip_ids = {'4od3_H', '4od1_H', '4odh_H', '4ocs_H'}
        matcher = re.compile(r"^antibody_number: (\w+)$")
        all_matches = set()
        process_count = 0
        for sequence in aligned_sequences:
            cnt = 0
            if sequence:
                for feature in sequence.features:
                    if 'note' in feature.qualifiers and isinstance(feature.qualifiers['note'], list):
                        match = matcher.match(feature.qualifiers['note'][0])
                        if match:
                            cnt += 1
                            pos = feature.location.start
                            all_matches.add(pos)
                            numbering_match = next((n for n in numbering if n.query_position == pos), None)
                            self.assertTrue(numbering_match)
                            self.assertEqual(len(feature.location), 1)
                            label = match.group(1)
                            self.assertEqual(numbering_match.label(), label)

                if sequence.id not in skip_ids:
                    self.assertTrue(cnt > 10)
                    process_count += 1

        # 4 or the first 5 sequences fail under anarci and we get 6 missing light chain positions so don't test on those
        missing = [n for n in numbering if n.query_position not in all_matches]
        if bad_input:
            self.assertLessEqual(len(missing), 6)
        else:
            self.assertLessEqual(len(missing), 0)

        return align_information

    # @skip("Abysis labelling is not reproducible over ANARCI versions")
    def test_abysis_small(self) -> None:
        def feature_location(name):
            feature = next(f for f in features if f.qualifiers['note'][0] == 'antibody_label: {}'.format(name))
            return feature.location.start, feature.location.end

        seq_file = os.path.join(os.path.dirname(__file__), 'resources/abysis-small.gb')
        antibody_sequences = list(SeqIO.parse(seq_file, 'gb'))
        labelled_sequences = [SeqRecord(s.seq, id=s.id, name=s.name, description=s.description,
                                        annotations={'molecule_type': 'protein'}) for s in antibody_sequences]
        label_antibody_sequences(labelled_sequences, 'chothia')

        features = labelled_sequences[0].features
        start, end = feature_location('HFR1')
        self.assertEqual(start, 0)
        self.assertEqual(end, 30)
        start, end = feature_location('CDR-H1')
        self.assertEqual(start, 30)
        self.assertEqual(end, 35)
        start, end = feature_location('HFR2')
        self.assertEqual(start, 35)
        self.assertEqual(end, 49)
        start, end = feature_location('CDR-H2')
        self.assertEqual(start, 49)
        self.assertEqual(end, 66)
        start, end = feature_location('HFR3')
        self.assertEqual(start, 66)
        self.assertEqual(end, 98)
        start, end = feature_location('CDR-H3')
        self.assertEqual(start, 98)
        self.assertEqual(end, 108)
        start, end = feature_location('HFR4')
        self.assertEqual(start, 108)
        self.assertEqual(end, 119)

        out_file = os.path.join(os.path.dirname(__file__), 'output', 'abysis-small.gb')
        with open(out_file, 'w') as fh:
            SeqIO.write(labelled_sequences, fh, 'gb')


if __name__ == "__main__":
    main()
