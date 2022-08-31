import os
from unittest import TestCase, main

from df.MultipleSequenceAlignment import MultipleSequenceAlignment
from df.data_transfer import DataFunctionRequest
from test_data_functions import clean_output_files


class TestMultipleSequenceAlignment(TestCase):

    def test_multiple_sequence_alignment1(self):
        json_file = os.path.join(os.path.dirname(__file__), 'resources', 'msa_error1.json')
        with open(json_file) as fh:
            request_json = fh.read()
        request = DataFunctionRequest.parse_raw(request_json)
        self.assertEqual(1, len(request.inputColumns))
        input_values = next(iter(request.inputColumns.values())).values
        self.assertEqual(84, len(input_values))
        self.assertIsNone(input_values[63])
        self.assertIsNone(input_values[74])
        test = [t for t in input_values if t]
        self.assertEqual(82, len(test))
        df = MultipleSequenceAlignment()
        response = df.execute(request)
        output_columns = response.outputColumns
        self.assertEqual(3, len(output_columns))
        for col in output_columns:
            output_values = col.values
            self.assertEqual(84, len(output_values))
            self.assertIsNone(output_values[63])
            self.assertIsNone(output_values[74])
            test = [t for t in output_values if t]
            self.assertEqual(82, len(test))

    @classmethod
    def tearDownClass(cls) -> None:
        clean_output_files()


if __name__ == "__main__":
    main()
