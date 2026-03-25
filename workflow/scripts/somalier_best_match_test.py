
import unittest
import pandas as pd
import tempfile
import os
import sys
from somalier_best_match import main


class TestSomalierBestMatch(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_dna_rna_match_logic(self):
        # Create a mock pairs file
        # DNA: _T, _N
        # RNA: _R
        input_data = {
            '#sample_a': ['P1_T', 'P1_T', 'P2_T', 'P1_R'],
            'sample_b':  ['P2_T', 'P1_R', 'P2_R', 'P2_R'],
            'relatedness': [0.99, 0.95, 0.94, 0.1],
            'ibs0': [0, 1, 1, 50],
            'ibs2': [200, 190, 195, 10],
            'n': [200, 200, 200, 200]
        }

        input_df = pd.DataFrame(input_data)
        input_path = os.path.join(self.temp_dir.name, "input.tsv")
        input_df.to_csv(input_path, sep='\t', index=False)

        output_path = os.path.join(self.temp_dir.name, "output.tsv")

        # Test with default cutoff (0.7)
        main(input_path, output_path, 0.7)

        output_df = pd.read_csv(output_path, sep='\t')

        # Only RNA samples should be in the 'Sample' column
        samples = output_df['Sample'].unique()
        for sample in samples:
            self.assertTrue(sample.endswith('_R'))

        # P1_R should match P1_T (0.95)
        p1_r_match = output_df[output_df['Sample'] == 'P1_R']
        self.assertEqual(p1_r_match['Best_Match'].values[0], 'P1_T')
        self.assertEqual(p1_r_match['Relatedness'].values[0], 0.95)
        self.assertEqual(p1_r_match['Match'].values[0], 'yes')
        self.assertEqual(p1_r_match['Variants_compared'].values[0], 200)

        # Test with a higher cutoff that would make it 'no'
        main(input_path, output_path, 0.96)
        output_df_high = pd.read_csv(output_path, sep='\t')
        p1_r_match_high = output_df_high[output_df_high['Sample'] == 'P1_R']
        self.assertEqual(p1_r_match_high['Match'].values[0], 'no')

    def test_empty_input(self):
        input_path = os.path.join(self.temp_dir.name, "empty.tsv")
        pd.DataFrame(
            columns=['#sample_a', 'sample_b', 'relatedness', 'ibs0', 'ibs2', 'n']
        ).to_csv(input_path, sep='\t', index=False)

        output_path = os.path.join(self.temp_dir.name, "output_empty.tsv")
        main(input_path, output_path, 0.7)

        output_df = pd.read_csv(output_path, sep='\t')
        self.assertTrue(output_df.empty)
        # Check header of empty output
        self.assertIn('Variants_compared', output_df.columns)
        self.assertIn('Match', output_df.columns)


if __name__ == '__main__':
    unittest.main()
