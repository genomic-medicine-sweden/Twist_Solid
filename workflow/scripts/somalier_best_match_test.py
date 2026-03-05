
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
        # P1_T (DNA) has:
        # P2_T (DNA) with 0.99 -> IGNORED (same type)
        # P1_R (RNA) with 0.95 -> BEST MATCH
        
        # P1_R (RNA) has:
        # P1_T (DNA) with 0.95 -> BEST MATCH
        # P2_R (RNA) with 0.1 -> IGNORED
        
        input_df = pd.DataFrame(input_data)
        input_path = os.path.join(self.temp_dir.name, "input.tsv")
        input_df.to_csv(input_path, sep='\t', index=False)
        
        output_path = os.path.join(self.temp_dir.name, "output.tsv")
        
        main(input_path, output_path)
        
        output_df = pd.read_csv(output_path, sep='\t')
        
        # P1_T should match P1_R (0.95)
        self.assertEqual(output_df[output_df['Sample'] == 'P1_T']['Best_Match'].values[0], 'P1_R')
        self.assertEqual(output_df[output_df['Sample'] == 'P1_T']['Relatedness'].values[0], 0.95)
        
        # P1_R should match P1_T (0.95)
        self.assertEqual(output_df[output_df['Sample'] == 'P1_R']['Best_Match'].values[0], 'P1_T')
        
        # P2_T should match P2_R (0.94)
        self.assertEqual(output_df[output_df['Sample'] == 'P2_T']['Best_Match'].values[0], 'P2_R')

    def test_empty_input(self):
        input_path = os.path.join(self.temp_dir.name, "empty.tsv")
        pd.DataFrame(columns=['#sample_a', 'sample_b', 'relatedness', 'ibs0', 'ibs2', 'n']).to_csv(input_path, sep='\t', index=False)
        
        output_path = os.path.join(self.temp_dir.name, "output_empty.tsv")
        main(input_path, output_path)
        
        output_df = pd.read_csv(output_path, sep='\t')
        self.assertTrue(output_df.empty)

if __name__ == '__main__':
    unittest.main()
