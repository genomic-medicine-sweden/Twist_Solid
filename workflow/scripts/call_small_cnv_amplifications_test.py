import tempfile
import os
import unittest


class TestUnitUtils(unittest.TestCase):
    def setUp(self):

        self.window_size = 3
        self.region_max_size = 5
        self.min_nr_stdev_diff = 8
        self.min_log_odds_diff = 0.4

        self.tempdir = tempfile.mkdtemp()

    def tearDown(self):
        pass

    def _test_call_small_cnv_amplifications_ok(self, test_results, variants):
        i = 0
        for variant in variants:
            columns = variant.strip().split("\t")
            try:
                self.assertEqual(test_results[i], columns[i])
            except AssertionError as e:
                print("Failed CNV amplification calling of: " + str(variant))
                raise e
            i += 1

    def _test_call_small_cnv_amplifications_filter(self, test_filter, filter):
        try:
            self.assertEqual(test_filter, filter)
        except AssertionError as e:
            print("Failed filter: " + str(test_filter))
            raise e

    def test_call_small_cnv_amplifications_ok(self):
        from call_small_cnv_amplifications import call_small_cnv_amplifications

        # True call
        cnv_data1 = ".tests/units/gatk_cnv/HD832.HES45_T.test2.tsv"
        regions_file = open(".tests/integration/DATA/cnv_amplification_genes.tsv")
        out_amplifications = open(os.path.join(self.tempdir, "HD832.HES45_T.amplifications.tsv"), "w")

        filter = call_small_cnv_amplifications(
            cnv_data1, regions_file, out_amplifications, self.window_size, self.region_max_size, self.min_nr_stdev_diff,
            self.min_log_odds_diff,
        )

        out_amplifications.close()
        regions_file.close()
        result_file = open(os.path.join(self.tempdir, "HD832.HES45_T.amplifications.tsv"))

        result = []
        next(result_file)
        for line in result_file:
            result.append(line)
            print(line)

        test_results = "CDKN2A,CDKN2B\t-0.53904725\t-0.785199\t-0.24615174999999997\t4\t-8.785756385043994\n"

        self._test_call_small_cnv_amplifications_ok(test_results.split("\t"), result)

    def test_call_small_cnv_amplifications_filter(self):
        from call_small_cnv_amplifications import call_small_cnv_amplifications

        # Filter 1 (Deletions)
        cnv_data = ".tests/units/gatk_cnv/HD832.HES45_T.filter1.tsv"
        regions_file = open(".tests/integration/DATA/cnv_amplification_genes.tsv")
        out_amplifications = open(os.path.join(self.tempdir, "HD832.HES45_T.amplifications.tsv"), "w")
        filter = call_small_cnv_amplifications(
            cnv_data, regions_file, out_amplifications, self.window_size, self.region_max_size, self.min_nr_stdev_diff,
            self.min_log_odds_diff,
        )
        self._test_call_small_cnv_amplifications_filter("Deletion", filter)
        out_amplifications.close()
        regions_file.close()

        # Filter 2 (Large amplifications)
        cnv_data = ".tests/units/gatk_cnv/HD832.HES45_T.filter2_amp.tsv"
        regions_file = open(".tests/integration/DATA/cnv_amplification_genes.tsv")
        out_amplifications = open(os.path.join(self.tempdir, "HD832.HES45_T.amplifications.tsv"), "w")
        filter = call_small_cnv_amplifications(
            cnv_data, regions_file, out_amplifications, self.window_size, self.region_max_size, self.min_nr_stdev_diff,
            self.min_log_odds_diff,
        )
        self._test_call_small_cnv_amplifications_filter("Too_large", filter)
        out_amplifications.close()
        regions_file.close()

        # Filter 4 (Filter amplifications not in the actual gene of interest)
        cnv_data = ".tests/units/gatk_cnv/HD832.HES45_T.filter4_amp.tsv"
        regions_file = open(".tests/integration/DATA/cnv_amplification_genes.tsv")
        out_amplifications = open(os.path.join(self.tempdir, "HD832.HES45_T.amplifications.tsv"), "w")
        filter = call_small_cnv_amplifications(
            cnv_data, regions_file, out_amplifications, self.window_size, self.region_max_size, self.min_nr_stdev_diff,
            self.min_log_odds_diff,
        )
        self._test_call_small_cnv_amplifications_filter("Not_in_gene", filter)
        out_amplifications.close()
        regions_file.close()

        # Filter 5 (Filter too short amplifications)
        self.window_size = 5
        cnv_data = ".tests/units/gatk_cnv/HD832.HES45_T.filter5_amp.tsv"
        regions_file = open(".tests/integration/DATA/cnv_amplification_genes.tsv")
        out_amplifications = open(os.path.join(self.tempdir, "HD832.HES45_T.amplifications.tsv"), "w")
        filter = call_small_cnv_amplifications(
            cnv_data, regions_file, out_amplifications, self.window_size, self.region_max_size, self.min_nr_stdev_diff,
            self.min_log_odds_diff,
        )
        self._test_call_small_cnv_amplifications_filter("Too_small", filter)
        out_amplifications.close()
        regions_file.close()

        # Filter 6 (Filter amplifications with positive median)
        self.window_size = 3
        cnv_data = ".tests/units/gatk_cnv/HD832.HES45_T.filter6_amp.tsv"
        regions_file = open(".tests/integration/DATA/cnv_amplification_genes.tsv")
        out_amplifications = open(os.path.join(self.tempdir, "HD832.HES45_T.amplifications.tsv"), "w")
        filter = call_small_cnv_amplifications(
            cnv_data, regions_file, out_amplifications, self.window_size, self.region_max_size, self.min_nr_stdev_diff,
            self.min_log_odds_diff,
        )
        self._test_call_small_cnv_amplifications_filter("Not_amplification", filter)
        out_amplifications.close()
        regions_file.close()

        # Filter 7 (Too low difference in log odds ratio)
        cnv_data = ".tests/units/gatk_cnv/HD832.HES45_T.filter7_amp.tsv"
        regions_file = open(".tests/integration/DATA/cnv_amplification_genes.tsv")
        out_amplifications = open(os.path.join(self.tempdir, "HD832.HES45_T.amplifications.tsv"), "w")
        filter = call_small_cnv_amplifications(
            cnv_data, regions_file, out_amplifications, self.window_size, self.region_max_size, self.min_nr_stdev_diff,
            self.min_log_odds_diff,
        )
        self._test_call_small_cnv_amplifications_filter("Low_abs_diff", filter)
        out_amplifications.close()
        regions_file.close()

        # Filter 8 (Too few standard deviations in difference)
        self.region_max_size = 6
        cnv_data = ".tests/units/gatk_cnv/HD832.HES45_T.filter8_amp.tsv"
        regions_file = open(".tests/integration/DATA/cnv_amplification_genes.tsv")
        out_amplifications = open(os.path.join(self.tempdir, "HD832.HES45_T.amplifications.tsv"), "w")
        filter = call_small_cnv_amplifications(
            cnv_data, regions_file, out_amplifications, self.window_size, self.region_max_size, 100,
            self.min_log_odds_diff,
        )
        self._test_call_small_cnv_amplifications_filter("Low_nr_std_diff", filter)
        out_amplifications.close()
        regions_file.close()

        # Filter 9 (Filter amplifications with too few data points outside to calculate median and standard deviations)
        cnv_data = ".tests/units/gatk_cnv/HD832.HES45_T.filter9_amp.tsv"
        regions_file = open(".tests/integration/DATA/cnv_amplification_genes.tsv")
        out_amplifications = open(os.path.join(self.tempdir, "HD832.HES45_T.amplifications.tsv"), "w")
        filter = call_small_cnv_amplifications(
            cnv_data, regions_file, out_amplifications, self.window_size, self.region_max_size, self.min_nr_stdev_diff,
            self.min_log_odds_diff,
        )
        self._test_call_small_cnv_amplifications_filter("Too_few_outside", filter)
        out_amplifications.close()
        regions_file.close()
