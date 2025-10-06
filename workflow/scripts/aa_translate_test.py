import tempfile
import os
import unittest


class TestUnitUtils(unittest.TestCase):
    def setUp(self):

        self.tempdir = tempfile.mkdtemp()

    def tearDown(self):
        pass

    def _test_translate_file_ok(self, test_results, results):
        i = 0
        for res_line in results:
            res_columns = res_line.split("\t")
            test_res_columns = test_results[i].split("\t")
            j = 0
            for res_column in res_columns:
                try:
                    self.assertEqual(res_columns[j], test_res_columns[j])
                except AssertionError as e:
                    print("Failed aa_translation of: " + str(res_columns[j]))
                    raise e
                j += 1
            i += 1

    def test_translate_file(self):
        from aa_translate import translate_file

        # True call
        input_file = open(".tests/units/cov_and_mut.tsv")
        output_file = open(os.path.join(self.tempdir, "cov_and_mut.aa_translate.tsv"), "w")
        translate_file(input_file, output_file)

        output_file.close()
        input_file.close()
        result_file = open(os.path.join(self.tempdir, "cov_and_mut.aa_translate.tsv"))

        result = []
        for line in result_file:
            result.append(line.strip())

        test_results = [
            "chr	start	stop	Gene	Consequence	Exon	AA_change	CDS_change	Accession_number	report	Comment	Min_read_depth200	gvcf_depth	ref	alt	DP	ref_depth	alt_depth	AF	rs	db1000G	GnomAD	clinvar	background_median	nr_std_from_background_median	nr_observed_in_normals	Accession_number_VEP",  # noqa
            "NC_000001.10	16464612	16464612	EPHA2	missense_variant	-	P350T	NM_004431.5:c.1048C>A	NM_004431.5	4-other	-	ok	1099	G	T	1829	1015	808	0.442	rs11543934	0.0004	0.003144	benign	0.00070	377.4	0,0,53	NM_004431.5",  # noqa
            "NC_000001.10	16464612	16464612	EPHA2	missense_variant	-	-	NM_004431.5:c.1048C>A	NM_004431.5	4-other	-	ok	1099	G	T	1829	1015	808	0.442	rs11543934	0.0004	0.003144	benign	0.00070	377.4	0,0,53	NM_004431.5",  # noqa
            "NC_000001.10	16464612	16464612	EPHA2	missense_variant	-	P350*	NM_004431.5:c.1048C>A	NM_004431.5	4-other	-	ok	1099	G	T	1829	1015	808	0.442	rs11543934	0.0004	0.003144	benign	0.00070	377.4	0,0,53	NM_004431.5",  # noqa
            "NC_000001.10	16464612	16464612	EPHA2	missense_variant	-	hmmm	NM_004431.5:c.1048C>A	NM_004431.5	4-other	-	ok	1099	G	T	1829	1015	808	0.442	rs11543934	0.0004	0.003144	benign	0.00070	377.4	0,0,53	NM_004431.5",  # noqa
            "NC_000001.10	16464612	16464612	EPHA2	missense_variant	-	P3500	NM_004431.5:c.1048C>A	NM_004431.5	4-other	-	ok	1099	G	T	1829	1015	808	0.442	rs11543934	0.0004	0.003144	benign	0.00070	377.4	0,0,53	NM_004431.5",  # noqa
            "NC_000001.10	16464612	16464612	EPHA2	missense_variant	-	P350	NM_004431.5:c.1048C>A	NM_004431.5	4-other	-	ok	1099	G	T	1829	1015	808	0.442	rs11543934	0.0004	0.003144	benign	0.00070	377.4	0,0,53	NM_004431.5",  # noqa
            "NC_000001.10	16464612	16464612	EPHA2	missense_variant	-	P35	NM_004431.5:c.1048C>A	NM_004431.5	4-other	-	ok	1099	G	T	1829	1015	808	0.442	rs11543934	0.0004	0.003144	benign	0.00070	377.4	0,0,53	NM_004431.5",  # noqa
            "NC_000001.10	16464612	16464612	EPHA2	missense_variant	-	P3A	NM_004431.5:c.1048C>A	NM_004431.5	4-other	-	ok	1099	G	T	1829	1015	808	0.442	rs11543934	0.0004	0.003144	benign	0.00070	377.4	0,0,53	NM_004431.5",  # noqa
            "NC_000001.10	16464612	16464612	EPHA2	missense_variant	-	A35H	NM_004431.5:c.1048C>A	NM_004431.5	4-other	-	ok	1099	G	T	1829	1015	808	0.442	rs11543934	0.0004	0.003144	benign	0.00070	377.4	0,0,53	NM_004431.5",  # noqa
            "NC_000001.10	16464612	16464612	EPHA2	missense_variant	-	P350Afs*40	NM_004431.5:c.1048C>A	NM_004431.5	4-other	-	ok	1099	G	T	1829	1015	808	0.442	rs11543934	0.0004	0.003144	benign	0.00070	377.4	0,0,53	NM_004431.5",  # noqa
            "NC_000001.10	16464612	16464612	EPHA2	missense_variant	-	P350Afs*	NM_004431.5:c.1048C>A	NM_004431.5	4-other	-	ok	1099	G	T	1829	1015	808	0.442	rs11543934	0.0004	0.003144	benign	0.00070	377.4	0,0,53	NM_004431.5",  # noqa
            "NC_000001.10	16464612	16464612	EPHA2	missense_variant	-	C184_T185insA	NM_004431.5:c.1048C>A	NM_004431.5	4-other	-	ok	1099	G	T	1829	1015	808	0.442	rs11543934	0.0004	0.003144	benign	0.00070	377.4	0,0,53	NM_004431.5",  # noqa
            "NC_000001.10	16464612	16464612	EPHA2	missense_variant	-	K319_T320insPVARIGQTGTPSVFSQRGKSRK	NM_004431.5:c.1048C>A	NM_004431.5	4-other	-	ok	1099	G	T	1829	1015	808	0.442	rs11543934	0.0004	0.003144	benign	0.00070	377.4	0,0,53	NM_004431.5",  # noqa
            "NC_000001.10	16464612	16464612	EPHA2	missense_variant	-	Q534_Q550del	NM_004431.5:c.1048C>A	NM_004431.5	4-other	-	ok	1099	G	T	1829	1015	808	0.442	rs11543934	0.0004	0.003144	benign	0.00070	377.4	0,0,53	NM_004431.5",  # noqa
            "NC_000001.10	16464612	16464612	EPHA2	missense_variant	-	Q5340del	NM_004431.5:c.1048C>A	NM_004431.5	4-other	-	ok	1099	G	T	1829	1015	808	0.442	rs11543934	0.0004	0.003144	benign	0.00070	377.4	0,0,53	NM_004431.5",  # noqa
            "NC_000001.10	16464612	16464612	EPHA2	missense_variant	-	T108delinsGP	NM_004431.5:c.1048C>A	NM_004431.5	4-other	-	ok	1099	G	T	1829	1015	808	0.442	rs11543934	0.0004	0.003144	benign	0.00070	377.4	0,0,53	NM_004431.5",  # noqa
            "NC_000001.10	16464612	16464612	EPHA2	missense_variant	-	K169_S192delins*	NM_004431.5:c.1048C>A	NM_004431.5	4-other	-	ok	1099	G	T	1829	1015	808	0.442	rs11543934	0.0004	0.003144	benign	0.00070	377.4	0,0,53	NM_004431.5",  # noqa
            "NC_000001.10	16464612	16464612	EPHA2	missense_variant	-	D988_M1015dup	NM_004431.5:c.1048C>A	NM_004431.5	4-other	-	ok	1099	G	T	1829	1015	808	0.442	rs11543934	0.0004	0.003144	benign	0.00070	377.4	0,0,53	NM_004431.5",  # noqa
            "NC_000001.10	16464612	16464612	EPHA2	missense_variant	-	D988dup	NM_004431.5:c.1048C>A	NM_004431.5	4-other	-	ok	1099	G	T	1829	1015	808	0.442	rs11543934	0.0004	0.003144	benign	0.00070	377.4	0,0,53	NM_004431.5",  # noqa
            "NC_000001.10	16464612	16464612	EPHA2	missense_variant	-	*541Kext*65	NM_004431.5:c.1048C>A	NM_004431.5	4-other	-	ok	1099	G	T	1829	1015	808	0.442	rs11543934	0.0004	0.003144	benign	0.00070	377.4	0,0,53	NM_004431.5",  # noqa
        ]

        self._test_translate_file_ok(test_results, result)
