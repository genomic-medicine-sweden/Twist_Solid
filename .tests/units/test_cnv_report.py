import os
import sys
import tempfile
import unittest
from dataclasses import dataclass

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPT_DIR = os.path.abspath(os.path.join(TEST_DIR, "../../workflow/scripts"))
sys.path.insert(0, SCRIPT_DIR)

from cnv_report import create_tsv_report  # noqa


class TestGetCaller(unittest.TestCase):
    def test_create_tsv_report(self):
        cnv = tempfile.mkdtemp() + "/tcnv.txt"
        out_additional_only = open(tempfile.mkdtemp() + "/out_additional_only.txt", "w")
        create_tsv_report(
            [".tests/units/vcf/test.cnv1.vcf"],
            [".tests/units/vcf/test.cnv2.vcf"],
            ".tests/units/vcf/test.deletions.tsv",
            ".tests/units/vcf/test.amplifications.tsv",
            2.0,
            cnv,
            out_additional_only,
            1.5,
            0.5,
            0.5,
        )

        @dataclass
        class TestCase:
            name: str
            expected: str

        testcases = [
                TestCase(
                    name="header row",
                    expected=("sample", "gene(s)", "chrom", "region", "callers", "freq_in_db", "copy_number")
                ),
                TestCase(
                    name="variant 1",
                    expected=("testSample_T", "FGFR1", "chr8", "34370199-43930232", "cnvkit", "0.01", "8.59")
                ),
                TestCase(
                    name="variant 2",
                    expected=("testSample_T", "FGFR1,MYC", "chr8", "35008818-146144253", "gatk_cnv", "0.01", "7.01")
                ),
                TestCase(
                    name="variant 3",
                    expected=("testSample_T", "MYC", "chr8", "46689525-146144003", "cnvkit", "0.09", "5.06")
                ),
                TestCase(
                    name="small deletion",
                    expected=("testSample_T", "CDKN2A,CDKN2B", "chr9", "21968207-22008972", "small_deletion", "NA", "-0.28")
                ),
                TestCase(
                    name="small amplification",
                    expected=("testSample_T", "MYCN", "chr2", "16968207-17008972", "small_amplification", "NA", "7.29")
                ),
        ]

        test_it = iter(testcases)
        with open(cnv) as reader:
            for line in reader:

                actual = tuple(line.rstrip().split("\t"))
                try:
                    case = next(test_it)
                    self.assertEqual(
                        case.expected,
                        actual,
                        "failed test '{}': expected {}, got {}".format(
                            case.name, case.expected, actual
                        ),
                    )
                except ValueError:
                    if case.expected == "ValueError":
                        assert True
                    else:
                        assert False
                except StopIteration:
                    print("More results than test cases")
                    assert False
