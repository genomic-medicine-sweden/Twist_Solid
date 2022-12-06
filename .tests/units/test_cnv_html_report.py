import json
import os
import pytest
import sys

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPT_DIR = os.path.abspath(os.path.join(TEST_DIR, "../../workflow/scripts"))
sys.path.insert(0, SCRIPT_DIR)

import cnv_json  # noqa


gatk_ratios = [
    "@HD\tVN=1.6\n",
    "@SQ\tSN:chr1\tLN:249250621\n",
    "@RG\tID:GATKCopyNumber\tSM:test\n",
    "CONTIG\tSTART\tEND\tLOG2_COPY_RATIO\n",
    "chr1\t1\t10\t1.2\n",
    "chr1\t10\t20\t1.1\n",
    "chr1\t20\t30\t3.4\n",
]

cnvkit_ratios = [
    "chromosome\tstart\tend\tgene\tdepth\tlog2\tweight\n",
    "chr1\t150500\t307183\tAntitarget\t1.5853\t0.504196\t0.948424\n",
    "chr1\t307183\t463867\tAntitarget\t13.1414\t0.494873\t0.847039\n",
    "chr1\t463867\t620550\tAntitarget\t0.841565\t0.722282\t0.637809\n",
]

cnvkit_segments = [
    "chromosome\tstart\tend\tgene\tlog2\tdepth\tprobes\tweight\tci_lo\tci_hi\n",
    "chr1\t150500\t61370310\tHES5_Exon2_NM_001010926.4,HES5_Exon1_NM_001010926.4\t-0.168018\t"
    "876.26\t862\t825.95\t-0.171634\t-0.163029\n",
    "chr1\t61370310\t62653728\t-\t4.3806\t22.4843\t8\t7.83195\t4.25645\t4.47999\n",
    "chr1\t62653728\t121635711\tJAK1_Exon26_NM_001320923.1,JAK1_Exon25_NM_001320923.1\t-0.0693903\t"
    "485.77\t564\t539.866\t-0.07\t48381\t-0.0629965\n",
    "chr1\t142388707\t188892112\tMCL1_Exon2_NM_182763.2,MCL1_Exon3_NM_001197320.1\t0.0677003\t"
    "584.127\t456\t435.344\t0.0627169\t0.0738394\n",
]


def test_parse_gatk_ratios(tmp_path):
    f = tmp_path / "gatk_ratios.txt"
    f.write_text("".join(gatk_ratios))

    ratios = cnv_json.PARSERS["gatk"]["ratios"](f)

    assert ratios[0]["chromosome"] == "chr1"
    assert ratios[0]["start"] == 1
    assert ratios[0]["end"] == 10
    assert ratios[0]["log2"] == pytest.approx(1.2)


def test_parse_cnvkit_ratios(tmp_path):
    f = tmp_path / "cnvkit_ratios.txt"
    f.write_text("".join(cnvkit_ratios))

    ratios = cnv_json.PARSERS["cnvkit"]["ratios"](f)

    assert ratios[0]["chromosome"] == "chr1"
    assert ratios[0]["start"] == 150500
    assert ratios[0]["end"] == 307183
    assert ratios[0]["log2"] == pytest.approx(0.504196)


def test_parse_cnvkit_segments(tmp_path):
    f = tmp_path / "cnvkit_segments.txt"
    f.write_text("".join(cnvkit_segments))

    ratios = cnv_json.PARSERS["cnvkit"]["segments"](f)

    assert ratios[0]["chromosome"] == "chr1"
    assert ratios[0]["start"] == 150500
    assert ratios[0]["end"] == 61370310
    assert ratios[0]["log2"] == pytest.approx(-0.168018)


def test_cnv_json(tmp_path):
    ratio_file = tmp_path / "cnvkit_ratios.txt"
    ratio_file.write_text("".join(cnvkit_ratios))
    segment_file = tmp_path / "cnvkit_segments.txt"
    segment_file.write_text("".join(cnvkit_segments))

    ratios = cnv_json.PARSERS["cnvkit"]["ratios"](ratio_file)
    segments = cnv_json.PARSERS["cnvkit"]["segments"](segment_file)

    json_str = cnv_json.to_json("cnvkit", ratios, segments)
    json_dict = json.loads(json_str)

    assert "caller" in json_dict
    assert "ratios" in json_dict
    assert "segments" in json_dict
    assert json_dict["caller"] == "cnvkit"
    assert len(json_dict["ratios"]) == 3
    assert len(json_dict["segments"]) == 4
