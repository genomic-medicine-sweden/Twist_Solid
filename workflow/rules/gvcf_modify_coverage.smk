# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2022, Jonas Almlöf"
__email__ = "jonas.almlof@igp.uu.se"
__license__ = "GPL-3"


rule gvcf_modify_coverage:
    input:
        coverage="qc/mosdepth_bed/{sample}_{type}.per-base.bed.gz",
        gvcf="snv_indels/mutect2_gvcf/{sample}_{type}.merged.gvcf.gz",
    output:
        gvcf="qc/gvcf_modify_coverage/{sample}_{type}.modified.gvcf.gz",
    container:
        config.get("gvcf_modify_coverage", {}).get("container", config["default_container"])
    script:
        "../scripts/gvcf_modify_coverage.py"
