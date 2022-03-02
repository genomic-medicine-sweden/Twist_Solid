# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2022, Jonas Almlöf"
__email__ = "jonas.almlof@igp.uu.se"
__license__ = "GPL-3"


rule add_mosdepth_coverage_to_gvcf:
    input:
        coverage="qc/mosdepth_bed/{sample}_{type}.per-base.bed.gz",
        gvcf="snv_indels/mutect2_gvcf/{sample}_{type}.merged.gvcf.gz",
    output:
        gvcf="qc/add_mosdepth_coverage_to_gvcf/{sample}_{type}.mosdepth.gvcf.gz",
    container:
        config.get("add_mosdepth_coverage_to_gvcf", {}).get("container", config["default_container"])
    script:
        "../scripts/add_mosdepth_coverage_to_gvcf.py"
