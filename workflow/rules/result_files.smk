# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas A"
__copyright__ = "Copyright 2021, Jonas A"
__email__ = "jonas.almlof@igp.uu.se"
__license__ = "GPL-3"

rule copy_bam:
    input:
        "alignment/merge_bam/{sample}_{type}.bam",
    output:
        "results/dna/bam/{sample}_{type}.bam",
    shell:
        "cp {input} {output}"
