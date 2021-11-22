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


rule copy_caller_vcf:
    input:
        "snv_indels/{caller}/{sample}_{type}.merged.vcf.gz",
    output:
        "results/dna/vcf/{caller}_{sample}_{type}.vcf.gz",
    shell:
        "cp {input} {output}"


rule copy_ensembled_vcf:
    input:
        "snv_indels/ensemble_vcf/{sample}_{type}.ensembled.vcf.gz",
    output:
        "results/dna/vcf/{sample}_{type}.ensembled.vcf.gz",
    shell:
        "cp {input} {output}"


rule copy_gvcf:
    input:
        "snv_indels/mutect2_gvcf/{sample}_{type}.merged.gvcf.gz",
    output:
        "results/dna/gvcf/{sample}_{type}.gvcf.gz",
    shell:
        "cp {input} {output}"


rule copy_vep_vcf:
    input:
        "annotation/ensemble_vcf/{sample}_{type}.ensembled.vep_annotated.vcf",
    output:
        "results/dna/vcf/{sample}_{type}.ensembled.vep_annotated.vcf",
    shell:
        "cp {input} {output}"
