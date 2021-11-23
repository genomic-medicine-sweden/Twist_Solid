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
        "snv_indels/ensemble_vcf/{sample}_{type}.ensembled.vep_annotated.vcf",
    output:
        "results/dna/vcf/{sample}_{type}.ensembled.vep_annotated.vcf",
    shell:
        "cp {input} {output}"


rule copy_filtered_vcf:
    input:
        "filtering/add_multi_snv_in_codon/{sample}_{type}.codon_snvs.sorted.vcf.gz",
    output:
        "results/dna/vcf/{sample}_{type}.ensembled.vep_annotated.filtered.codon_snvs.vcf.gz",
    shell:
        "cp {input} {output}"


rule copy_fastqc:
    input:
        "qc/fastqc/{sample}_{type}_{read}_fastqc.html",
    output:
        "results/dna/qc/{sample}_{type}_{read}_fastqc.html",
    shell:
        "cp {input} {output}"


rule copy_picard_duplication_metrics:
    input:
        "qc/picard_duplication_metrics/{sample}_{type}.duplication_metrics.txt",
    output:
        "results/dna/qc/{sample}_{type}.duplication_metrics.txt",
    shell:
        "cp {input} {output}"


rule copy_picard_alignment_summary_metrics:
    input:
        "qc/picard_alignment_summary_metrics/{sample}_{type}.alignment_summary_metrics.txt",
    output:
        "results/dna/qc/{sample}_{type}.alignment_summary_metrics.txt",
    shell:
        "cp {input} {output}"


rule copy_picard_hs_metrics:
    input:
        "qc/picard_hs_metrics/{sample}_{type}.HsMetrics.txt",
    output:
        "results/dna/qc/{sample}_{type}.HsMetrics.txt",
    shell:
        "cp {input} {output}"


rule copy_picard_insert_size:
    input:
        "qc/picard_insert_size/{sample}_{type}.insert_size_metrics.txt",
    output:
        "results/dna/qc/{sample}_{type}.insert_size_metrics.txt",
    shell:
        "cp {input} {output}"


rule copy_samtools_stats:
    input:
        "qc/samtools_stats/{sample}_{type}.samtools-stats.txt",
    output:
        "results/dna/qc/{sample}_{type}.samtools-stats.txt",
    shell:
        "cp {input} {output}"


rule copy_hotspot_info:
    input:
        "qc/hotspot_info/{sample}_{type}.hotspot_info.tsv",
    output:
        "results/dna/hotspot_info/{sample}_{type}.hotspot_info.tsv",
    shell:
        "cp {input} {output}"


rule copy_msisensor_pro:
    input:
        "biomarker/msisensor_pro/{sample}_{type}",
    output:
        "results/dna/msi/{sample}_{type}.msisensor_pro.tsv",
    shell:
        "cp {input} {output}"


rule copy_tmb:
    input:
        "biomarker/tmb/{sample}_{type}.TMB.txt",
    output:
        "results/dna/tmb/{sample}_{type}.TMB.txt",
    shell:
        "cp {input} {output}"
