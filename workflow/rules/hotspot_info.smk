# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Patrik Smeds"
__copyright__ = "Copyright 2021, Patrik Smeds"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "GPL3"


def get_flowcell(units, wildcards):
    flowcells = set([u.flowcell for u in get_units(units, wildcards)])
    if len(flowcells) > 1:
        raise ValueError("Sample type combination from different sequence flowcells")
    return flowcells.pop()


rule hotspot_info:
    input:
        bam="alignment/merge_bam/{sample}_{type}.bam",
        bai="alignment/merge_bam/{sample}_{type}.bam.bai",
        vcf="snv_indels/ensemble_vcf/{sample}_{type}.ensembled.vep_annotated.vcf",
        hotspots=config.get("hotspot_info", {}).get("hotspot_mutations", ""),
        background_panel=config.get("reference", {}).get("background", ""),
        background_run=lambda wildcards: "annotation/calculate_seqrun_background/%s_seqrun_background.tsv"
        % get_flowcell(units, wildcards),
        gvcf="qc/add_mosdepth_coverage_to_gvcf/{sample}_{type}.mosdepth.gvcf.gz",
        gvcf_index="qc/add_mosdepth_coverage_to_gvcf/{sample}_{type}.mosdepth.gvcf.gz.tbi",
    output:
        low_coverage=temp("qc/hotspot_info/{sample}_{type}.hotspot_low_coverage.txt"),
        hotspot_info=temp("qc/hotspot_info/{sample}_{type}.hotspot_coverage_info.tsv"),
    params:
        min_coverage=config.get("hotspot_info", {}).get("min_coverage", 200),
    log:
        "qc/hotspot_info/{sample}_{type}.log",
    threads: config.get("hotspot_info", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("hotspot_info", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("hotspot_info", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("hotspot_info", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("hotspot_info", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("hotspot_info", {}).get("partition", config["default_resources"]["partition"]),
    benchmark:
        repeat(
            "qc/hotspot_info/{sample}_{type}.fastq.gz.fastp_trimming.benchmark.tsv",
            config.get("hotspot_info", {}).get("benchmark_repeats", 1),
        )
    conda:
        "../envs/hotspot_info.yaml"
    container:
        config.get("hotspot_info", {}).get("container", config["default_container"])
    message:
        "{rule}: Annotate hotspot file and report low coverage: qc/hotspot_info/{wildcards.sample}_{wildcards.type}"
    script:
        "../scripts/hotspot_info.py"
