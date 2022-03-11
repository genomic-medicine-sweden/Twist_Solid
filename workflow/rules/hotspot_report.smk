# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Patrik Smeds"
__copyright__ = "Copyright 2022, Patrik Smeds"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "GPL-3"


rule hotspot_report:
    input:
        hotspots=config['hotspot_report']['hotspot_mutations'],
        vcf="filtering/add_multi_snv_in_codon/{sample}_{type}.codon_snvs.sorted.vep_annotated.vcf.gz",
        vcf_index="filtering/add_multi_snv_in_codon/{sample}_{type}.codon_snvs.sorted.vep_annotated.vcf.gz.tbi",
        gvcf="qc/add_mosdepth_coverage_to_gvcf/{sample}_{type}.mosdepth.gvcf.gz",
        gvcf_index="qc/add_mosdepth_coverage_to_gvcf/{sample}_{type}.mosdepth.gvcf.gz.tbi",
    output:
        report=temp("qc/hotspot_report/{sample}_{type}.output.tsv"),
    params:
        levels=config.get("hotspot_report", {}).get("levels", []),
        sample_name=lambda wildcards: "%s_%s" % (wildcards.sample, wildcards.type),
        report_config=config.get("hotspot_report", {})["report_config"],
        chr_translation_file=config.get("hotspot_report", {})["chr_translation_file"],
        extra=config.get("hotspot_report", {}).get("extra", ""),
    log:
        "qc/hotspot_report/{sample}_{type}.output.tsv.log",
    benchmark:
        repeat(
            "qc/hotspot_report/{sample}_{type}.output.benchmark.tsv", config.get("hotspot_report", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("hotspot_report", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("hotspot_report", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("hotspot_report", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("hotspot_report", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("hotspot_report", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("hotspot_report", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("hotspot_report", {}).get("container", config["default_container"])
    conda:
        "../envs/hotspot_report.yaml"
    message:
       "{rule}: Do stuff on twist_dna_solid_uppsala/{rule}/{wildcards.sample}_{wildcards.type}.input"
    script:
        "../scripts/hotspot_report.py"
