# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Patrik Smeds"
__copyright__ = "Copyright 2022, Patrik Smeds"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "GPL-3"


rule hotspot_report:
    input:
        hotspots=lambda wildcards: config["hotspot_report"]["hotspot_mutations"][wildcards.tag],
        vcf=lambda wildcards: "%s%s" % (get_hotspot_report_vcf_input(wildcards), ".vep_annotated.vcf.gz"),
        vcf_index=lambda wildcards: "%s%s" % (get_hotspot_report_vcf_input(wildcards), ".vep_annotated.vcf.gz.tbi"),
        vcf_file_wo_pick=lambda wildcards: "%s%s" % (get_hotspot_report_vcf_input(wildcards), ".vep_annotated_wo_pick.vcf.gz"),
        vcf_file_wo_pick_index=lambda wildcards: "%s%s"
        % (get_hotspot_report_vcf_input(wildcards), ".vep_annotated_wo_pick.vcf.gz"),
        gvcf="qc/add_mosdepth_coverage_to_gvcf/{sample}_{type}.mosdepth.g.vcf.gz",
        gvcf_index="qc/add_mosdepth_coverage_to_gvcf/{sample}_{type}.mosdepth.g.vcf.gz.tbi",
    output:
        report=temp("qc/hotspot_report/{sample}_{type}.{tag}.output.tsv"),
    params:
        levels=config.get("hotspot_report", {}).get("levels", []),
        sample_name=lambda wildcards: "%s_%s" % (wildcards.sample, wildcards.type),
        report_config=config.get("hotspot_report", {})["report_config"],
        chr_translation_file=config.get("hotspot_report", {})["chr_translation_file"],
        extra=config.get("hotspot_report", {}).get("extra", ""),
    log:
        "qc/hotspot_report/{sample}_{type}.{tag}.output.tsv.log",
    benchmark:
        repeat(
            "qc/hotspot_report/{sample}_{type}.{tag}.output.benchmark.tsv",
            config.get("hotspot_report", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("hotspot_report", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("hotspot_report", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("hotspot_report", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("hotspot_report", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("hotspot_report", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("hotspot_report", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("hotspot_report", {}).get("container", config["default_container"])
    message:
        "{rule}: Make a coverage and mutations report: {output.report}"
    script:
        "../scripts/hotspot_report.py"


rule hotspot_report_aa_translate:
    input:
        report="qc/hotspot_report/{sample}_{type}.{tag}.output.tsv",
    output:
        report=temp("qc/hotspot_report/{sample}_{type}.{tag}.aa_translated.tsv"),
    log:
        "twist_solid/hotspot_report_aa_translate/{sample}_{type}.{tag}.aa_translated.tsv.log",
    benchmark:
        repeat(
            "twist_solid/hotspot_report_aa_translate/{sample}_{type}.{tag}.aa_translated.tsv.benchmark.tsv",
            config.get("hotspot_report_aa_translate", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("hotspot_report_aa_translate", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("hotspot_report_aa_translate", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("hotspot_report_aa_translate", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("hotspot_report_aa_translate", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("hotspot_report_aa_translate", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("hotspot_report_aa_translate", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("hotspot_report_aa_translate", {}).get("container", config["default_container"])
    message:
        "{rule}: Translate aa three letter codes into one letter code in {input.report}"
    script:
        "../scripts/aa_translate.py"
