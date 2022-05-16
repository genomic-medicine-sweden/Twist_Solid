# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Patrik Smeds"
__copyright__ = "Copyright 2021, Patrik Smeds"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "GPL3"


rule filter_cnv:
    input:
        vcf="{file}.vcf",
    output:
        vcf="{file}.include_cnv.{tag}.vcf",
    params:
        filter_config=lambda wildcards: config["filter_cnv"][wildcards.tag],
    log:
        "{file}.include_cnv.{tag}.log",
    threads: config.get("filter_cnv", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("filter_cnv", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("filter_cnv", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("filter_cnv", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("filter_cnv", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("filter_cnv", {}).get("partition", config["default_resources"]["partition"]),
    benchmark:
        repeat(
            "{file}.include_cnv.{tag}.benchmark.tsv",
            config.get("filter_cnv", {}).get("benchmark_repeats", 1),
        )
    conda:
        "../envs/filter_cnv.yaml"
    container:
        config.get("filter_cnv", {}).get("container", config["default_container"])
    message:
        "{rule}: Filter cnv based on gene bed file: {output.vcf}"
    script:
        "../scripts/filter_cnv.py"
