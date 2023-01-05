__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2022, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule fix_vcf_ad_for_qci:
    input:
        vcf="annotation/add_multi_snv_in_codon/{sample}_{type}.background_annotation.include.exon.filter.snv_hard_filter.codon_snvs.sorted.vep_annotated.vcf",
    output:
        vcf="annotation/fix_vcf_ad_for_qci/{sample}_{type}.background_annotation.include.exon.filter.snv_hard_filter.codon_snvs.sorted.vep_annotated.qci.vcf",
    log:
        "annotation/fix_vcf_ad_for_qci/{sample}_{type}.background_annotation.include.exon.filter.snv_hard_filter.codon_snvs.sorted.vep_annotated.qci.vcf.log",
    benchmark:
        repeat(
            "annotation/fix_vcf_ad_for_qci/{sample}_{type}.background_annotation.include.exon.filter.snv_hard_filter.codon_snvs.sorted.vep_annotated.qci.vcf.benchmark.tsv",
            config.get("fix_vcf_ad_for_qci", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("fix_vcf_ad_for_qci", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("fix_vcf_ad_for_qci", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("fix_vcf_ad_for_qci", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("fix_vcf_ad_for_qci", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("fix_vcf_ad_for_qci", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("fix_vcf_ad_for_qci", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("fix_vcf_ad_for_qci", {}).get("container", config["default_container"])
    conda:
        "../envs/fix_vcf_ad_for_qci.yaml"
    message:
        "{rule}: Correct AD so the correct AF is shown in QCI for vcf {input.vcf}"
    script:
        "../scripts/fix_vcf_ad_for_qci.py"
