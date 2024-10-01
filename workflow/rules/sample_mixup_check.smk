__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2024, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule sample_mixup_check:
    input:
        id_snp_vcf_dna=[
            "snv_indels/bcftools_id_snps/%s_%s.id_snps.vcf" % (sample, unit_type)
            for sample in samples.index
            for unit_type in get_unit_types(units, sample)
            if unit_type == "T"
        ],
        id_snp_vcf_rna=[
            "snv_indels/bcftools_id_snps/%s_%s.id_snps.vcf" % (sample, unit_type)
            for sample in samples.index
            for unit_type in get_unit_types(units, sample)
            if unit_type == "R"
        ],
    output:
        mixup_report="qc/sample_mixup_check/sample_mixup_check.tsv",
    params:
        extra=config.get("sample_mixup_check", {}).get("extra", ""),
        match_cutoff=config.get("sample_mixup_check", {}).get("match_cutoff", "0.7"),
    log:
        "twist_solid/sample_mixup_check/sample_mixup_check.tsv.log",
    benchmark:
        repeat(
            "twist_solid/sample_mixup_check/sample_mixup_check.tsv.benchmark.tsv",
            config.get("sample_mixup_check", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("sample_mixup_check", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("sample_mixup_check", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("sample_mixup_check", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("sample_mixup_check", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("sample_mixup_check", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("sample_mixup_check", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("sample_mixup_check", {}).get("container", config["default_container"])
    message:
        "{rule}: find most similar sample in rna ID-snps compared to dna ID-snps"
    script:
        "../scripts/sample_mixup_check.py"
