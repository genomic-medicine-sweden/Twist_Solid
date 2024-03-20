__author__ = "Patrik Smeds"
__copyright__ = "Copyright 2022, Patrik Smeds"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "GPL-3"


rule cnv_tsv_report:
    input:
        vcfs=[
            "cnv_sv/svdb_query/{sample}_{type}.{tc_method}.svdb_query.annotate_cnv.cnv_amp_genes.filter.cnv_hard_filter_amp.vcf",
            "cnv_sv/svdb_query/{sample}_{type}.{tc_method}.svdb_query.annotate_cnv.{tag}.filter.cnv_hard_filter_loh.vcf",
        ],
        org_vcfs=[
            "cnv_sv/svdb_query/{sample}_{type}.{tc_method}.svdb_query.annotate_cnv.cnv_amp_genes.vcf.gz",
            "cnv_sv/svdb_query/{sample}_{type}.{tc_method}.svdb_query.annotate_cnv.{tag}.vcf.gz",
        ],
        org_vcfs_tbi=[
            "cnv_sv/svdb_query/{sample}_{type}.{tc_method}.svdb_query.annotate_cnv.cnv_amp_genes.vcf.gz.tbi",
            "cnv_sv/svdb_query/{sample}_{type}.{tc_method}.svdb_query.annotate_cnv.{tag}.vcf.gz.tbi",
        ],
        deletions="cnv_sv/call_small_cnv_deletions/{sample}_{type}.deletions.tsv",
        amplifications="cnv_sv/call_small_cnv_amplifications/{sample}_{type}.amplifications.tsv",
        tc_file=get_tc_file,
    output:
        tsv=temp("cnv_sv/svdb_query/{sample}_{type}.{tc_method}.{tag}.cnv_report.tsv"),
        tsv_additional_only=temp("cnv_sv/svdb_query/{sample}_{type}.{tc_method}.{tag}.cnv_additional_variants_only.tsv"),
    params:
        call_small_amplifications_cn_limit=config.get("cnv_tsv_report", {}).get("amp_cn_limit", "6"),
        del_1p19q_cn_limit=config.get("cnv_tsv_report", {}).get("del_1p19q_cn_limit", "2"),
        del_1p19q_chr_arm_fraction=config.get("cnv_tsv_report", {}).get("del_1p19q_chr_arm_fraction", "0"),
        tc=get_tc,
    log:
        "cnv_sv/svdb_query/{sample}_{type}.{tc_method}.{tag}.cnv_report.tsv.log",
    benchmark:
        repeat(
            "cnv_sv/svdb_query/{sample}_{type}.{tc_method}.{tag}.cnv_report.tsv.benchmark.tsv",
            config.get("cnv_tsv_report", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnv_tsv_report", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cnv_tsv_report", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnv_tsv_report", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnv_tsv_report", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cnv_tsv_report", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnv_tsv_report", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cnv_tsv_report", {}).get("container", config["default_container"])
    message:
        "{rule}: Convert cnv vcf to a tsv file: {output.tsv}"
    script:
        "../scripts/cnv_report.py"
