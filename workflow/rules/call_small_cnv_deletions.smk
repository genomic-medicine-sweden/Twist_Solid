__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2022, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule call_small_cnv_deletions:
    input:
        cnv_data="cnv_sv/gatk_denoise_read_counts/{sample}_{type}.clean.denoisedCR.tsv",
        regions_file=config.get("call_small_cnv_deletions", {}).get("regions_file", ""),
    output:
        deletions=temp("cnv_sv/call_small_cnv_deletions/{sample}_{type}.deletions.tsv"),
    params:
        window_size=config.get("call_small_cnv_deletions", {}).get("window_size", 4),
        region_max_size=config.get("call_small_cnv_deletions", {}).get("region_max_size", 30),
        min_nr_stdev_diff=config.get("call_small_cnv_deletions", {}).get("min_nr_stdev_diff", 2.5),
        min_log_odds_diff=config.get("call_small_cnv_deletions", {}).get("min_log_odds_diff", 0.3),
    log:
        "cnv_sv/call_small_cnv_deletions/{sample}_{type}.deletions.tsv.log",
    benchmark:
        repeat(
            "cnv_sv/call_small_cnv_deletions/{sample}_{type}.deletions.tsv.benchmark.tsv",
            config.get("call_small_cnv_deletions", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("call_small_cnv_deletions", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("call_small_cnv_deletions", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("call_small_cnv_deletions", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("call_small_cnv_deletions", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("call_small_cnv_deletions", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("call_small_cnv_deletions", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("call_small_cnv_deletions", {}).get("container", config["default_container"])
    conda:
        "../envs/call_small_cnv_deletions.yaml"
    message:
        "{rule}: call small deletions in cnv data into {output.deletions}"
    script:
        "../scripts/call_small_cnv_deletions.py"
