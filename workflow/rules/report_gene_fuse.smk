__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2022, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL3"


rule report_gene_fuse:
    input:
        fusions="fusions/gene_fuse/{sample}_{type}_gene_fuse_fusions.txt",
    output:
        report=temp("fusions/report_gene_fuse/{sample}_{type}.gene_fuse_report.tsv"),
    params:
        min_unique_reads=config.get("report_gene_fuse", {}).get("min_unique_reads", 6),
    log:
        "qc/report_gene_fuse/{sample}_{type}.gene_fuse_report.tsv.log",
    threads: config.get("report_gene_fuse", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("report_gene_fuse", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("report_gene_fuse", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("report_gene_fuse", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("report_gene_fuse", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("report_gene_fuse", {}).get("partition", config["default_resources"]["partition"]),
    benchmark:
        repeat(
            "qc/report_gene_fuse/{sample}_{type}.gene_fuse_report.tsv.tsv",
            config.get("report_gene_fuse", {}).get("benchmark_repeats", 1),
        )
    conda:
        "../envs/report_gene_fuse.yaml"
    container:
        config.get("report_gene_fuse", {}).get("container", config["default_container"])
    message:
        "{rule}: Collect and filter gene fuse dna fusions and create report: {output.fusions}"
    script:
        "../scripts/report_gene_fuse.py"
