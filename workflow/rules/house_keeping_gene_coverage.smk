__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2022, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule house_keeping_gene_coverage:
    input:
        bam="fusions/star_fusion/{sample}_{type}/Aligned.out.sorted.bam",
        bai="fusions/star_fusion/{sample}_{type}/Aligned.out.sorted.bam.bai",
        bed=config.get("reference", {}).get("design_bed_rna", ""),
    output:
        result=temp("qc/house_keeping_gene_coverage/{sample}_{type}.house_keeping_gene_coverage.tsv"),
    log:
        "fusions/house_keeping_gene_coverage/{sample}_{type}.house_keeping_gene_coverage.tsv.log",
    benchmark:
        repeat(
            "fusions/house_keeping_gene_coverage/{sample}_{type}.house_keeping_gene_coverage.tsv.benchmark.tsv",
            config.get("house_keeping_gene_coverage", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("house_keeping_gene_coverage", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("house_keeping_gene_coverage", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("house_keeping_gene_coverage", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("house_keeping_gene_coverage", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("house_keeping_gene_coverage", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("house_keeping_gene_coverage", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("house_keeping_gene_coverage", {}).get("container", config["default_container"])
    conda:
        "../envs/house_keeping_gene_coverage.yaml"
    message:
        "{rule}: Find intergenic fusions and report them in {output.result}"
    script:
        "../scripts/house_keeping_gene_coverage.py"
