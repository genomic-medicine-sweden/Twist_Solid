__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2022, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL3"


rule report_fusions:
    input:
        arriba = "fusions/arriba/{sample}_{type}.fusions.tsv",
        bed = config.get("report_fusions", {}).get("design_bed", ""),
        starfusion = "fusions/star_fusion/{sample}_{type}/star-fusion.fusion_predictions.abridged.tsv",
        fusioncatcher = "fusions/fusioncatcher/{sample}_{type}/final-list_candidate-fusion-genes.hg19.txt",
        bam = "fusions/star_fusion/{sample}_{type}/Aligned.out.bam",
        bai = "fusions/star_fusion/{sample}_{type}/Aligned.out.bam.bai",
    output:
        fusions = "fusions/report_fusions/{sample}_{type}.fusion_report.tsv",
        coverage = "fusions/report_fusions/{sample}_{type}.coverage.tsv",
    params:
        star_fusion=config.get("report_fusions", {}).get("filtering", []),
    log:
        "qc/report_fusions/{sample}_{type}.fusion_report.tsv.log",
    threads: config.get("report_fusions", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("report_fusions", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("report_fusions", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("report_fusions", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("report_fusions", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("report_fusions", {}).get("partition", config["default_resources"]["partition"]),
    benchmark:
        repeat(
            "qc/report_fusions/{sample}_{type}.fusion_report.tsv",
            config.get("report_fusions", {}).get("benchmark_repeats", 1),
        )
    conda:
        "../envs/report_fusions.yaml"
    container:
        config.get("report_fusions", {}).get("container", config["default_container"])
    message:
        "{rule}: Collect and filter fusions and create report: {output.fusions}"
    script:
        "../scripts/report_fusions.py"
