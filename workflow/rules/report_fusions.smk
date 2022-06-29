__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2022, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL3"


rule report_fusions:
    input:
        arriba="fusions/arriba/{sample}_{type}.fusions.tsv",
        bam="fusions/star_fusion/{sample}_{type}/Aligned.out.sorted.bam",
        bai="fusions/star_fusion/{sample}_{type}/Aligned.out.sorted.bam.bai",
        bed=config.get("reference", {}).get("design_bed_rna", ""),
        bed_extra_annotation=config.get("report_fusions", {}).get("annotation_bed", ""),
        fusioncatcher="fusions/fusioncatcher/{sample}_{type}/final-list_candidate-fusion-genes.hg19.txt",
        star_fusion="fusions/star_fusion/{sample}_{type}/star-fusion.fusion_predictions.abridged.coding_effect.tsv",
    output:
        fusions=temp("fusions/report_fusions/{sample}_{type}.fusion_report.tsv"),
    params:
        fusioncather_flag_low_support=config.get("report_fusions", {}).get("fusioncather", {}).get("flag_low_support", 15),
        fusioncather_low_support=config.get("report_fusions", {}).get("fusioncather", {}).get("low_support", 3),
        fusioncather_low_support_fp_genes=config.get("report_fusions", {}).get("fusioncather", {}).get("low_support_fp_genes", 20),
        fusioncather_low_support_inframe=config.get("report_fusions", {}).get("fusioncather", {}).get("low_support_inframe", 6),
        star_fusion_flag_low_support=config.get("report_fusions", {}).get("star_fusion", {}).get("flag_low_support", 15),
        star_fusion_low_support=config.get("report_fusions", {}).get("star_fusion", {}).get("low_support", 2),
        star_fusion_low_support_fp_genes=config.get("report_fusions", {}).get("star_fusion", {}).get("low_support_fp_genes", 20),
        star_fusion_low_support_inframe=config.get("report_fusions", {}).get("star_fusion", {}).get("low_support_inframe", 6),
    log:
        "qc/report_fusions/{sample}_{type}.fusion_report.tsv.log",
    threads: config.get("report_fusions", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("report_fusions", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("report_fusions", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("report_fusions", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("report_fusions", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("report_fusions", {}).get("time", config["default_resources"]["time"]),
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
