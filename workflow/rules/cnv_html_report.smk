rule cnv_json:
    input:
        amp_bed=config.get("annotate_cnv", {}).get("cnv_amp_genes", []),
        cnr="cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cnr",
        cns="cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cns",
        gatk_ratios="cnv_sv/gatk_denoise_read_counts/{sample}_{type}.clean.denoisedCR.tsv",
        gatk_segments="cnv_sv/gatk_model_segments/{sample}_{type}.clean.cr.seg",
        germline_vcf="snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vep_annotated.filter.germline.vcf",
        fai=config.get("reference").get("fai"),
        loh_bed=config.get("annotate_cnv", {}).get("cnv_loh_genes", []),
        svdb_vcfs=[
            "cnv_sv/svdb_query/{sample}_{type}.svdb_query.annotate_cnv.cnv_amp_genes.filter.cnv_hard_filter_amp.vcf.gz",
            "cnv_sv/svdb_query/{sample}_{type}.svdb_query.annotate_cnv.cnv_loh_genes.filter.cnv_hard_filter_loh.vcf.gz",
        ],
        svdb_tbis=[
            "cnv_sv/svdb_query/{sample}_{type}.svdb_query.annotate_cnv.cnv_amp_genes.filter.cnv_hard_filter_amp.vcf.gz.tbi",
            "cnv_sv/svdb_query/{sample}_{type}.svdb_query.annotate_cnv.cnv_loh_genes.filter.cnv_hard_filter_loh.vcf.gz.tbi",
        ],
    output:
        json=temp("cnv_sv/cnv_html_report/{sample}_{type}.cnv.json"),
    params:
        skip_chromosomes=config.get("reference", {}).get("skip_chrs", None),
    log:
        "cnv_sv/cnv_html_report/{sample}_{type}.cnv.json.log",
    benchmark:
        repeat(
            "cnv_sv/cnv_html_report/{sample}_{type}.cnv.json.benchmark.tsv",
            config.get("cnv_json", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnv_json", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cnv_json", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnv_json", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnv_json", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cnv_json", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnv_json", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cnv_json", {}).get("container", config["default_container"])
    conda:
        "../envs/cnv_json.yaml"
    message:
        "{rule}: Create a JSON representation of the cnv results: {output.json}"
    script:
        "../scripts/cnv_json.py"


rule cnv_html_report:
    input:
        json="cnv_sv/cnv_html_report/{sample}_{type}.cnv.json",
        template=config.get("cnv_html_report", {}).get("template", ""),
    output:
        html=temp("cnv_sv/cnv_html_report/{sample}_{type}.cnv.html"),
    log:
        "cnv_sv/cnv_html_report/{sample}_{type}.cnv.html.log",
    benchmark:
        repeat(
            "cnv_sv/cnv_html_report/{sample}_{type}.cnv.json.benchmark.tsv",
            config.get("cnv_html_report", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnv_html_report", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cnv_html_report", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnv_html_report", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnv_html_report", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cnv_html_report", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnv_html_report", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cnv_html_report", {}).get("container", config["default_container"])
    conda:
        "../envs/cnv_html_report.yaml"
    message:
        "{rule}: Generate an interactive HTML report: {output.html}"
    script:
        "../scripts/cnv_html_report.py"
