rule cnvkit_json:
    input:
        amp_bed=config.get("annotate_cnv", {}).get("cnv_amp_genes", []),
        cnr="cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cnr",
        cns="cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cns",
        fai=config.get("reference").get("fai"),
        loh_bed=config.get("annotate_cnv", {}).get("cnv_loh_genes", []),
        vcf="snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vep_annotated.filter.germline.vcf",
    output:
        json=temp("cnv_sv/cnvkit_report/{sample}_{type}.cnvkit.json"),
    params:
        skip_chromosomes=config.get("reference", {}).get("skip_chrs", None),
    log:
        "cnv_sv/cnvkit_report/{sample}_{type}.cnvkit.json.log",
    benchmark:
        repeat(
            "cnv_sv/cnvkit_report/{sample}_{type}.cnvkit.json.benchmark.tsv",
            config.get("cnvkit_json", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnvkit_json", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cnvkit_json", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvkit_json", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvkit_json", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cnvkit_json", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvkit_json", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cnvkit_json", {}).get("container", config["default_container"])
    conda:
        "../envs/cnvkit_json.yaml"
    message:
        "{rule}: Create a JSON representation of the cnvkit results: {output.json}"
    script:
        "../scripts/cnvkit_json.py"


rule cnvkit_html_report:
    input:
        json="cnv_sv/cnvkit_report/{sample}_{type}.cnvkit.json",
        template=config.get("cnvkit_html_report", {}).get("template", ""),
    output:
        html=temp("cnv_sv/cnvkit_report/{sample}_{type}.cnvkit.html"),
    log:
        "cnv_sv/cnvkit_report/{sample}_{type}.cnvkit.html.log",
    benchmark:
        repeat(
            "cnv_sv/cnvkit_report/{sample}_{type}.cnvkit.json.benchmark.tsv",
            config.get("cnvkit_html_report", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnvkit_html_report", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cnvkit_html_report", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvkit_html_report", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvkit_html_report", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cnvkit_html_report", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvkit_html_report", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cnvkit_html_report", {}).get("container", config["default_container"])
    conda:
        "../envs/cnvkit_html_report.yaml"
    message:
        "{rule}: Generate an interactive HTML report: {output.html}"
    script:
        "../scripts/cnvkit_html_report.py"
