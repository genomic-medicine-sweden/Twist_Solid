rule cnv_json:
    input:
        ratios=get_cnv_ratio_file,
        segments=get_cnv_segment_file,
    output:
        json=temp("cnv_sv/cnv_html_report/{sample}_{type}.{caller}.json"),
    log:
        "cnv_sv/cnv_html_report/{sample}_{type}.{caller}.cnv_json.log",
    benchmark:
        repeat(
            "cnv_sv/cnv_html_report/{sample}_{type}.{caller}.cnv_json.benchmark.tsv",
            config.get("merge_json", {}).get("benchmark_repeats", 1),
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
        "{rule}: Create JSON representation of CNV results: {output.json}"
    script:
        "../scripts/cnv_json.py"


rule merge_json:
    input:
        annotation_bed=list(config.get("annotate_cnv", {}).values()),
        fai=config.get("reference", {}).get("fai", ""),
        germline_vcf="snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vep_annotated.filter.germline.vcf",
        json=expand(
            "cnv_sv/cnv_html_report/{{sample}}_{{type}}.{caller}.json",
            caller=config.get("svdb_merge", {}).get("cnv_callers", []),
        ),
        svdb_vcfs=expand(
            "cnv_sv/svdb_query/{{sample}}_{{type}}.svdb_query.annotate_cnv.{tag}.vcf",
            tag=list(config.get("annotate_cnv", {}).keys()),
        ),
    output:
        json=temp("cnv_sv/cnv_html_report/{sample}_{type}.merged.json"),
    params:
        skip_chromosomes=config.get("reference", {}).get("skip_chrs", None),
    log:
        "cnv_sv/cnv_html_report/{sample}_{type}.merge_json.log",
    benchmark:
        repeat(
            "cnv_sv/cnv_html_report/{sample}_{type}.merge_json.benchmark.tsv",
            config.get("merge_json", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("merge_json", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("merge_json", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("merge_json", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("merge_json", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("merge_json", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("merge_json", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("merge_json", {}).get("container", config["default_container"])
    conda:
        "../envs/merge_json.yaml"
    message:
        "{rule}: Merge JSON representations from multiple callers: {output.json}"
    script:
        "../scripts/merge_json.py"


rule cnv_html_report:
    input:
        json="cnv_sv/cnv_html_report/{sample}_{type}.merged.json",
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
