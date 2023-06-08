__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2023, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule purecn_modify_vcf:
    input:
        vcf="snv_indels/gatk_mutect2/{sample}_{type}.normalized.sorted.vep_annotated.filter.snv_hard_filter_purecn.bcftools_annotated_purecn.vcf",
    output:
        vcf="cnv_sv/purecn_modify_vcf/{sample}_{type}.normalized.sorted.vep_annotated.filter.snv_hard_filter_purecn.bcftools_annotated_purecn.mbq.vcf",
    params:
        extra=config.get("purecn_modify_vcf", {}).get("extra", ""),
    log:
        "cnv_sv/purecn_modify_vcf/{sample}_{type}.normalized.sorted.vep_annotated.filter.snv_hard_filter_purecn.bcftools_annotated_purecn.mbq.vcf.log",
    benchmark:
        repeat(
            "cnv_sv/purecn_modify_vcf/{sample}_{type}.normalized.sorted.vep_annotated.filter.snv_hard_filter_purecn.bcftools_annotated_purecn.mbq.vcf.benchmark.tsv",
            config.get("purecn_modify_vcf", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("purecn_modify_vcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("purecn_modify_vcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("purecn_modify_vcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("purecn_modify_vcf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("purecn_modify_vcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("purecn_modify_vcf", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("purecn_modify_vcf", {}).get("container", config["default_container"])
    message:
        "{rule}: modify mbq in {input.vcf}"
    shell:
        # new lines are used to tell the linter that it is not an absolut path
        "(sed -r 's/(.*)(\\MBQ=[0-9]+,)([0-9]+)(.*)/echo \"\\1\\2$((\\3+5))\\4\"/"
        "ge' {input.vcf} | "
        "sed -r 's/(.*)(\\MBQ=)([0-9]+)(.*)/echo \"\\1\\2$((\\3+5))\\4\"/"
        "ge' > "
        "{output.vcf}) &> {log}"
