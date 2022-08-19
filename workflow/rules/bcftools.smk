__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2022, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule bcftools_id_snps:
    input:
        bam="fusions/star_fusion/{sample}_{type}/Aligned.out.sorted.bam",
        bai="fusions/star_fusion/{sample}_{type}/Aligned.out.sorted.bam.bai",
        bed=config.get("bcftools_id_snps", {}).get("snps_bed", "missing bcftools_id_snps snps_bed file"),
        ref=config.get("reference", {}).get("fasta_rna", "missing reference fasta_rna"),
    output:
        vcf=temp("snv_indels/bcftools_id_snps/{sample}_{type}.id_snps.vcf"),
    log:
        "snv_indels/bcftools_id_snps/{sample}_{type}.id_snps.vcf.log",
    benchmark:
        repeat(
            "snv_indels/bcftools_id_snps/{sample}_{type}.id_snps.vcf.benchmark.tsv",
            config.get("bcftools_id_snps", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("bcftools_id_snps", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("bcftools_id_snps", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcftools_id_snps", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bcftools_id_snps", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bcftools_id_snps", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("bcftools_id_snps", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("bcftools_id_snps", {}).get("container", config["default_container"])
    conda:
        "../envs/bcftools.yaml"
    message:
        "{rule}: call id SNPs in the RNA data {output.vcf}"
    shell:
        "(bcftools mpileup "
        "-R {input.bed} "
        "-O u "
        "-f {input.ref} "
        "-d 1000000 "
        "| bcftools call "
        "--skip-variants indels "
        "-m -O v "
        "-o {output.vcf}) "
        "&> {log}"
