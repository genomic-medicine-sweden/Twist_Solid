__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2022, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"

"bcftools mpileup -r chr19:49469088 -O v -f /data/ref_genomes/hg19/bwa/BWA_0.7.10_refseq/hg19.with.mt.fasta R21-259_Aligned.sortedByCoord.out.bam -d 20000 | bcftools call -c"

rule bcftools_id_snps:
    input:
        bam="fusions/star_fusion/{sample}_{type}.Aligned.out.sorted.bam",
        bed=config.get("bcftools_id_snps", {}).get("snps_bed", ""),
    output:
        result=temp("fusions/exon_skipping/{sample}_{type}.results.tsv"),
    log:
        "fusions/exon_skipping/{sample}_{type}.results.tsv.log",
    benchmark:
        repeat(
            "fusions/exon_skipping/{sample}_{type}.results.tsv.benchmark.tsv",
            config.get("exon_skipping", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("exon_skipping", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("exon_skipping", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("exon_skipping", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("exon_skipping", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("exon_skipping", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("exon_skipping", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("exon_skipping", {}).get("container", config["default_container"])
    conda:
        "../envs/exon_skipping.yaml"
    message:
        "{rule}: Find intergenic fusions and report them in {output.result}"
    script:
        "../scripts/exon_skipping.py"
