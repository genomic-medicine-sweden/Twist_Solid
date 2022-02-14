# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas A"
__copyright__ = "Copyright 2021, Jonas A"
__email__ = "jonas.almlof@igp.uu.se"
__license__ = "GPL-3"


def compile_result_file_list(wildcards):
    files = {
        "results/dna/bam": {
            "output": [
                ".bam",
                ".bam.bai",
            ],
            "input": [
                ["alignment/merge_bam", ".bam"],
                ["alignment/merge_bam", ".bam.bai"],
            ],"
        },
        "results/dna/vcf": {
            "output": [
                ".ensembled.vcf.gz",
                ".ensembled.vep_annotated.vcf",
                ".ensembled.vep_annotated.filtered.codon_snvs.vcf.gz",
                ".ensembled.vep_annotated.filtered.codon_snvs.nocnv.vcf.gz",
                ".ensembled.vep_annotated.filtered.codon_snvs.exon_only.vcf.gz",
            ],
            "input": [
                ["snv_indels/ensemble_vcf", ".ensembled.vcf.gz"],
                ["snv_indels/ensemble_vcf", ".ensembled.vep_annotated.vcf"],
                ["filtering/add_multi_snv_in_codon", ".codon_snvs.sorted.vcf.gz"],
                ["filtering/add_multi_snv_in_codon", ".codon_snvs.sorted.include.nocnv.vcf.gz"],
                ["filtering/add_multi_snv_in_codon", ".codon_snvs.sorted.include.exon.vcf.gz"],
            ],
        },
        "results/dna/gvcf": {
            "output": [".gvcf.gz"],
            "input": [["snv_indels/mutect2_gvcf", ".gvcf.gz"]],
        },
        "results/dna/qc":
            "output": [
                ".duplication_metrics.txt",
                ".alignment_summary_metrics.txt",
                ".HsMetrics.txt",
                ".insert_size_metrics.txt",
                ".samtools-stats.txt",
            ],
            "input": [
                ["qc/picard_collect_duplication_metrics", ".duplication_metrics.txt"],
                ["qc/picard_collect_duplication_metrics", ".alignment_summary_metrics.txt"],
                ["qc/picard_collect_hs_metrics", ".HsMetrics.txt"],
                ["qc/picard_collect_insert_size_metrics", ".insert_size_metrics.txt"],
                ["qc/samtools_stats", ".samtools-stats.txt"],
            ],
        },
        "results/dna/msi": {
            "output": [".msisensor_pro.tsv"],
            "input": [["biomarker/msisensor_pro", ".msisensor_pro.tsv"]],
        },
        "results/dna/tmb": {
            "output": [".TMB.txt"],
            "input": [["biomarker/tmb", ".TMB.txt"]],
        },
        "results/dna/hrd": {
            "output": [".hrd_score.txt"],
            "input": [["biomarker/hrd", ".hrd_score.txt"]],
        },
        "results/dna/fusions": {
            "output": ["_gene_fuse_fusions.txt"],
            "input": [["fusions/gene_fuse", "_gene_fuse_fusions.txt"]],
        },
        "results/dna/cnv": {
            "output": [
                ".cnvkit_loh.cns",
                ".gatk_cnv.seg",
                ".gatk_cnv.vcf",
                ".cnvkit.vcf",
                ".merged.vcf",
            ],
            "input": [
                ["cnv_sv/cnvkit_call", ".cnvkit_loh.cns"],
                ["cnv_sv/gatk_cnv_call_copy_ratio_segments", ".gatk_cnv.seg"],
                ["cnv_sv/gatk_cnv_vcf", ".gatk_cnv.vcf"],
                ["cnv_sv/cnvkit_vcf", ".cnvkit.vcf"],
                ["cnv_sv/svdb_merge", ".merged.vcf"],
            ],
        },
    }
    output_files = [
        "%s/%s_%s%s" % (prefix, sample, unit_type, suffix)
        for prefix in files.keys()
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
        for suffix in files[prefix][output]
    ]
    input_files = [
        "%s/%s_%s%s" % (prefix, sample, unit_type, suffix)
        for files[prefix][input][0] in files.keys()
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
        for suffix in files[prefix][input][0]
    ]
    output_files.append(
        [
            "results/dna/vcf/%s_%s_%s.vcf.gz" % (caller, sample, t)
            for caller in ["mutect2", "vardict"]
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
        ]
    )
    input_files.append(
        [
            "snv_indels/%s/%s_%s.merged.vcf.gz" % (caller, sample, t)
            for caller in ["mutect2", "vardict"]
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
        ]
    )
    output_files.append(
        [
            "results/dna/qc/%s_%s_%s_fastqc.html" % (sample, t, read)
            for read in ["fastq1", "fastq2"]
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
        ]
    )
    input_files.append(
        [
            "qc/fastqc/%s_%s_%s_fastqc.html" % (sample, t, read)
            for read in ["fastq1", "fastq2"]
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
        ]
    )
    output_files.append("results/dna/qc/MultiQC.html")
    input_files.append("qc/multiqc/MultiQC.html")
    return input_files, output_files


input_files, output_files = compile_result_file_list()


rule copy_results_files:
    input:
        input_files,
    output:
        output_files,
    run:
        import subprocess
        i = 0
        for file in input[0] :
            subprocess.run(["cp", file, output[0][i]]
            i += 1
