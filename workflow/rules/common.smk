# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas A"
__copyright__ = "Copyright 2021, Jonas A"
__email__ = "jonas.almlof@igp.uu.se"
__license__ = "GPL-3"

import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *

min_version("6.10.0")

### Set and validate config file


configfile: "config/config.yaml"


validate(config, schema="../schemas/config.schema.yaml")
config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")


### Read and validate samples file

samples = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

### Read and validate units file

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "type"], drop=False)
validate(units, schema="../schemas/units.schema.yaml")

### Set wildcard constraints


wildcard_constraints:
    sample="|".join(samples.index),
    unit="N|T|R",


def compile_result_file_list():
    files = {
        "results/dna/bam": [
            {"in": ["alignment/merge_bam", ".bam"], "out": ".bam"},
            {"in": ["alignment/merge_bam", ".bam.bai"], "out": ".bam.bai"},
        ],
        "results/dna/vcf": [
            {"in": ["snv_indels/ensemble_vcf", ".ensembled.vcf.gz"], "out": ".ensembled.vcf.gz"},
            {"in": ["snv_indels/ensemble_vcf", ".ensembled.vep_annotated.vcf"], "out": ".ensembled.vep_annotated.vcf"},
            {
                "in": ["filtering/add_multi_snv_in_codon", ".codon_snvs.sorted.vcf.gz"],
                "out": ".ensembled.vep_annotated.filtered.codon_snvs.vcf.gz"
            },
            {
                "in": ["filtering/add_multi_snv_in_codon", ".codon_snvs.sorted.include.nocnv.vcf.gz"],
                "out": ".ensembled.vep_annotated.filtered.codon_snvs.nocnv.vcf.gz"
            },
            {
                "in": ["filtering/add_multi_snv_in_codon", ".codon_snvs.sorted.include.exon.vcf.gz"],
                "out": ".ensembled.vep_annotated.filtered.codon_snvs.exon_only.vcf.gz"
            },
        ],
        "results/dna/gvcf": [
            {"in": ["snv_indels/mutect2_gvcf", ".gvcf.gz"], "out": ".gvcf.gz"},
        ],
        "results/dna/qc": [
            {"in": ["qc/picard_collect_duplication_metrics", ".duplication_metrics.txt"], "out": ".duplication_metrics.txt"},
            {
                "in": ["qc/picard_collect_duplication_metrics", ".alignment_summary_metrics.txt"],
                "out": ".alignment_summary_metrics.txt"
            },
            {"in": ["qc/picard_collect_hs_metrics", ".HsMetrics.txt"], "out": ".HsMetrics.txt"},
            {"in": ["qc/picard_collect_insert_size_metrics", ".insert_size_metrics.txt"], "out": ".insert_size_metrics.txt"},
            {"in": ["qc/samtools_stats", ".samtools-stats.txt"], "out": ".samtools-stats.txt"},
        ],
        "results/dna/msi": [
            {"in": ["biomarker/msisensor_pro", ".msisensor_pro.tsv"], "out": ".msisensor_pro.tsv"},
        ],
        "results/dna/tmb": [
            {"in": ["biomarker/tmb", ".TMB.txt"], "out": ".TMB.txt"},
        ],
        "results/dna/hrd": [
            {"in": ["biomarker/hrd", ".hrd_score.txt"], "out": ".hrd_score.txt"},
        ],
        "results/dna/fusions": [
            {"in": ["fusions/gene_fuse", "_gene_fuse_fusions.txt"], "out": "_gene_fuse_fusions.txt"},
        ],
        "results/dna/cnv": [
            {"in": ["cnv_sv/cnvkit_call", ".cnvkit_loh.cns"], "out": ".cnvkit_loh.cns"},
            {"in": ["cnv_sv/gatk_cnv_call_copy_ratio_segments", ".gatk_cnv.seg"], "out": ".gatk_cnv.seg"},
            {"in": ["cnv_sv/gatk_cnv_vcf", ".gatk_cnv.vcf"], "out": ".gatk_cnv.vcf"},
            {"in": ["cnv_sv/cnvkit_vcf", ".cnvkit.vcf"], "out": ".cnvkit.vcf"},
            {"in": ["cnv_sv/svdb_merge", ".merged.vcf"], "out": ".merged.vcf"},
        ],
    }
    output_files = [
        "%s/%s_%s%s" % (prefix, sample, unit_type, suffix)
        for prefix in files.keys()
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
        for suffix in [f_dict["out"] for p_out in files.keys() for f_dict in files[p_out]]
    ]
    input_files = [
        "%s/%s_%s%s" % (prefix, sample, unit_type, suffix)
        for prefix in [f_dict["in"][0] for p_out in files.keys() for f_dict in files[p_out]]
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
        for suffix in [f_dict["in"][1] for p_out in files.keys() for f_dict in files[p_out]]
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
    print(input_files)
    print(output_files)
    return input_files, output_files


def compile_output_list(wildcards):
    return output_files


input_files, output_files = compile_result_file_list()
