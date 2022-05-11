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
from hydra_genetics import min_version as hydra_min_version

hydra_min_version("0.10.0")

min_version("6.10.0")

### Set and validate config file


configfile: "config/config.yaml"


validate(config, schema="../schemas/config.schema.yaml")
config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")

if workflow.use_singularity is True:
    validate(config, schema="../schemas/singularity.schema.yaml")

### Read and validate samples file

samples = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

### Read and validate units file

units = (
    pandas.read_table(config["units"], dtype=str)
    .set_index(["sample", "type", "flowcell", "lane", "barcode"], drop=False)
    .sort_index()
)
validate(units, schema="../schemas/units.schema.yaml")

### Set wildcard constraints


wildcard_constraints:
    barcode="[A-Z+]+",
    chr="[^_]+",
    flowcell="[A-Z0-9]+",
    lane="L[0-9]+",
    sample="|".join(get_samples(samples)),
    type="N|T|R",


def compile_result_file_list():
    files = [
        #{"in": ["alignment/samtools_merge_bam", ".bam"], "out": ["results/dna/bam", ".bam"]},
        #{"in": ["alignment/samtools_merge_bam", ".bam.bai"], "out": ["results/dna/bam", ".bam.bai"]},
        {
            "in": ["snv_indels/bcbio_variation_recall_ensemble", ".ensembled.vcf.gz"],
            "out": ["results/dna/vcf", ".ensembled.vcf.gz"],
        },
        {
            "in": ["annotation/background_annotation", ".background_annotation.vcf.gz"],
            "out": ["results/dna/vcf", ".annotated.vcf.gz"],
        },
        {
            "in": ["annotation/background_annotation", ".background_annotation.include.nocnv.vcf.gz"],
            "out": ["results/dna/vcf", ".annotated.nocnv.vcf.gz"],
        },
        {
            "in": ["annotation/background_annotation", ".background_annotation.include.exon.soft_filter.vcf"],
            "out": ["results/dna/vcf", ".annotated.exon_only.soft_filter.vcf"],
        },
        {
            "in": ["annotation/background_annotation", ".background_annotation.include.exon.hard_filter.vcf"],
            "out": ["results/dna/vcf", ".annotated.exon_only.hard_filter.vcf"],
        },
        {
            "in": ["qc/picard_collect_duplication_metrics", ".duplication_metrics.txt"],
            "out": ["results/dna/qc", ".duplication_metrics.txt"],
        },
        {
            "in": ["qc/picard_collect_alignment_summary_metrics", ".alignment_summary_metrics.txt"],
            "out": ["results/dna/qc", ".alignment_summary_metrics.txt"],
        },
        {"in": ["qc/picard_collect_hs_metrics", ".HsMetrics.txt"], "out": ["results/dna/qc", ".HsMetrics.txt"]},
        {
            "in": ["qc/picard_collect_insert_size_metrics", ".insert_size_metrics.txt"],
            "out": ["results/dna/qc", ".insert_size_metrics.txt"],
        },
        {"in": ["qc/samtools_stats", ".samtools-stats.txt"], "out": ["results/dna/qc", ".samtools-stats.txt"]},
        {"in": ["qc/add_mosdepth_coverage_to_gvcf", ".mosdepth.g.vcf.gz"], "out": ["results/dna/gvcf", ".mosdepth.g.vcf.gz"]},
        {"in": ["qc/hotspot_report", ".output.tsv"], "out": ["results/dna/qc", ".hotspot.tsv"]},
        {"in": ["qc/hotspot_info", ".hotspot_coverage_info.tsv"], "out": ["results/dna/qc", ".hotspot_coverage_info.tsv"]},
        {"in": ["biomarker/msisensor_pro", ""], "out": ["results/dna/msi", ".msisensor_pro.score.tsv"]},
        {"in": ["biomarker/tmb", ".TMB.txt"], "out": ["results/dna/tmb", ".TMB.txt"]},
        {"in": ["biomarker/hrd", ".hrd_score.txt"], "out": ["results/dna/hrd", ".hrd_score.txt"]},
        {"in": ["fusions/gene_fuse", "_gene_fuse_fusions.txt"], "out": ["results/dna/fusions", ".gene_fuse_fusions.txt"]},
        {"in": ["cnv_sv/cnvkit_call", ".loh.cns"], "out": ["results/dna/cnv", ".cnvkit.loh.cns"]},
        {
            "in": ["cnv_sv/gatk_cnv_call_copy_ratio_segments", ".clean.calledCNVs.seg"],
            "out": ["results/dna/cnv", ".gatk_cnv.seg"],
        },
        {"in": ["cnv_sv/gatk_cnv_vcf", ".vcf"], "out": ["results/dna/cnv", ".gatk_cnv.vcf"]},
        {"in": ["cnv_sv/cnvkit_vcf", ".vcf"], "out": ["results/dna/cnv", ".cnvkit.vcf"]},
        {"in": ["cnv_sv/cnvkit_scatter", ".png"], "out": ["results/dna/cnv", ".cnvkit.scatter.png"]},
        {"in": ["cnv_sv/cnvkit_diagram", ".pdf"], "out": ["results/dna/cnv", ".cnvkit.diagram.pdf"]},
        {"in": ["cnv_sv/svdb_merge", ".merged.vcf"], "out": ["results/dna/cnv", ".merged.vcf"]},
        {"in": ["cnv_sv/svdb_query", ".svdb_query.vcf"], "out": ["results/dna/cnv", ".svdb_query.vcf"]},
    ]
    output_files = [
        "%s/%s_%s%s" % (file_info["out"][0], sample, unit_type, file_info["out"][1])
        for file_info in files
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
    ]
    input_files = [
        "%s/%s_%s%s" % (file_info["in"][0], sample, unit_type, file_info["in"][1])
        for file_info in files
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
    ]
    output_files += [
        "results/dna/vcf/%s_%s_%s.vcf.gz" % (caller, sample, t)
        for caller in ["mutect2", "vardict"]
        for sample in get_samples(samples)
        for t in get_unit_types(units, sample)
    ]
    input_files += [
        "snv_indels/%s/%s_%s.merged.vcf.gz" % (caller, sample, t)
        for caller in ["mutect2", "vardict"]
        for sample in get_samples(samples)
        for t in get_unit_types(units, sample)
    ]
    output_files += [
        "results/dna/cnv/%s_%s.manta_tumorSV.vcf.gz" % (sample, t)
        for sample in get_samples(samples)
        for t in get_unit_types(units, sample)
        if t == "T"
    ]
    input_files += [
        "cnv_sv/manta_run_workflow_t/%s/results/variants/tumorSV.vcf.gz" % (sample)
        for sample in get_samples(samples)
        for t in get_unit_types(units, sample)
        if t == "T"
    ]
    # output_files += [
    #     "results/dna/optitype/%s_%s.hla_type_result.tsv" % (sample, t)
    #     for sample in get_samples(samples)
    #     for t in get_unit_types(units, sample)
    # ]
    # input_files += [
    #     "biomarker/optitype/%s_%s/%s_%s_hla_type_result.tsv" % (sample, t, sample, t)
    #     for sample in get_samples(samples)
    #     for t in get_unit_types(units, sample)
    # ]
    output_files.append("results/dna/qc/multiqc.html")
    input_files.append("qc/multiqc/multiqc.html")
    return input_files, output_files


def compile_output_list(wildcards):
    return output_files


def get_flowcell(units, wildcards):
    flowcells = set([u.flowcell for u in get_units(units, wildcards)])
    if len(flowcells) > 1:
        raise ValueError("Sample type combination from different sequence flowcells")
    return flowcells.pop()


input_files, output_files = compile_result_file_list()
