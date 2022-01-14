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


def compile_output_list(wildcards: snakemake.io.Wildcards):
    output_files = [
        "results/dna/bam/%s_%s.bam" % (sample, type)
        for sample in get_samples(samples)
        for type in get_unit_types(units, sample)
    ]
    output_files.append(
        [
            "results/dna/vcf/%s_%s_%s.vcf.gz" % (caller, sample, t)
            for caller in ["mutect2", "vardict"]
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
        ]
    )
    output_files.append(
        [
            "results/dna/vcf/%s_%s.ensembled.vcf.gz" % (sample, t)
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
        ]
    )
    output_files.append(
        [
            "results/dna/gvcf/%s_%s.gvcf.gz" % (sample, t)
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
        ]
    )
    output_files.append(
        [
            "results/dna/vcf/%s_%s.ensembled.vep_annotated.vcf" % (sample, t)
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
        ]
    )
    output_files.append(
        [
            "results/dna/vcf/%s_%s.ensembled.vep_annotated.filtered.codon_snvs.vcf.gz" % (sample, t)
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
    output_files.append(
        [
            "results/dna/qc/%s_%s.duplication_metrics.txt" % (sample, t)
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
        ]
    )
    output_files.append(
        [
            "results/dna/qc/%s_%s.alignment_summary_metrics.txt" % (sample, t)
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
        ]
    )
    output_files.append(
        [
            "results/dna/qc/%s_%s.HsMetrics.txt" % (sample, t)
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
        ]
    )
    output_files.append(
        [
            "results/dna/qc/%s_%s.insert_size_metrics.txt" % (sample, t)
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
        ]
    )
    output_files.append(
        [
            "results/dna/qc/%s_%s.samtools-stats.txt" % (sample, t)
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
        ]
    )
    output_files.append(
        [
            "results/dna/hotspot_info/%s_%s.hotspot_info.tsv" % (sample, t)
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
        ]
    )
    output_files.append(
        [
            "results/dna/msi/%s_%s.msisensor_pro.tsv" % (sample, t)
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
        ]
    )
    output_files.append(
        [
            "results/dna/tmb/%s_%s.TMB.txt" % (sample, t)
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
        ]
    )
    output_files.append(
        [
            "results/dna/hrd/%s_%s.hrd_score.txt" % (sample, t)
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
        ]
    )
    output_files.append(
        [
            "results/dna/fusions/%s_%s_gene_fuse_fusions.txt" % (sample, t)
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
        ]
    )
    output_files.append(
        [
            "results/dna/cnv/%s_%s.cnvkit_loh.cns" % (sample, t)
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
        ]
    )
    output_files.append(
        [
            "results/dna/cnv/%s_%s.gatk_cnv.seg" % (sample, t)
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
        ]
    )
    output_files.append(
        ["results/dna/cnv/%s_%s.gatk_cnv.vcf" % (sample, t) for sample in get_samples(samples) for t in get_unit_types(units, sample)]
    )
    output_files.append(
        ["results/dna/cnv/%s_%s.cnvkit.vcf" % (sample, t) for sample in get_samples(samples) for t in get_unit_types(units, sample)]
    )
    output_files.append(
        ["results/dna/cnv/%s_%s.png" % (sample, t) for sample in get_samples(samples) for t in get_unit_types(units, sample)]
    )
    output_files.append(
        [
            "results/dna/cnv/%s_%s.merged.vcf" % (sample, t)
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
        ]
    )
    output_files.append(
        [
            "results/dna/cnv/%s_%s.svdb_query.vcf" % (sample, t)
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
        ]
    )
    output_files.append("results/dna/qc/MultiQC.html")
    return output_files
