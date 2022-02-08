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


def compile_output_list(wildcards):
    files = {
        "results/dna/vcf": [
            ".ensembled.vcf.gz",
            ".ensembled.vep_annotated.vcf",
            ".ensembled.vep_annotated.filtered.codon_snvs.vcf.gz",
            ".ensembled.vep_annotated.filtered.codon_snvs.nocnv.vcf.gz",
            ".ensembled.vep_annotated.filtered.codon_snvs.exon_only.vcf.gz",
        ],
        "results/dna/gvcf": [".gvcf.gz"],
        "results/dna/qc": [
            ".duplication_metrics.txt",
            ".alignment_summary_metrics.txt",
            ".HsMetrics.txt",
            ".insert_size_metrics.txt",
            ".samtools-stats.txt",
        ],
        "results/dna/msi": [".msisensor_pro.tsv"],
        "results/dna/tmb": [".TMB.txt"],
        "results/dna/hrd": [".hrd_score.txt"],
        "results/dna/fusions": ["_gene_fuse_fusions.txt"],
        "results/dna/cnv": [
            ".cnvkit_loh.cns",
            ".gatk_cnv.seg",
            ".gatk_cnv.vcf",
            ".cnvkit.vcf",
            ".merged.vcf",
        ],
    }
    output_files = [
        "%s/%s_%s%s" % (prefix, sample, unit_type, suffix)
        for prefix in files.keys()
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
        for suffix in files[prefix]
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
            "results/dna/qc/%s_%s_%s_fastqc.html" % (sample, t, read)
            for read in ["fastq1", "fastq2"]
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
        ]
    )
    output_files.append("results/dna/qc/MultiQC.html")
    return output_files
