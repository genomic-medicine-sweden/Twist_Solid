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

min_version("6.8.0")

### Set and validate config file


configfile: "config.yaml"


validate(config, schema="../schemas/config.schema.yaml")
config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")


### Read and validate samples file

samples = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

### Read and validate units file

units = pandas.read_table(config["units"], dtype=str).set_index(["sample", "type"], drop=False)
validate(units, schema="../schemas/units.schema.yaml")

### Set wildcard constraints


def get_bams(units: pandas.DataFrame) -> typing.List[str]:
    """
    function used to extract all bam files found in units.tsv
    Args:
        units: DataFrame generate by importing a file following schema defintion
               found in pre-alignment/workflow/schemas/units.schema.tsv
    Returns:
        List of strings with all bam file names and path
    """
    return list(set([unit.bam for unit in units.itertuples()]))


def get_gvcfs(units: pandas.DataFrame) -> typing.List[str]:
    """
    function used to extract all gvcf files found in units.tsv
    Args:
        units: DataFrame generate by importing a file following schema defintion
               found in pre-alignment/workflow/schemas/units.schema.tsv
    Returns:
        List of strings with all gvcf file names and path
    """
    return list(set([unit.gvcf for unit in units.itertuples()]))


def get_vcfs(units: pandas.DataFrame) -> typing.List[str]:
    """
    function used to extract all vcf files found in units.tsv
    Args:
        units: DataFrame generate by importing a file following schema defintion
               found in pre-alignment/workflow/schemas/units.schema.tsv
    Returns:
        List of strings with all vcf file names and path
    """
    return list(set([unit.vcf for unit in units.itertuples()]))


wildcard_constraints:
    sample="|".join(samples.index),
    unit="N|T|R",


def compile_output_list(wildcards: snakemake.io.Wildcards):
    return [
        "results/dna/bam/%s_%s.bam" % (sample, type)
        for sample in get_samples(samples)
        for type in get_unit_types(units, sample)
    ]
