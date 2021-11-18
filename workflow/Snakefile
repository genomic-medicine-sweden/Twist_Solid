# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas A"
__copyright__ = "Copyright 2021, Jonas A"
__email__ = "jonas.almlof@igp.uu.se"
__license__ = "GPL-3"


include: "rules/common.smk"
include: "rules/result_files.smk"


report: "report/workflow.rst"

module prealignment:
   snakefile: github("hydra-genetics/prealignment", path="workflow/Snakefile", tag="develop")
   config: config

use rule * from prealignment as prealignment_*

module alignment:
   snakefile: github("hydra-genetics/alignment", path="workflow/Snakefile", tag="develop")
   config: config

use rule * from alignment as alignment_*

rule all:
    input:
        unpack(compile_output_list),
