# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas A"
__copyright__ = "Copyright 2021, Jonas A"
__email__ = "jonas.almlof@igp.uu.se"
__license__ = "GPL-3"


include: "rules/common_references.smk"
include: "rules/result_files_references.smk"

rule all:
    input:
        unpack(compile_output_list),

report: "report/workflow.rst"

module references:
   snakefile: github("hydra-genetics/references", path="workflow/Snakefile", tag="svdb")
   config: config

use rule * from references as references_*
