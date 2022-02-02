# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas A"
__copyright__ = "Copyright 2021, Jonas A"
__email__ = "jonas.almlof@igp.uu.se"
__license__ = "GPL-3"


include: "rules/common.smk"
include: "rules/result_files.smk"

rule all:
    input:
        unpack(compile_output_list),

report: "report/workflow.rst"

module prealignment:
   snakefile: github("hydra-genetics/prealignment", path="workflow/Snakefile", tag="2bccef0")
   config: config

use rule * from prealignment as prealignment_*

module alignment:
   snakefile: github("hydra-genetics/alignment", path="workflow/Snakefile", tag="daf3747")
   config: config

use rule * from alignment as alignment_*

module snv_indels:
   snakefile: github("hydra-genetics/snv_indels", path="workflow/Snakefile", tag="36b9dfe")
   config: config

use rule * from snv_indels as snv_indels_*

module annotation:
   snakefile: github("hydra-genetics/annotation", path="workflow/Snakefile", tag="3953402")
   config: config

use rule * from annotation as annotation_*

module filtering:
   snakefile: github("hydra-genetics/filtering", path="workflow/Snakefile", tag="897af95")
   config: config

use rule * from filtering as filtering_*

module qc:
   snakefile: github("hydra-genetics/qc", path="workflow/Snakefile", tag="4ec562d")
   config: config

use rule * from qc as qc_*

module biomarker:
   snakefile: github("hydra-genetics/biomarker", path="workflow/Snakefile", tag="1ec1b3b")
   config: config

use rule * from biomarker as biomarker_*

module fusions:
   snakefile: github("hydra-genetics/fusions", path="workflow/Snakefile", tag="12a1f00")
   config: config

use rule * from fusions as fusions_*

module cnv_sv:
   snakefile: github("hydra-genetics/cnv_sv", path="workflow/Snakefile", tag="2462507")
   config: config

use rule * from cnv_sv as cnv_sv_*
