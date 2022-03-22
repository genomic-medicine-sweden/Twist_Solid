# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas A"
__copyright__ = "Copyright 2021, Jonas A"
__email__ = "jonas.almlof@igp.uu.se"
__license__ = "GPL-3"


include: "rules/common.smk"
include: "rules/hotspot_report.smk"
include: "rules/result_files.smk"

rule all:
    input:
        unpack(compile_output_list),

ruleorder: snv_indels_tabix_vcf > misc_tabix
ruleorder: snv_indels_mutect2 > misc_tabix
ruleorder: snv_indels_mutect2_gvcf > misc_tabix
ruleorder: snv_indels_bgzip_vcf > misc_bgzip
ruleorder: annotation_bgzip_vcf > misc_bgzip
ruleorder: annotation_tabix_vcf > misc_tabix
ruleorder: annotation_add_mosdepth_coverage_to_gvcf > misc_bgzip

report: "report/workflow.rst"

module prealignment:
   snakefile: github("hydra-genetics/prealignment", path="workflow/Snakefile", tag="develop")
   config: config

use rule * from prealignment as prealignment_*

module alignment:
   snakefile: github("hydra-genetics/alignment", path="workflow/Snakefile", tag="develop")
   config: config

use rule * from alignment as alignment_*

module snv_indels:
   snakefile: github("hydra-genetics/snv_indels", path="workflow/Snakefile", tag="develop")
   config: config

use rule * from snv_indels as snv_indels_*

module annotation:
   snakefile: github("hydra-genetics/annotation", path="workflow/Snakefile", tag="develop")
   config: config

use rule * from annotation as annotation_*

module filtering:
   snakefile: github("hydra-genetics/filtering", path="workflow/Snakefile", tag="develop")
   config: config

use rule * from filtering as filtering_*

module qc:
   snakefile: github("hydra-genetics/qc", path="workflow/Snakefile", tag="develop")
   config: config

use rule * from qc as qc_*

module biomarker:
   snakefile: github("hydra-genetics/biomarker", path="workflow/Snakefile", tag="develop")
   config: config

use rule * from biomarker as biomarker_*

module fusions:
   snakefile: github("hydra-genetics/fusions", path="workflow/Snakefile", tag="develop")
   config: config

use rule * from fusions as fusions_*

module cnv_sv:
   snakefile: github("hydra-genetics/cnv_sv", path="workflow/Snakefile", tag="develop")
   config: config

use rule * from cnv_sv as cnv_sv_*

module misc:
   snakefile: github("hydra-genetics/misc", path="workflow/Snakefile", tag="develop")
   config: config

use rule tabix from misc as misc_tabix

use rule bgzip from misc as misc_bgzip
