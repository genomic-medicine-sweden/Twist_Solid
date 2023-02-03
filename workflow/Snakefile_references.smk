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


module references:
    snakefile:
        github("hydra-genetics/references", path="workflow/Snakefile", tag="f6051df")
    config:
        config


use rule * from references as references_*


module misc:
    snakefile:
        get_module_snakefile(config, "hydra-genetics/misc", path="workflow/Snakefile", tag="v0.2.0")
    config:
        config


use rule tabix from misc as misc_tabix


use rule bgzip from misc as misc_bgzip
