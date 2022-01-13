# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas A"
__copyright__ = "Copyright 2021, Jonas A"
__email__ = "jonas.almlof@igp.uu.se"
__license__ = "GPL-3"

output_files = [
    "references/cnvkit_build_normal_reference/cnvkit.PoN.cnn",
    "references/create_read_count_panel_of_normals/gatk_cnv_panel_of_normal.hdf5",
    "references/msisensor_pro_baseline/Msisensor_pro_reference.list_baseline",
    "references/create_background_file/background_panel.tsv",
    "references/create_artifact_file/artifact_panel.tsv",
    "references/svdb_export/svdb_cnv.vcf",
]
output_files = [
    "results/cnvkit.PoN.cnn",
    "results/gatk_cnv_panel_of_normal.hdf5",
    "results/Msisensor_pro_reference.list_baseline",
    "results/background_panel.tsv",
    "results/artifact_panel.tsv",
    "results/svdb_cnv.vcf",
]

rule copy_cnvkit:
    input:
        "references/cnvkit_build_normal_reference/cnvkit.PoN.cnn",
    output:
        "results/cnvkit.PoN.cnn",
    shell:
        "cp {input} {output}"


rule copy_gatk_cnv:
    input:
        "references/create_read_count_panel_of_normals/gatk_cnv_panel_of_normal.hdf5",
    output:
        "results/gatk_cnv_panel_of_normal.hdf5",
    shell:
        "cp {input} {output}"


rule copy_msisensor_pro:
    input:
        "references/msisensor_pro_baseline/Msisensor_pro_reference.list_baseline",
    output:
        "results/Msisensor_pro_reference.list_baseline",
    shell:
        "cp {input} {output}"


rule copy_background:
    input:
        "references/create_background_file/background_panel.tsv",
    output:
        "results/background_panel.tsv",
    shell:
        "cp {input} {output}"


rule copy_artifact:
    input:
        "references/create_artifact_file/artifact_panel.tsv",
    output:
        "results/artifact_panel.tsv",
    shell:
        "cp {input} {output}"


rule copy_svdb_cnv:
    input:
        "references/svdb_export/svdb_cnv.vcf",
    output:
        "results/svdb_cnv.vcf",
    shell:
        "cp {input} {output}"
