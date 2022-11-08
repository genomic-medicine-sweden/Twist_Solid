# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas A"
__copyright__ = "Copyright 2021, Jonas A"
__email__ = "jonas.almlof@igp.uu.se"
__license__ = "GPL-3"

import json
import os
import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *
from hydra_genetics.utils.misc import extract_chr
from hydra_genetics import min_version as hydra_min_version

hydra_min_version("0.14.1")

min_version("7.13.0")

### Set and validate config file

if not workflow.overwrite_configfiles:
    sys.exit("At least one config file must be passed using --configfile/--configfiles, by command line or a profile!")


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
with open(config["output"]) as output:
    output_json = json.load(output)


wildcard_constraints:
    barcode="[A-Z+]+",
    chr="[^_]+",
    flowcell="[A-Z0-9]+",
    lane="L[0-9]+",
    sample="|".join(get_samples(samples)),
    type="N|T|R",


<<<<<<< HEAD
=======
def compile_result_file_list():
    dna_files = [
        {"in": ["alignment/samtools_merge_bam", ".bam"], "out": ["bam_dna", ".bam"]},
        {"in": ["alignment/samtools_merge_bam", ".bam.bai"], "out": ["bam_dna", ".bam.bai"]},
        {"in": ["snv_indels/gatk_mutect2_merge", ".bam"], "out": ["bam_dna/mutect2_indel_bam", ".bam"]},
        {"in": ["snv_indels/gatk_mutect2_merge", ".bam.bai"], "out": ["bam_dna/mutect2_indel_bam", ".bam.bai"]},
        {
            "in": ["annotation/background_annotation", ".background_annotation.vcf.gz"],
            "out": ["results/dna/vcf", ".annotated.vcf.gz"],
        },
        {
            "in": ["annotation/background_annotation", ".background_annotation.include.exon.filter.snv_soft_filter.vcf"],
            "out": ["results/dna/vcf", ".annotated.exon_only.filter.soft_filter.vcf"],
        },
        {
            "in": ["annotation/background_annotation", ".background_annotation.include.exon.filter.snv_hard_filter.vcf"],
            "out": ["results/dna/vcf", ".annotated.exon_only.filter.hard_filter.vcf"],
        },
        {
            "in": [
                "annotation/add_multi_snv_in_codon",
                ".background_annotation.include.exon.filter.snv_hard_filter.codon_snvs.sorted.vep_annotated.vcf",
            ],
            "out": ["results/dna/vcf", ".annotated.exon_only.filter.hard_filter.codon_snv.vcf"],
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
        {"in": ["qc/add_mosdepth_coverage_to_gvcf", ".mosdepth.g.vcf.gz"], "out": ["gvcf_dna", ".mosdepth.g.vcf.gz"]},
        {"in": ["qc/hotspot_report", ".output.tsv"], "out": ["results/dna/qc", ".coverage_and_mutations.tsv"]},
        {"in": ["qc/gatk_calculate_contamination", ".contamination.table"], "out": ["results/dna/qc", ".contamination.table"]},
        {"in": ["biomarker/msisensor_pro", ""], "out": ["results/dna/msi", ".msisensor_pro.score.tsv"]},
        {"in": ["biomarker/tmb", ".TMB.txt"], "out": ["results/dna/tmb", ".TMB.txt"]},
        {"in": ["biomarker/scarhrd", ".scarhrd_cnvkit_score.txt"], "out": ["results/dna/hrd", ".scarhrd_cnvkit_score.txt"]},
        {"in": ["fusions/gene_fuse", "_gene_fuse_fusions.txt"], "out": ["results/dna/fusions", ".gene_fuse_fusions.txt"]},
        {"in": ["fusions/report_gene_fuse", ".gene_fuse_report.tsv"], "out": ["results/dna/fusions", ".gene_fuse_report.tsv"]},
    ]
    output_files = [
        "%s/%s_%s%s" % (file_info["out"][0], sample, unit_type, file_info["out"][1])
        for file_info in dna_files
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
        if unit_type != "R"
    ]
    input_files = [
        "%s/%s_%s%s" % (file_info["in"][0], sample, unit_type, file_info["in"][1])
        for file_info in dna_files
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
        if unit_type != "R"
    ]
    dna_files2 = [
        {"in": ["cnv_sv/cnvkit_call", ".loh.cns"], "out": ["results/dna/cnv", ".cnvkit.loh.cns"]},
        {
            "in": ["cnv_sv/gatk_call_copy_ratio_segments", ".clean.calledCNVs.seg"],
            "out": ["results/dna/cnv", ".gatk_cnv.seg"],
        },
        {"in": ["cnv_sv/cnv_html_report", ".cnv.html"], "out": ["results/dna/cnv", ".cnv.html"]},
        {"in": ["cnv_sv/cnvkit_scatter", ".png"], "out": ["results/dna/cnv", ".cnvkit.scatter.png"]},
        {"in": ["cnv_sv/cnvkit_diagram", ".pdf"], "out": ["results/dna/cnv", ".cnvkit.diagram.pdf"]},
        {"in": ["cnv_sv/svdb_query", ".svdb_query.vcf"], "out": ["results/dna/cnv", ".svdb_query.vcf"]},
        {
            "in": ["cnv_sv/svdb_query", ".svdb_query.annotate_cnv.cnv_amp_genes.filter.cnv_hard_filter_amp.vcf"],
            "out": ["results/dna/cnv", ".cnv_hard_filter_amp.vcf"],
        },
        {
            "in": ["cnv_sv/svdb_query", ".svdb_query.annotate_cnv.cnv_loh_genes.filter.cnv_hard_filter_loh.vcf"],
            "out": ["results/dna/cnv", ".cnv_hard_filter_loh.vcf"],
        },
        {"in": ["cnv_sv/svdb_query", ".cnv_report.tsv"], "out": ["results/dna/cnv", ".cnv_report.tsv"]},
    ]
    output_files = [
        "%s/%s_%s/%s_%s%s" % (file_info["out"][0], sample, unit_type, sample, unit_type, file_info["out"][1])
        for file_info in dna_files2
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
        if unit_type != "R"
    ]
    input_files = [
        "%s/%s_%s%s" % (file_info["in"][0], sample, unit_type, file_info["in"][1])
        for file_info in dna_files2
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
        if unit_type != "R"
    ]
    output_files += [
        "results/dna/vcf/%s_%s_%s.vcf.gz" % (caller, sample, unit_type)
        for caller in ["gatk_mutect2", "vardict"]
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
        if unit_type != "R"
    ]
    input_files += [
        "snv_indels/%s/%s_%s.merged.vcf.gz" % (caller, sample, unit_type)
        for caller in ["gatk_mutect2", "vardict"]
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
        if unit_type != "R"
    ]
    output_files += [
        "results/dna/cnv/%s_%s/%s_%s.manta_tumorSV.vcf.gz" % (sample, unit_type, sample, unit_type)
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
        if unit_type == "T"
    ]
    input_files += [
        "cnv_sv/manta_run_workflow_t/%s/results/variants/tumorSV.vcf.gz" % (sample)
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
        if unit_type == "T"
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

    rna_files = [
        {"in": ["fusions/arriba", ".fusions.tsv"], "out": ["results/rna/fusion", ".arriba.fusions.tsv"]},
        {"in": ["fusions/arriba_draw_fusion", ".pdf"], "out": ["results/rna/fusion", ".arriba.fusions.pdf"]},
        {"in": ["fusions/report_fusions", ".fusion_report.tsv"], "out": ["results/rna/fusion", ".fusion_report.tsv"]},
        {"in": ["fusions/exon_skipping", ".results.tsv"], "out": ["results/rna/fusion", ".exon_skipping.tsv"]},
        {
            "in": ["qc/house_keeping_gene_coverage", ".house_keeping_gene_coverage.tsv"],
            "out": ["results/rna/qc", ".house_keeping_gene_coverage.tsv"],
        },
        {"in": ["snv_indels/bcftools_id_snps", ".id_snps.vcf"], "out": ["results/rna/id_snps", ".id_snps.vcf"]},
    ]
    output_files += [
        "%s/%s_%s%s" % (file_info["out"][0], sample, unit_type, file_info["out"][1])
        for file_info in rna_files
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
        if unit_type == "R"
    ]
    input_files += [
        "%s/%s_%s%s" % (file_info["in"][0], sample, unit_type, file_info["in"][1])
        for file_info in rna_files
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
        if unit_type == "R"
    ]

    rna_files2 = [
        {
            "in": ["fusions/star_fusion", "star-fusion.fusion_predictions.tsv"],
            "out": ["results/rna/fusion", ".star-fusion.fusion_predictions.tsv"],
        },
        {"in": ["fusions/star_fusion", "Aligned.out.sorted.bam"], "out": ["bam_rna", ".star_fusion.bam"]},
        {"in": ["fusions/star_fusion", "Aligned.out.sorted.bam.bai"], "out": ["bam_rna", ".star_fusion.bam.bai"]},
        {
            "in": ["fusions/fusioncatcher", "final-list_candidate-fusion-genes.hg19.txt"],
            "out": ["results/rna/fusion", ".fusioncatcher.fusion_predictions.txt"],
        },
    ]
    output_files += [
        "%s/%s_%s%s" % (file_info["out"][0], sample, unit_type, file_info["out"][1])
        for file_info in rna_files2
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
        if unit_type == "R"
    ]
    input_files += [
        "%s/%s_%s/%s" % (file_info["in"][0], sample, unit_type, file_info["in"][1])
        for file_info in rna_files2
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
        if unit_type == "R"
    ]

    types = set([unit.type for unit in units.itertuples()])
    if "R" in types:
        output_files.append("results/rna/qc/multiqc_RNA.html")
        input_files.append("qc/multiqc/multiqc_RNA.html")
    if not set(["N", "T"]).isdisjoint(types):
        output_files.append("results/dna/qc/multiqc_DNA.html")
        input_files.append("qc/multiqc/multiqc_DNA.html")
    return input_files, output_files


>>>>>>> feat: added merged mutect2 bam files and bai. Also put cnv output in one folder per sample.
def compile_output_list(wildcards):
    output_files = []
    types = set([unit.type for unit in units.itertuples()])
    for output in output_json:
        output_files += set(
            [
                output.format(sample=sample, type=unit_type, caller=caller)
                for sample in get_samples(samples)
                for unit_type in get_unit_types(units, sample)
                if unit_type in set(output_json[output]["types"]).intersection(types)
                for caller in config["bcbio_variation_recall_ensemble"]["callers"]
            ]
        )
    return list(set(output_files))


def get_flowcell(units, wildcards):
    flowcells = set([u.flowcell for u in get_units(units, wildcards)])
    if len(flowcells) > 1:
        raise ValueError("Sample type combination from different sequence flowcells")
    return flowcells.pop()


def generate_copy_code(workflow, output_json):
    code = ""
    for result, values in output_json.items():
        if values["file"] is not None:
            input_file = values["file"]
            output_file = result
            rule_name = values["name"]
            mem_mb = config.get('_copy', {}).get("mem_mb", config["default_resources"]["mem_mb"])
            mem_per_cpu = config.get('_copy', {}).get("mem_mb", config["default_resources"]["mem_mb"])
            partition = config.get("_copy", {}).get("partition", config["default_resources"]["partition"])
            threads = config.get("_copy", {}).get("threads", config["default_resources"]["threads"])
            time = config.get("_copy", {}).get("time", config["default_resources"]["time"])
            copy_container = config.get("_copy", {}).get("container", config["default_container"])
            result_file = os.path.basename(output_file)
            code += f'@workflow.rule(name="{rule_name}")\n'
            code += f'@workflow.input("{input_file}")\n'
            code += f'@workflow.output("{output_file}")\n'
            code += f'@workflow.log("logs/{rule_name}_{result_file}.log")\n'
            code += f'@workflow.container("{copy_container}")\n'
            code += f'@workflow.conda("../env/copy_result.yaml")\n'
            code += f'@workflow.resources(time = "{time}", threads = {threads}, mem_mb = {mem_mb}, mem_per_cpu = {mem_per_cpu}, partition = "{partition}")\n'
            code += '@workflow.shellcmd("cp {input} {output}")\n\n'
            code += "@workflow.run\n"
            code += (
                f"def __rule_{rule_name}(input, output, params, wildcards, threads, resources, log, version, rule, "
                "conda_env, container_img, singularity_args, use_singularity, env_modules, bench_record, jobid, is_shell, "
                "bench_iteration, cleanup_scripts, shadow_dir, edit_notebook, conda_base_path, basedir, runtime_sourcecache_path, "
                "__is_snakemake_rule_func=True):\n"
                '\tshell ( "(cp {input[0]} {output[0]}) &> {log}" , bench_record=bench_record, bench_iteration=bench_iteration)\n\n'
            )
    exec(compile(code, "result_to_copy", "exec"), workflow.globals)


generate_copy_code(workflow, output_json)
