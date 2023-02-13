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

from hydra_genetics.utils.misc import get_module_snakefile
from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *
from hydra_genetics.utils.misc import extract_chr
from hydra_genetics import min_version as hydra_min_version

hydra_min_version("0.15.0")

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


def get_cnv_ratio_file(wildcards):
    caller = wildcards.get("caller", "")
    if caller == "cnvkit":
        return "cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cnr".format(**wildcards)
    elif caller == "gatk":
        return "cnv_sv/gatk_denoise_read_counts/{sample}_{type}.clean.denoisedCR.tsv".format(**wildcards)
    else:
        raise NotImplementedError(f"not implemented for caller {caller}")


def get_cnv_segment_file(wildcards):
    caller = wildcards.get("caller", "")
    if caller == "cnvkit":
        return "cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cns".format(**wildcards)
    elif caller == "gatk":
        return "cnv_sv/gatk_model_segments/{sample}_{type}.clean.cr.seg".format(**wildcards)
    else:
        raise NotImplementedError(f"not implemented for caller {caller}")


def get_json_for_merge_json(wildcards):
    json_dict = {}
    for v in config.get("svdb_merge", {}).get("tc_method"):
        tc_method = v["name"]
        callers = v["cnv_caller"]
        for caller in callers:
            if tc_method in json_dict:
                json_dict[tc_method].append(
                    f"cnv_sv/cnv_html_report/{wildcards.sample}_{wildcards.type}.{caller}.{tc_method}.json"
                )
            else:
                json_dict[tc_method] = [f"cnv_sv/cnv_html_report/{wildcards.sample}_{wildcards.type}.{caller}.{tc_method}.json"]
    return json_dict[wildcards.tc_method]


def get_filtered_cnv_vcfs_for_merge_json(wildcards):
    cnv_vcfs = []
    tags = config.get("cnv_html_report", {}).get("cnv_vcf", [])
    for t in tags:
        cnv_vcfs.append(
            f"cnv_sv/svdb_query/{wildcards.sample}_{wildcards.type}.{wildcards.tc_method}.svdb_query."
            f"annotate_cnv.{t['annotation']}.filter.{t['filter']}.vcf.gz"
        )
    return sorted(cnv_vcfs)


def get_tc(wildcards):
    tc_method = wildcards.tc_method
    if tc_method == "pathology":
        return get_sample(samples, wildcards)["tumor_content"]
    else:
        tc_file = f"cnv_sv/{tc_method}_purity_file/{wildcards.sample}_{wildcards.type}.purity.txt"
        if not os.path.exists(tc_file):
            return -1
        else:
            with open(tc_file) as f:
                return f.read()


def get_unfiltered_cnv_vcfs_for_merge_json(wildcards):
    cnv_vcfs = []
    tags = config.get("cnv_html_report", {}).get("cnv_vcf", [])
    for t in tags:
        cnv_vcfs.append(
            f"cnv_sv/svdb_query/{wildcards.sample}_{wildcards.type}.{wildcards.tc_method}.svdb_query."
            f"annotate_cnv.{t['annotation']}.vcf.gz"
        )
    return sorted(cnv_vcfs)


def generate_copy_code(workflow, output_json):
    code = ""
    for result, values in output_json.items():
        if values["file"] is not None:
            input_file = values["file"]
            output_file = result
            rule_name = values["name"]
            mem_mb = config.get("_copy", {}).get("mem_mb", config["default_resources"]["mem_mb"])
            mem_per_cpu = config.get("_copy", {}).get("mem_mb", config["default_resources"]["mem_mb"])
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
