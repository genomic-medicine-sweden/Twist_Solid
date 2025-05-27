__author__ = "Jonas A"
__copyright__ = "Copyright 2021, Jonas A"
__email__ = "jonas.almlof@igp.uu.se"
__license__ = "GPL-3"

import json
import os
import pandas as pd
import re
from datetime import datetime
from snakemake.utils import validate
from snakemake.utils import min_version
import yaml

from hydra_genetics.utils.misc import get_module_snakefile
from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *
from hydra_genetics.utils.misc import extract_chr
from hydra_genetics.utils.misc import replace_dict_variables
from hydra_genetics import min_version as hydra_min_version

from hydra_genetics.utils.misc import export_config_as_file
from hydra_genetics.utils.software_versions import add_version_files_to_multiqc
from hydra_genetics.utils.software_versions import add_software_version_to_config
from hydra_genetics.utils.software_versions import export_pipeline_version_as_file
from hydra_genetics.utils.software_versions import export_software_version_as_file
from hydra_genetics.utils.software_versions import get_pipeline_version
from hydra_genetics.utils.software_versions import touch_pipeline_version_file_name
from hydra_genetics.utils.software_versions import touch_software_version_file
from hydra_genetics.utils.software_versions import use_container


hydra_min_version("3.0.0")

min_version("7.13.0")

### Set and validate config file

if not workflow.overwrite_configfiles:
    sys.exit("At least one config file must be passed using --configfile/--configfiles, by command line or a profile!")


config = replace_dict_variables(config)

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

with open(config["output"]) as output:
    if output.name.endswith(".json"):
        output_spec = json.load(output)
    elif output.name.endswith(".yaml") or output.name.endswith(".yml"):
        output_spec = yaml.safe_load(output)
    else:
        raise ValueError(f"output specification should be JSON or YAML: {output.name}")

validate(output_spec, schema="../schemas/output_files.schema.yaml")

## get version information on pipeline, containers and software

pipeline_name = "Twist_Solid"
pipeline_version = get_pipeline_version(workflow, pipeline_name=pipeline_name)
version_files = touch_pipeline_version_file_name(
    pipeline_version, date_string=pipeline_name, directory="results/versions/software"
)
if use_container(workflow):
    version_files.append(touch_software_version_file(config, date_string=pipeline_name, directory="results/versions/software"))
add_version_files_to_multiqc(config, version_files)


onstart:
    export_pipeline_version_as_file(pipeline_version, date_string=pipeline_name, directory="results/versions/software")
    if use_container(workflow):
        update_config, software_info = add_software_version_to_config(config, workflow, False)
        export_software_version_as_file(software_info, date_string=pipeline_name, directory="results/versions/software")
    export_config_as_file(update_config, date_string=pipeline_name, directory="results/versions")


### Set wildcard constraints
wildcard_constraints:
    barcode="[A-Z+]+",
    chr="[^_]+",
    flowcell="[A-Z0-9]+",
    lane="L[0-9]+",
    sample="|".join(get_samples(samples)),
    type="N|T|R",


if config.get("subsample", None) == "seqtk":
    merged_input_arriba = lambda wildcards: expand(
        "prealignment/seqtk_subsample_arriba/{{sample}}_{{type}}_{flowcell_lane_barcode}_{{read}}.ds.fastq.gz",
        flowcell_lane_barcode=[
            "{}_{}_{}".format(unit.flowcell, unit.lane, unit.barcode) for unit in get_units(units, wildcards, wildcards.type)
        ],
    )
else:
    merged_input_arriba = lambda wildcards: expand(
        "prealignment/fastp_pe_arriba/{{sample}}_{{type}}_{flowcell_lane_barcode}_{{read}}.fastq.gz",
        flowcell_lane_barcode=[
            "{}_{}_{}".format(unit.flowcell, unit.lane, unit.barcode) for unit in get_units(units, wildcards, wildcards.type)
        ],
    )


def compile_output_list(wildcards):
    output_files = []
    types = set([unit.type for unit in units.itertuples()])
    dedup_types = set([])
    for filedef in output_spec["files"]:
        output_files += set(
            [
                filedef["output"].format(sample=sample, type=unit_type, caller=caller)
                for sample in samples.index
                if "deduplication" not in filedef or samples.loc[sample].get("deduplication", "") in filedef["deduplication"]
                for unit_type in get_unit_types(units, sample)
                if unit_type in set(filedef["types"]).intersection(types)
                for caller in config["bcbio_variation_recall_ensemble"]["callers"]
            ]
        )
    return list(set(output_files))


def get_deduplication_bam_input(wildcards):
    sample = get_sample(samples, wildcards)
    if sample.get("deduplication", "") == "umi":
        return "alignment/bwa_mem_realign_consensus_reads/{sample}_{type}.umi.bam"
    else:
        return "alignment/samtools_merge_bam/{sample}_{type}.bam"


def get_deduplication_bam_input_manta(wildcards):
    sample = get_sample(samples, wildcards)
    if sample.get("deduplication", "") == "umi":
        return "alignment/bwa_mem_realign_consensus_reads/{sample}_T.umi.bam"
    else:
        return "alignment/samtools_merge_bam/{sample}_T.bam"


def get_deduplication_bam_chr_input(wildcards):
    sample = get_sample(samples, wildcards)
    if sample.get("deduplication", "") == "umi":
        return "alignment/samtools_extract_reads_umi/{sample}_{type}_{chr}.umi.bam"
    else:
        return "alignment/picard_mark_duplicates/{sample}_{type}_{chr}.bam"


def get_vardict_min_af(wildcards):
    sample = get_sample(samples, wildcards)
    if sample.get("deduplication", "") == "umi":
        return config.get("vardict", {}).get("allele_frequency_threshold_umi", "0.001")
    else:
        return config.get("vardict", {}).get("allele_frequency_threshold", "0.01")


def get_flowcell(units, wildcards):
    flowcells = set([u.flowcell for u in get_units(units, wildcards)])
    if len(flowcells) > 1:
        raise ValueError("Sample type combination from different sequence flowcells")
    return flowcells.pop()


def get_tc(wildcards):
    tc_method = wildcards.tc_method
    if tc_method == "pathology_purecn":
        tc = ""
        tc_file = f"cnv_sv/purecn_purity_file/{wildcards.sample}_{wildcards.type}.purity.txt"
        if os.path.exists(tc_file):
            with open(tc_file) as f:
                tc = f.read().strip()
        if tc == "" or float(tc) < 0.35:
            return get_sample(samples, wildcards)["tumor_content"]
        else:
            return tc
    elif tc_method == "pathology":
        return get_sample(samples, wildcards)["tumor_content"]
    else:
        tc_file = f"cnv_sv/purecn_purity_file/{wildcards.sample}_{wildcards.type}.purity.txt"
        if not os.path.exists(tc_file):
            return -1
        else:
            with open(tc_file) as f:
                tc = f.read().strip()
                if tc == "":
                    return "0.2"
                else:
                    return tc


def get_tc_file(wildcards):
    tc_method = wildcards.tc_method
    if tc_method == "pathology_purecn":
        return [f"cnv_sv/purecn_purity_file/{wildcards.sample}_{wildcards.type}.purity.txt", "samples.tsv"]
    elif tc_method == "pathology":
        return "samples.tsv"
    else:
        return f"cnv_sv/purecn_purity_file/{wildcards.sample}_{wildcards.type}.purity.txt"


def generate_star_read_group(wildcards):
    return "-R '@RG\\tID:{}\\tSM:{}\\tPL:{}\\tPU:{}\\tLB:{}' -v 1 ".format(
        "{}_{}".format(wildcards.sample, wildcards.type),
        "{}_{}".format(wildcards.sample, wildcards.type),
        "Illumina",
        "{}_{}".format(wildcards.sample, wildcards.type),
        "{}_{}".format(wildcards.sample, wildcards.type),
    )


def generate_copy_code(workflow, output_spec):
    code = ""
    for filedef in output_spec["files"]:
        if filedef["input"] is None:
            continue

        input_file = filedef["input"]
        output_file = filedef["output"]
        rule_name = "_copy_{}".format("_".join(re.sub(r"[\"'-.,]", "", filedef["name"].strip().lower()).split()))
        result_file = os.path.basename(filedef["output"])

        mem_mb = config.get("_copy", {}).get("mem_mb", config["default_resources"]["mem_mb"])
        mem_per_cpu = config.get("_copy", {}).get("mem_mb", config["default_resources"]["mem_mb"])
        partition = config.get("_copy", {}).get("partition", config["default_resources"]["partition"])
        threads = config.get("_copy", {}).get("threads", config["default_resources"]["threads"])
        time = config.get("_copy", {}).get("time", config["default_resources"]["time"])
        copy_container = config.get("_copy", {}).get("container", config["default_container"])

        code += f'@workflow.rule(name="{rule_name}")\n'
        code += f'@workflow.input("{input_file}")\n'
        code += f'@workflow.output("{output_file}")\n'
        code += f'@workflow.log("logs/{rule_name}_{result_file}.log")\n'
        code += f'@workflow.container("{copy_container}")\n'
        code += f'@workflow.resources(time = "{time}", threads = {threads}, mem_mb = {mem_mb}, mem_per_cpu = {mem_per_cpu}, partition = "{partition}")\n'
        code += '@workflow.shellcmd("cp --preserve=timestamps {input} {output}")\n\n'
        code += "@workflow.run\n"
        code += (
            f"def __rule_{rule_name}(input, output, params, wildcards, threads, resources, log, version, rule, "
            "conda_env, container_img, singularity_args, use_singularity, env_modules, bench_record, jobid, is_shell, "
            "bench_iteration, cleanup_scripts, shadow_dir, edit_notebook, conda_base_path, basedir, runtime_sourcecache_path, "
            "__is_snakemake_rule_func=True):\n"
            '\tshell ( "(cp --preserve=timestamps {input[0]} {output[0]}) &> {log}" , bench_record=bench_record, bench_iteration=bench_iteration)\n\n'
        )

    exec(compile(code, "result_to_copy", "exec"), workflow.globals)


generate_copy_code(workflow, output_spec)
