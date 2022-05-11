# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@igp.uu.se"
__license__ = "GPL-3"


rule copy_results_files:
    input:
        input_files,
    output:
        output_files,
    log:
        "logs/copy_result.log",
    resources:
        threads=config.get("copy_results_files", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("copy_results_files", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("copy_results_files", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("copy_results_files", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("copy_results_files", {}).get("partition", config["default_resources"]["partition"]),
    run:
        import subprocess

        i = 0
        for file in input:
            subprocess.run(["rsync", "--update", "-a", file, output[i]])
            i += 1
