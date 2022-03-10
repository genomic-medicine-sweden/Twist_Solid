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
    run:
        import subprocess
        i = 0
        for file in input:
            subprocess.run(["cp", file, output[i]])
            i += 1
