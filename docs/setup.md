# Pipeline setup and configuration
There are a number of main files that governs how the pipeline is executed listed below:

* Snakefile
* common.smk
* config.yaml
* resources.yaml
* samples.tsv and units.tsv

There is more general information about the content of these files in hydra-genetics documentation in [code standards](https://hydra-genetics.readthedocs.io/en/latest/standards/), [config](https://hydra-genetics.readthedocs.io/en/latest/config/) and [Snakefile](https://hydra-genetics.readthedocs.io/en/latest/import/).

## Snakefile
The `Snakefile` is located in workflow/ and imports hydra-genetics modules and rules as well as modifies these rules when needed. It also imports pipeline specific rules and define rule orders. Finally, this is where the rule all is defined.

## common.smk
The `common.smk` is located under workflow/rules/. This is a general rule taking care of any actions that are not directly connected with running a specific program. It includes version checks, import of config, resources, tsv-files and validations using schemas. Functions used by pipeline specific rules are also defined here as well as the output files using the function **compile_output_list** which programmatically generates a list of all necessary output files for the module to be targeted in the all rule defined in the `Snakemake` file. See further [Result files](https://hydra-genetics.readthedocs.io/en/latest/results/).

## config.yaml
The `config.yaml` is located under config/. The file ties all file and other dependencies as well as parameters for different rules together.
See further [pipeline configuration](https://hydra-genetics.readthedocs.io/en/latest/config/).

## resources.yaml
The `resources.yaml` is located under config/. The file declares default resources used by rules as well as resources for specific rules that needs more resources than allocated by default. See further [pipeline configuration](https://hydra-genetics.readthedocs.io/en/latest/config/).

## samples.tsv and units.tsv
The `samples.tsv` and `units.tsv` are input files that must be generated before running the pipeline and should in general be located in the base folder of the analysis folder even if this can be set in the config.yaml. See further [running the pipeline](running.md) and [create input files](https://hydra-genetics.readthedocs.io/en/latest/create_sample_files/).
