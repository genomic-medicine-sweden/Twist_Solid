# Running the pipeline

## Requirements
**Recommended hardware**

- CPU: >10 cores per sample
- Memory: 6GB per core
- Storage: >75GB per sample

**Note**: Running the pipeline with less resources may work, but has not been tested.

**Software**

- [python](https://www.python.org/), version 3.9 or newer
- [pip3](https://pypi.org/project/pip/)
- [virtuelenv](https://docs.python.org/3/library/venv.html)
- [singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html)

**Nice to have**

- DRMAA compatible scheduler

## Installation
A list of releases of the Twist Solid pipeline can be found at: [Releases](https://github.com/genomic-medicine-sweden/Twist_Solid/releases).

### Clone the Twist Solid git repo
```bash
VERSION="v0.4.0"
git clone --branch ${VERSION} https://github.com/genomic-medicine-sweden/Twist_Solid.git
```

### Create python environment
To run the Twist Solid pipeline a python virtual environment is needed.
```bash
# Create a new virtual environment
python3 -m venv /path/to/new/virtual/environment
```

### Install pipeline requirements
Activate the virtual environment and install pipeline requirements specified in `requirements.txt`.
```bash
cd path/to/Twist_Solid/git/clone/
source /path/to/new/virtual/environment/bin/activate
pip install -r requirements.txt
pip install hydra_genetics
```

## Input sample files
The pipeline uses sample input files (`samples.tsv` and `units.tsv`) with information regarding sample information, sequencing meta information as well as the location of the fastq-files. Specification for the input files can be found at [Twist Solid schemas](https://github.com/genomic-medicine-sweden/Twist_Solid/blob/develop/workflow/schemas/). Using the python virtual environment created above it is possible to generate these files automatically using [hydra-genetics create-input-files](https://hydra-genetics.readthedocs.io/en/latest/create_sample_files/):
```bash
hydra-genetics create-input-files -d path/to/fastq-files/
```

## Run command
Using the activated python virtual environment created above, this is a basic command for running the pipeline:
```bash
snakemake --profile profiles/ -s workflow/Snakefile
```
The are many additional [snakemake running options](https://snakemake.readthedocs.io/en/stable/executing/cli.html#) some of which is listed below. However, options that are always used should be put in the [profile](https://hydra-genetics.readthedocs.io/en/latest/profile/).

* --notemp - Saves all intermediate files. Good for development and testing different options.
* --until <rule> - Runs only rules dependent on the specified rule.
