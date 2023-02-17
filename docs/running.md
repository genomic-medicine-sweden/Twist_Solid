# Running the pipeline

## Installation
### Clone the Twist Solid git repo
```bash
git clone https://github.com/genomic-medicine-sweden/Twist_Solid.git
```

### Create python environment
To run the Twist Solid pipeline a python virtual environment is needed.
Create a new [python environment](https://docs.python.org/3/library/venv.html):
```bash
python3 -m venv /path/to/new/virtual/environment
```

### Install hydra-genetics
Activate the virtual environment and install hydra-genetics tools and other requirements specified in `requirements.txt`.
```bash
cd path/to/Twist_Solid/git/clone/
source /path/to/new/virtual/environment/bin/activate
pip install -r requirements.txt
pip install hydra_genetics
```

## Dependencies
The pipeline assumes that Snakemake and support for Singularity containers are installed.

## Input sample files
The pipeline uses sample input files (`samples.tsv` and `units.tsv`) with information regarding sample information, sequencing meta information as well as the location of the fastq-files. Using the python virtual environment created above it is possible to generate these files automatically using [hydra-genetics create-input-files](https://hydra-genetics.readthedocs.io/en/latest/create_sample_files/):
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
