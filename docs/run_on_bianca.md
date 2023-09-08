Follow the instructions [RUNNING THE PIPELINE -> Closed System](run_on_closed_env.md)

# Setup environment

1. Copy updated config folder to the working directory
2. Update bianca profile to match your projects
3. Create input files using: hydra-genetics create-input-files
4. Run the pipeline

**NOTE**: it's recommended to run your job using a slurm script, since the login node my close down due to inactivity.
# SBATCH

*run_pipeline.sh*
```bash
#!/bin/bash -l

#SBATCH -p core -n 1 -t 24:00:00
#SBATCH -A sensXXXX -J run_twist_solid

PATH_TO_ENV={PATH_TO_ENV}

PATH_TO_ANALYSIS_FOLDER={PATH_TO_ANALYSIS}

PATH_TO_TWIST_SOLID={PATH_TO_FOLDER_WITH_PIPELINE}

PATH_TO_UPDATED_PROFILE={PATH_TO_PROFILE}

cd $PATH_TO_ANALYSIS_FOLDER

module load slurm-drmaa 

source  ${PATH_TO_EXTRACTED_CONDA_ENV}/venv/bin/activate

snakemake -s ${PATH_TO_FOLDER_WITH_PIPELINE}/Twist_Solid/workflow/Snakefile  --profile ${PATH_TO_UPDATED_PROFILE}/bianca

```

# Example setup

## Folder structure

```bash
proj/sensYYYYXXX/nobackup/username
|---design_and_ref_files
|   |---GMS560
|   |   |---Artifact/
|   |   |---Background/
|   |   |---design/
|   |   |---PoN/
|---ref_data
|   |---arriba_v2.3.0
|   |---fuseq_wes
|   |---fusioncatcher
|   |---GNOMAD
|   |---hg19
|   |---refGene
|   |---star
|   |---star-fusion
|   |---vep
|---Twist_Solid_env
|   |---venv/ # conda enc
|   |---hydra-genetics/ # Modules
|   |---snakemake-wrappers/ #Wrappers
|   |---Twist Solid/  # Pipeline
|       |---profiles/bianca ¤ Update profile
|       |---config/ # Updated config
|       |---workflow/
|       |---Snakefile
|---analysis
|   |---samples.tsv
|   |---units.tsv
|   |---configs/ # Copied config
|---singularity_cache
```

Point to uploaded reference files
```yaml
# config/config.data.hg19.yaml
# Update the following lines:
PROJECT_DESIGN_DATA: "/proj/sensYYYYXXX/nobackup/username/design_and_ref_files"
PROJECT_PON_DATA: "/proj/sensYYYYXXX/nobackup/username/design_and_ref_files"
PROJECT_REF_DATA: "/proj/sensYYYYXXX/nobackup/username/design_and_ref_files"
```

**config.yaml**
```yaml
# Update the following line
hydra_local_path: "/proj/sensYYYYXXX/nobackup/username/Twist_Solid_env/hydra-genetics"

# All containers need to point to your singularity cache
# Ex
default_container: /proj/sensYYYYXXX/nobackup/username/singularity_cache/hydragenetics_common_0.1.9.sif
```

**profiles/bianaca/config.yaml**
```yaml
wrapper-prefix="git+file://proj/sensYYYXXX/nobackup/username/Twist_Solid_env/snakemake-wrappers"

drmaa: " -A sensYYYXXX -N 1-1 -t {resources.time} -n {resources.threads} --mem={resources.mem_mb} --mem-per-cpu={resources.mem_per_cpu} --mem-per-cpu={resources.mem_per_cpu} --partition={resources.partition} -J {rule} -e slurm_out/{rule}_%j.err -o slurm_out/{rule}_%j.out"
```