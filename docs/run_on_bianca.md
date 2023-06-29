Bianca is a cluster without access to internet which causes some problems for pipelines relying on resources found online. But it's easy to handle=)

# Preperations

```bash
# Set Twist Solid version
TAG_OR_BRANCH="v0.6.1"

# Clone selected version
git clone --branch ${VERSION} https://github.com/genomic-medicine-sweden/Twist_Solid.git
cd Twist_Solid
virtualenv venv && source venv/bin/activate
pip install -r requirements.txt
```

## Fetach resources

### Download reference files

```bash
# NextSeq
 hydra-genetics references download -o design_and_ref_files -v config/references/design_files.hg19.yaml -v config/references/nextseq.hg19.pon.yaml -v config/references/references.hg19.yaml

 #NovaSeq, not all files are prepare for novaseq
 hydra-genetics references download -o design_and_ref_files -v config/references/design_files.hg19.yaml -v config/references/novaseq.hg19.pon.yaml -v config/references/references.hg19.yaml
```

## Download  Containers
```bash
# NOTE: singularity command need to be available for this step
hydra-genetics prepare-environment create-singularity-files -c config/config.yaml -o singularity_cache 
```

## Environment

Create an environment, on a computer/server with access to internet, that can be moved to bianca.

### Conda version

Requires:

 - conda

```bash
# Build compressed file containing, named Twist_Solid_{TAG_OR_BRANCH}.tar.gz
# - Twist Solid Pipeline
# - snakemake-wrappers
# - hydra-genetics modules
# - conda env
TAG_OR_BRANCH="vX.Y.X" bash build/build_conda.sh
```


### Singularity version
Requires:

 - docker
 - singularity

```bash
# Build compressed file containing
# - Twist Solid Pipeline
# - snakemake-wrappers
# - hydra-genetics modules
# - conda env
TAG_OR_BRANCH="vX.Y.X" bash build/build_container.sh
```
 
---

# Bianca
The following file/folders need to be uploaded to bianca, with sftp:

1. design_and_ref_files
2. Twist_Solid_{TAG_OR_BRANCH}.tar.gz: for conda
3. Twist_Solid_{TAG_OR_BRANCH}.sif: for singularity
3. singularity_cache 

---

# On bianca

## Setup environment
### Conda

#### Unpack environment and activate
```bash
# Extract tar.
TAG_OR_BRANCH=develop
tar -xvf Twist_Solid_${TAG_OR_BRANCH}.tar.gz
cd Twist_Solid_{TAG_OR_BRANCH}
mkdir venv && tar xvf env.tar.gz -C venv/
source venv/bin/activate

# Variable that will be used lated
PATH_TO_ENV=${PWD}
PATH_TO_HYDRA_MODULES=${PATH_TO_EXTRACTED_TAR_GZ}/hydra-genetics
```

#### Create config and profile

*config files*

```bash
PATH_TO_design_and_ref_files="PATH_TO_UPLOADED_DESIGN_AND_REF_FILES"
PATH_TO_singularity_cache="PATH_TO_UPLOADED_SINGULARITY_CACHE"
# Conda environment still need to be active
# Prepare config
cd $PATH_TO_EXTRACTED_TAR_GZ/Twist_Solid
cp config/config.yaml config/config.yaml.copy
cp config/bianca/config.hg19.yaml config/bianca/config.hg19.yaml.copy

# Update design and ref files location and hydra-genetics module location
# Update hydra module location
# Make sure the environment still is active
hydra-genetics prepare-environment reference-path-update  -c config/bianca/config.hg19.yaml.copy -n config/bianca/config.hg19.yaml --reference-path /PROJECT_DATA:${PATH_TO_design_and_ref_files} --reference-path PATH_TO_REPO:${PATH_TO_HYDRA_MODULES}

# Make use of local singularities: ex /proj/sens2022566/nobackup/patriksm/singularity_cache
hydra-genetics prepare-environment container-path-update -c config/config.yaml.copy -n config/config.yaml -p ${PATH_TO_singularity_cache}

# Ex: updating configs
#  - ref files att /proj/sens2022566/nobackup/patriksm/design_and_ref_files
#  - module at /proj/sens2022566/nobackup/patriksm/Twist_Solid_{TAG_OR_BRANCH}/hydra-genetics
# command: hydra-genetics prepare-environment reference-path-update -c config.bianca.hg19.yaml.copy -n config.bianca.hg19.v1.yaml --reference-path /PROJECT_DATA:/proj/sens2022566/nobackup/patriksm/design_and_ref_files --reference-path PATH_TO_REPO:/proj/sens2022566/nobackup/patriksm/Twist_Solid_add-validation-ref-yaml/hydra-genetics
#  - singularity cache at 
# command: hydra-genetics prepare-environment container-path-update -c config.bianca.hg19.v1.yaml -n config.bianca.hg19.v2.yaml -p /proj/sens2022566/nobackup/patriksm/singularity_cache
```
<br />
*profile*

Edit bianca profile Twist_Solid_${TAG_OR_BRANCH}/Twist_Solid/profiles/bianca/config.yaml
```yaml
# Found at Twist_Solid_{TAG_OR_BRANCH}/snakemake-wrappers, use absolute_path with 'git+file:/'
wrapper-prefix="PATH_TO_WRAPPERS"
# ex: wrapper-prefix="git+file://proj/sens2022566/nobackup/patriksm/Twist_Solid_add-{TAG_OR_BRANCH}/snakemake-wrappers/"

# Update account info, change ADD_YOUR_ACCOUNT to your bianca project id
drmaa: " -A ADD_YOUR_ACCOUNT -N 1-1 -t {resources.time} -n {resources.threads} --mem={resources.mem_mb} --mem-per-cpu={resources.mem_per_cpu} --mem-per-cpu={resources.mem_per_cpu} --partition={resources.partition} -J {rule} -e slurm_out/{rule}_%j.err -o slurm_out/{rule}_%j.out"
```

#### Validate config files

```bash
# This will make sure that all design and reference files exists and haven't changed
# Warnings for possible file PATH/hydra-genetics and missing tbi files in config can be ignored
hydra-genetics --debug references validate -c config/config.yaml -c config/bianca/config.hg19.yaml -v config/references/design_files.hg19.yaml -v config/references/references.bianca.hg19.yaml -v config/references/nextseq.hg19.pon.yaml -v config/references/references.hg19.yaml  -p ${PATH_TO_design_and_ref_files} 
```

### Singularity

#### Setup env, snakemake-wrappers and Hydra-Genetics modules

```bash
# Extract tar.
singularity run --app copy-twist-solid-env twist_solid_add-validation-ref-yaml.sif Pipeline
singularity run --app copy-twist-pipeline twist_solid_add-validation-ref-yaml.sif Pipeline
singularity run --app copy-twist-solid-env twist_solid_add-validation-ref-yaml.sif Pipeline
singularity run --app copy-hydra-modules twist_solid_add-validation-ref-yaml.sif Pipeline
singularity run --app copy-snakemake-wrappers twist_solid_add-validation-ref-yaml.sif Pipeline

source Pipeline/twist_solid_venv/bin/activate

```

#### Create config and profile
```bash
# Create config folder with updated paths
singularity run --app create-config-folder-bianca twist_solid.sif /PATH/singularity_cache /PATH/design_and_ref_files /PATH/Pipeline

# Create profile with project id and path to snakemake-wrappers
singularity run --app create-profile-bianca twist_solid_add-validation-ref-yaml.sif sensXXXX /PATH/Pipeline/snakemake-wrappers

```

#### Validate config
```bash
singularity run --app validate-config-hg19-bianca twist_solid_add-validation-ref-yaml.sif  -c /PATH_TO_UPDATE/config/config.yaml -c /PATH_TO_UPDATE/config/bianca/config.hg19.yaml -p /proj/sensXXXX/nobackup/USER/design_and_ref_files
```

## Run Pipeline

```bash
# Create analysis
mkdir analysis
# Enter folder
cd analysis
# Copy config files
cp -r PATH_TO_UPDATED_CONFIGS/config .

# Create samples.tsv and units.tsv
```

### Manually

```bash
module load slurm-drmaa

# Enter runfolder
cd PATH_TO_ANALYSIS_FOLDER

# For conda
source /{PATH_TO_ENV}/venv/bin/activate
snakemake -s /{PATH_TO_PIPELINE}/Twist_Solid/workflow/Snakefile  --profile ${PATH_TO_UPDATED_PROFILE}/Twist_Solid/profiles/bianca

# Note that bianca may close your session before the workflow is done
```

### SBATCH

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