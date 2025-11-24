#Building version v0.24.0

Note! Before running the build script conda needs to be available. On HPC-clusters conda can usually be accessed through the 
module system, e.g. with ```module load miniconda3```. 

```bash
TAG_OR_BRANCH="v0.24.0" CONFIG_VERSION="v1.18.0" PIPELINE_NAME="Twist_Solid" PYTHON_VERSION="3.9" \
PIPELINE_GITHUB_REPO="https://github.com/genomic-medicine-sweden/Twist_Solid.git" \
CONFIG_GITHUB_REPO="https://github.com/clinical-genomics-uppsala/GMS560_config.git" \
CONFIG_NAME="GMS560_config" \
PATH_TO_apptainer_cache="path/to/.sif/on/remote" \
APPT_CACHE_STATUS="build" \
bash build_conda.sh config/references/design_files.hg19.yaml config/references/nextseq.hg19.pon.yaml \
config/references/references.hg19.yaml GMS560_config/config/reference_files.yaml
```


