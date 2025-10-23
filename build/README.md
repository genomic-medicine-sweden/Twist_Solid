#Building version v0.24.0

```bash
TAG_OR_BRANCH="v0.24.0" CONFIG_VERSION="v1.15.0" PIPELINE_NAME="Twist_Solid" PYTHON_VERSION="3.9" \
PIPELINE_GITHUB_REPO="https://github.com/genomic-medicine-sweden/Twist_Solid.git" \
CONFIG_GITHUB_REPO="https://github.com/clinical-genomics-uppsala/GMS560_config.git" \
PATH_TO_apptainer_cache="path/to/.sif/on/remote" \
bash build_conda.sh config/references/design_files.hg19.yaml config/references/nextseq.hg19.pon.yaml \
config/references/references.hg19.yaml GMS560_config/config/reference_files.yaml
```


