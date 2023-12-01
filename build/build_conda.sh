#!/usr/bin/env bash
set -e

eval "$(conda shell.bash hook)"

TAG_OR_BRANCH="${TAG_OR_BRANCH:-develop}"

conda create --name twist_solid_${TAG_OR_BRANCH} python=3.9 -y

conda activate twist_solid_${TAG_OR_BRANCH}

conda install -c conda-forge pip -y

if [ -d Twist_Solid_${TAG_OR_BRANCH} ];
then
    rm -fr Twist_Solid_${TAG_OR_BRANCH}
fi

mkdir Twist_Solid_${TAG_OR_BRANCH}

git clone --branch ${TAG_OR_BRANCH} https://github.com/genomic-medicine-sweden/Twist_Solid.git Twist_Solid_${TAG_OR_BRANCH}/Twist_Solid

pip install -r Twist_Solid_${TAG_OR_BRANCH}/Twist_Solid/requirements.txt 

conda deactivate

conda pack -n twist_solid_${TAG_OR_BRANCH} -o Twist_Solid_${TAG_OR_BRANCH}/env.tar.gz

mkdir -p Twist_Solid_${TAG_OR_BRANCH}/hydra-genetics

git clone https://github.com/snakemake/snakemake-wrappers.git Twist_Solid_${TAG_OR_BRANCH}/snakemake-wrappers

git clone https://github.com/hydra-genetics/prealignment.git Twist_Solid_${TAG_OR_BRANCH}/hydra-genetics/prealignment
git clone https://github.com/hydra-genetics/alignment.git Twist_Solid_${TAG_OR_BRANCH}/hydra-genetics/alignment
git clone https://github.com/hydra-genetics/snv_indels.git Twist_Solid_${TAG_OR_BRANCH}/hydra-genetics/snv_indels
git clone https://github.com/hydra-genetics/annotation.git Twist_Solid_${TAG_OR_BRANCH}/hydra-genetics/annotation
git clone https://github.com/hydra-genetics/filtering.git Twist_Solid_${TAG_OR_BRANCH}/hydra-genetics/filtering
git clone https://github.com/hydra-genetics/qc.git Twist_Solid_${TAG_OR_BRANCH}/hydra-genetics/qc
git clone https://github.com/hydra-genetics/biomarker.git Twist_Solid_${TAG_OR_BRANCH}/hydra-genetics/biomarker
git clone https://github.com/hydra-genetics/fusions.git Twist_Solid_${TAG_OR_BRANCH}/hydra-genetics/fusions
git clone https://github.com/hydra-genetics/cnv_sv.git Twist_Solid_${TAG_OR_BRANCH}/hydra-genetics/cnv_sv
git clone https://github.com/hydra-genetics/compression.git Twist_Solid_${TAG_OR_BRANCH}/hydra-genetics/compression
git clone https://github.com/hydra-genetics/misc.git Twist_Solid_${TAG_OR_BRANCH}/hydra-genetics/misc
git clone https://github.com/hydra-genetics/reports.git Twist_Solid_${TAG_OR_BRANCH}/hydra-genetics/reports

tar -zcvf Twist_Solid_${TAG_OR_BRANCH}.tar.gz Twist_Solid_${TAG_OR_BRANCH}

if [ -d Twist_Solid_${TAG_OR_BRANCH} ];
then
    rm -fr Twist_Solid_${TAG_OR_BRANCH}
fi
