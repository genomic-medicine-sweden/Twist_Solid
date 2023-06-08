#!/usr/bin/env bash
set -e

eval "$(conda shell.bash hook)"

TWIST_SOLID_VERSION="${TAG_OR_BRANCH:-develop}"

conda create --name twist_solid_${TWIST_SOLID_VERSION} python=3.9 -y

conda activate twist_solid_${TWIST_SOLID_VERSION}

conda install -c conda-forge pip -y

if [ -d Twist_Solid_${TWIST_SOLID_VERSION} ];
then
    rm -fr Twist_Solid_${TWIST_SOLID_VERSION}
fi

mkdir Twist_Solid_${TWIST_SOLID_VERSION}

git clone --branch ${TWIST_SOLID_VERSION} https://github.com/genomic-medicine-sweden/Twist_Solid.git Twist_Solid_${TWIST_SOLID_VERSION}/Twist_Solid

pip install -r Twist_Solid_${TWIST_SOLID_VERSION}/Twist_Solid/requirements.txt 

conda deactivate

conda pack -n twist_solid_${TWIST_SOLID_VERSION} -o Twist_Solid_${TWIST_SOLID_VERSION}/env.tar.gz

mkdir -p Twist_Solid_${TWIST_SOLID_VERSION}/hydra-genetics

git clone https://github.com/snakemake/snakemake-wrappers.git Twist_Solid_${TWIST_SOLID_VERSION}/snakemake-wrappers

git clone https://github.com/hydra-genetics/prealignment.git Twist_Solid_${TWIST_SOLID_VERSION}/hydra-genetics/prealignment
git clone https://github.com/hydra-genetics/alignment.git Twist_Solid_${TWIST_SOLID_VERSION}/hydra-genetics/alignment
git clone https://github.com/hydra-genetics/snv_indels.git Twist_Solid_${TWIST_SOLID_VERSION}/hydra-genetics/snv_indels
git clone https://github.com/hydra-genetics/annotation.git Twist_Solid_${TWIST_SOLID_VERSION}/hydra-genetics/annotation
git clone https://github.com/hydra-genetics/filtering.git Twist_Solid_${TWIST_SOLID_VERSION}/hydra-genetics/filtering
git clone https://github.com/hydra-genetics/qc.git Twist_Solid_${TWIST_SOLID_VERSION}/hydra-genetics/qc
git clone https://github.com/hydra-genetics/biomarker.git Twist_Solid_${TWIST_SOLID_VERSION}/hydra-genetics/biomarker
git clone https://github.com/hydra-genetics/fusions.git Twist_Solid_${TWIST_SOLID_VERSION}/hydra-genetics/fusions
git clone https://github.com/hydra-genetics/cnv_sv.git Twist_Solid_${TWIST_SOLID_VERSION}/hydra-genetics/cnv_sv
git clone https://github.com/hydra-genetics/compression.git Twist_Solid_${TWIST_SOLID_VERSION}/hydra-genetics/compression
git clone https://github.com/hydra-genetics/misc.git Twist_Solid_${TWIST_SOLID_VERSION}/hydra-genetics/misc
git clone https://github.com/hydra-genetics/reports.git Twist_Solid_${TWIST_SOLID_VERSION}/hydra-genetics/reports

tar -zcvf Twist_Solid_${TWIST_SOLID_VERSION}.tar.gz Twist_Solid_${TWIST_SOLID_VERSION}

if [ -d Twist_Solid_${TWIST_SOLID_VERSION} ];
then
    rm -fr Twist_Solid_${TWIST_SOLID_VERSION}
fi
