FROM continuumio/miniconda3:4.10.3
#FROM centos:7

ARG TAG_OR_BRANCH

RUN apt-get update &&  apt-get install -y build-essential

RUN . /root/.bashrc && conda init bash \
    && conda install -c conda-forge conda-pack  \
    && conda create --name twist_solid_${TAG_OR_BRANCH} python=3.9 -y \
    && conda activate twist_solid_${TAG_OR_BRANCH}


RUN mkdir -p /Twist_Solid/hydra-genetics
RUN echo ${TAG_OR_BRANCH} > /Twist_Solid/version.txt
RUN git clone https://github.com/snakemake/snakemake-wrappers.git
WORKDIR /Twist_Solid
RUN git clone https://github.com/snakemake/snakemake-wrappers.git
WORKDIR /Twist_Solid/hydra-genetics
RUN git config --global --add safe.directory /Twist_Solid/snakemake-wrappers
RUN git clone https://github.com/hydra-genetics/prealignment.git
RUN git config --global --add safe.directory /Twist_Solid/hydra-genetics/prealignment
RUN git clone https://github.com/hydra-genetics/alignment.git
RUN git config --global --add safe.directory /Twist_Solid/hydra-genetics/alignment
RUN git clone https://github.com/hydra-genetics/snv_indels.git
RUN git config --global --add safe.directory /Twist_Solid/hydra-genetics/snv_indels
RUN git clone https://github.com/hydra-genetics/annotation.git
RUN git config --global --add safe.directory /Twist_Solid/hydra-genetics/annotation
RUN git clone https://github.com/hydra-genetics/filtering.git
RUN git clone https://github.com/hydra-genetics/qc.git
RUN git config --global --add safe.directory /Twist_Solid/hydra-genetics/filtering
RUN git clone https://github.com/hydra-genetics/biomarker.git
RUN git config --global --add safe.directory /Twist_Solid/hydra-genetics/biomarker
RUN git clone https://github.com/hydra-genetics/fusions.git
RUN git config --global --add safe.directory /Twist_Solid/hydra-genetics/fusions
RUN git clone https://github.com/hydra-genetics/cnv_sv.git
RUN git config --global --add safe.directory /Twist_Solid/hydra-genetics/cnv_sv
RUN git clone https://github.com/hydra-genetics/compression.git
RUN git config --global --add safe.directory /Twist_Solid/hydra-genetics/compression
RUN git clone https://github.com/hydra-genetics/reports.git  
RUN git config --global --add safe.directory /Twist_Solid/hydra-genetics/reports
RUN git clone https://github.com/hydra-genetics/misc.git 
RUN git config --global --add safe.directory /Twist_Solid/hydra-genetics/misc

WORKDIR /Twist_Solid/
RUN git clone --branch ${TAG_OR_BRANCH} https://github.com/genomic-medicine-sweden/Twist_Solid.git
RUN git config --global --add safe.directory '/Twist_Solid/Twist_Solid'
RUN pip install -r Twist_Solid/requirements.txt

RUN . /root/.bashrc && conda init bash \ 
    && conda activate twist_solid_${TAG_OR_BRANCH} \
    && pip install -r Twist_Solid/requirements.txt \
    && pip install drmaa

RUN . /root/.bashrc && conda init bash \ 
    && conda activate twist_solid_${TAG_OR_BRANCH} \
    && conda pack -n twist_solid_${TAG_OR_BRANCH} -o /Twist_Solid/env.tar.gz \
    && chmod 555 /Twist_Solid/env.tar.gz
