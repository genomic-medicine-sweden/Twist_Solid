FROM python:3.8-slim

ARG TAG_OR_BRANCH=develop
#ARG SINGULARITY_VERSION=3.8.0


RUN apt-get update \
 #&& apt-get install -y --no-install-recommends gosu build-essential libssl-dev uuid-dev libgpgme11-dev squashfs-tools libseccomp-dev wget pkg-config git cryptsetup libslurm-dev slurm-client \
 && apt-get install -y --no-install-recommends git libslurm-dev slurm-client\
 && apt-get purge -y --auto-remove && rm -rf /var/lib/apt/lists/*
 
#RUN git clone https://github.com/natefoo/slurm-drmaa.git /slurm-drmaa
#WORKDIR /slurm-drmaa
#RUN ./configure && make && make install

#RUN export VERSION=1.13 OS=linux ARCH=amd64 && wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz \
#    && tar -C /usr/local -xzvf go$VERSION.$OS-$ARCH.tar.gz && rm go$VERSION.$OS-$ARCH.tar.gz

#RUN echo 'export GOPATH=${HOME}/go' >> ~/.bashrc && echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc && bash ~/.bashrc

#ENV GOPATH=${HOME}/go
#ENV PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin

#RUN wget https://github.com/sylabs/singularity/releases/download/v${SINGULARITY_VERSION}/singularity-ce-${SINGULARITY_VERSION}.tar.gz
#RUN tar -xzf singularity-ce-${SINGULARITY_VERSION}.tar.gz
#WORKDIR singularity-ce-${SINGULARITY_VERSION}
#RUN ./mconfig  --without-suid  && make -C builddir && make -C builddir install
#RUN cd ..  && rm -r singularity-ce-${SINGULARITY_VERSION}

RUN mkdir -p /Twist_Solid/hydra-genetics
RUN echo ${TAG_OR_BRANCH} > /Twist_Solid/version.txt
WORKDIR /Twist_Solid/hydra-genetics
RUN git clone https://github.com/snakemake/snakemake-wrappers.git
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
#RUN hydra-genetics references download -o design_and_ref_files -v Twist_Solid/config/references/design_files.hg19.yaml -v Twist_Solid/config/references/nextseq.hg19.pon.yaml -v Twist_Solid/config/references/references.hg19.yaml

#COPY ./docker/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN sed 's/#hydra_local_path: "PATH_TO_REPO"/hydra_local_path: \/Twist_Solid/' -i  /Twist_Solid/Twist_Solid/config/config.yaml
#RUN sed 's/resources: "/resources: "\/Twist_Solid\/Twist_Solid\//' -i  /Twist_Solid/Twist_Solid/config/config.yaml
#RUN sed 's/singularity_schema: "/singularity_schema: "\/Twist_Solid\/Twist_Solid\//' -i  /Twist_Solid/Twist_Solid/config/config.yaml


#RUN chmod +x /usr/local/bin/entrypoint.sh
#ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
