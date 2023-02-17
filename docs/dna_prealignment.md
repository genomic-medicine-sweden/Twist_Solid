# Prealignment
See the [prealignment hydra-genetics module](https://prealignment.readthedocs.io/en/latest/) documentation for more details on the softwares.

## Trimming
Trimming of fastq files is performed by **fastp** v0.20.1.  

**Resources**  

* threads: 5  
* mem_mb: 30720  
* mem_per_cpu: 6144  

## Merging
Merging of fastq files belonging to the same sample are performed by simply concatenating the files with **cat**.
