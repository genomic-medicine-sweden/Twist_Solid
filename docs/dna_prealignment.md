# Prealignment
See the [prealignment hydra-genetics module](https://prealignment.readthedocs.io/en/latest/) documentation for more details on the softwares.

## Trimming
Trimming of fastq files is performed by **[fastp](https://github.com/OpenGene/fastp)** v0.20.1.  

**Resources**

| **Options** | **Value** |
|-------------|-|
| mem_mb | 30720 |
| mem_per_cpu | 6144 |
| threads | 5 |

## Merging
Merging of fastq files belonging to the same sample are performed by simply concatenating the files with **cat**.
