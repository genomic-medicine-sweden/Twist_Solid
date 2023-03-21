# Prealignment
See the [prealignment hydra-genetics module](https://prealignment.readthedocs.io/en/latest/) documentation for more details on the softwares. Default hydra-genetics settings/resources are used if no configuration is specfied.

<br />
![dag plot](images/prealignment.png)

## Pipeline output files:
Only temporary intermediate files created.

## Trimming
Trimming of fastq files is performed by **[fastp](https://github.com/OpenGene/fastp)** v0.20.1.  

### Configuration


**Resources**

| **Options** | **Value** |
|-------------|-|
| mem_mb | 30720 |
| mem_per_cpu | 6144 |
| threads | 5 |

## Merging
Merging of fastq files belonging to the same sample are performed by simply concatenating the files with **cat**.

<br />