# ID SNPs
The RNA design includes a number of probes covering SNPs that can be used to check if check that the RNA-sample is the same as the DNA-sample and thereby avoid sample swaps. The ID SNPs are called using **[bcftools mpileup](https://samtools.github.io/bcftools/bcftools.html#mpileup)** and **[bcftools call](https://samtools.github.io/bcftools/bcftools.html#call)**v1.15 resulting in a vcf file.

## Configuration
**References**

* ID SNP bed file: [`ID_SNPs.bed`](references.md#bcftools_id_snps)
* [RNA fasta](references.md#star_fusion) reference genome from Star-Fusion: `GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa` - (see [references](references.md#star-fusion))

<br />
**Software settings for bcftools mpileup (hard coded)**

| **Option** | **Description** |
|-------------|-|
| -O u | Output uncompressed BCF file (piped to bcftools call) |
| -d 1000000 | max read depth to consider |

<br />
**Software settings for bcftools call (hard coded)**

| **Option** | **Description** |
|-------------|-|
| --skip-variants indels | only look at SNPs |
| -m | multiallelic caller |
| -O v | Output uncompressed vcf file (piped to bcftools call) |

## Result file

* `results/rna/{sample}_{type}/id_snps/{sample}_{type}.id_snps.vcf`

<br />
