# CNV calling merging and purity estimation
See the [cnv hydra-genetics module](https://snv_indels.readthedocs.io/en/latest/) documentation for more details on the softwares for cnv calling.

**Result files**

* `results/dna/cnv/{sample}_{type}/{sample}_{type}.purecn.cnv.html`
* `results/dna/cnv/{sample}_{type}/{sample}_{type}.purecn.cnv_report.tsv`
* `results/dna/cnv/{sample}_{type}/{sample}_{type}.pathology.cnv.html`
* `results/dna/cnv/{sample}_{type}/{sample}_{type}.pathology.cnv_report.tsv`
* `results/dna/cnv/{sample}_{type}/{sample}_{type}.deletions.tsv`
* `results/dna/cnv/{sample}_{type}/{sample}_{type}.manta_tumorSV.vcf.gz`

## CNV calling
CNVs are called using both CNVkit and GATK CNV, followed by merging and annotation by SVDB and finally filtered and visualized into an html report. Additionally, smaller CNV deletions are called in specific genes using an in-house script.

### CNVkit
CNV segmentation is performed by **[CNVkit](https://cnvkit.readthedocs.io/en/stable/)** v0.9.9 on BWA-mem aligned and merged bam-files. CNVkit uses a panel of normal (see [references](dna_references.md) on how the PoN was created) and a germline filtered vcf file. To call the final CNVs the program uses the estimated tumor purity, which can be from pathologists estimates or from purecn estimates.  
The following steps are in included in the CNVkit calling:

* **CNVkit batch** - Make segmentation
    - method: "hybrid"
    - Panel of normal
* **CNVkit call** - Call CNVs
    - Germline vcf with variants called in the sample (see [Germline vcf](#germline-vcf) for details)
    - Tumor purity
* **CNVkit diagram** - Pdf with chromosome overview
* **CNVkit scatter** - Scatterplot of segmentation and BAF

### GATK CNV
CNV segmentation is performed by **[GATK CNV](https://gatk.broadinstitute.org/hc/en-us/articles/360035535892-Somatic-copy-number-variant-discovery-CNVs-)** v4.1.9.0 on BWA-mem aligned and merged bam-files. GATK CNV uses a panel of normal (see [references](dna_references.md) on how the PoN was created) and a precompiled germline vcf file. GATK do not use the estimated tumor purity so the copy number levels are instead adjusted when the segmentation is coverted to vcf file.
The following steps are in included in the GATK CNV calling:

* **gatk CollectReadCounts** - Coverage
    - Design bed file in picards interval file format
* **gatk CollectAllelicCounts** - Allele frequencies of SNPs
    - Fasta reference genome
    - SNP file with SNPs found in GnomAD with population frequency above 0.1%.
* **gatk DenoiseReadCounts** - Normalize coverage against panel of normal
    - Uses a panel of normal
* **gatk ModelSegments** - Make segmentation
* **gatk CallCopyRatioSegments** - Call CNVs

## CNV seqmentation file conversion to vcf


## CNV merging and annotation using SVDB

## CNV filtering

## CNV report

## CNV visulization

## Germline vcf
The germline vcf used by CNVkit and the CNV html report is based on the [VEP annotated vcf](dna_snv_indels.md#vep) file from the SNV and INDEL calling. Annotated vcfs are hard filtered first by removing black listed regions with noisy germline VAFs in normal samples and then filtered by a number of filtering criteria described below. See the [filtering hydra-genetics module](https://filtering.readthedocs.io/en/latest/) for additional information.

### Exclude exonic regions
Use **[bcftools filter -T](https://samtools.github.io/bcftools/bcftools.html)** v1.15 to exclude variants overlapping blacklisted regions defined in a bed file.

**References**

* Bed file with blacklisted regions

### Filter vcf
The germline vcf file are filtered using the **[hydra-genetics filtering](https://filtering.readthedocs.io/en/latest/)** functionality included in v0.15.0. The filters are specified in the config file `config_hard_filter_germline.yaml` and consists of the following filters:

* germline - Filter germline variants when the GnomAD global population allele frequency is below 0.1%
* allele observations - Filter variants with fewer than 50 supporting reads in both alleles
* variant type - Filter variants that are not SNVs

## Purity estimation using PureCN
**[PureCN](https://github.com/lima1/PureCN)**

## Manta
