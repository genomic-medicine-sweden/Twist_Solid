# CNV calling merging and purity estimation
See the [cnv hydra-genetics module](https://snv_indels.readthedocs.io/en/latest/) documentation for more details on the softwares for cnv calling.

**Result files**

* `results/dna/cnv/{sample}_{type}/{sample}_{type}.pathology.cnv.html`
* `results/dna/cnv/{sample}_{type}/{sample}_{type}.pathology.cnv_report.tsv`
* `results/dna/cnv/{sample}_{type}/{sample}_{type}.purecn.cnv.html`
* `results/dna/cnv/{sample}_{type}/{sample}_{type}.purecn.cnv_report.tsv`
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

## CNV call file conversion to vcf
The CNVs call files from CNVkit and GATK CNV are converted to vcf format using **in-house scripts**. In both cases a vcf header is added followed by the following annotation for each field in the vcf file:

* CHROM - Chromosome
* POS - CNV start position
* ID - .
* QUAL - .
* FILTER - .
* REF - N
* ALT - <DEL\> / <DUP\> / <COPY_NORMAL\>
* INFO
    - SVTYPE - DEL / DUP / COPY_NORMAL
    - END - CNV end position
    - SVLEN - CNV length
    - LOG_ODDS_RATIO - CNV log odds copy ratio (0 = 2 copies)
    - CALLER - cnvkit / gatk
    - CN - Copy number (NA for cnvkit)
    - CORR_CN - Purity corrected copy number (calculated for GATK CNV)
    - PROBES - Number of probes within the cnv
    - BAF - SNP minor allele frequency
* FORMAT / DATA
    - GT - Genotype: 1/1 for homozygous deletion, 1/0 for heterozygous deletion and or duplication, and 0/0 for copy neutral
    - CN - Copy number (corrected for CNVkit and uncorrected for GATK CNV)
    - CCN - Corrected copy number (GATK CNV only)
    - CNQ - Number of probes in CNV
    - DP - Average coverage over region (cnvkit only)
    - BAF - SNP minor allele frequency
    - BAFQ - Number of SNPs for BAF (GATK CNV only)

**Options**

* sample_name: {sample}_{type}
* hom_del_limit: 0.5 - Used for setting GT to 1/1
* het_del_limit: 1.5 - Used for determining FORMAT:GT, ALT, and INFO:SVTYPE
* dup_limit: 2.5 - Used for determining FORMAT:GT, ALT, and INFO:SVTYPE

## CNV merging and annotation using SVDB
Merging the vcfs with CNV calls from the two callers and then annotating the calls with the frequency in a panel of normals is done by **[SVDB](https://github.com/J35P312/SVDB)** v2.6.0 using the following steps:

* **SVDB --merge**
    - --pass_only - Merge the two vcf files without actually merging overlapping regions
* **SVDB --query**
    - panel of normal - Use a SVDB panel of normal to annotate the frequency of overlapping regions  with at least 60 % overlap and at max 10000 bases from the breakpoints (see [references](dna_references.md) on how the PoN was created)
    - --out_frq Twist_AF - frequency annotation name in the INFO field
    - --out_occ Twist_OCC - occurrence annotation name in the INFO field

## CNV gene annotation
CNV regions that overlap with clinically relevant genes for amplifications (`cnv_amp_genes.bed`) and deletions (`cnv_loh_genes.bed`) are annotated separately and put into two different files by an **in-house script**. The relevant genes are annotated in the INFO field with the "Genes" tag.

**Amplification genes**

MTOR, NTRK1, MYCN, ALK, FGFR3, PDGFRA, KIT, FGFR4, EGFR, CDK6, MET, BRAF, FGFR1, MYC, CDKN2A, NTRK2, PTEN, FGFR2, KRAS, CDK4, NTRK3, ERBB2, AR

**Deletion genes**

BAP1, FAT1, CHD1, MCPH1, CDKN2A, CDKN2B, TSC1, PTEN, ATM, BRCA2, RB1, SPRED1, TSC2, PALB2, CDH1, FANCA, TP53, NF1, BRCA1, RAD51C

## CNV filtering
Filtering the CNV amplifications and deletions are performed by the [filtering hydra-genetics module](https://filtering.readthedocs.io/en/latest/).

### Amplification filtering
Genes and filtering criteria specified in `config_hard_filter_cnv_amp.yaml` are listed below:

* Filter cnvs with under 6 copies
* Filter cnvs found in more than 15% of the normal samples
* Filter cnvs not annotated in the INFO:Genes tag

### Deletion filtering
Genes and filtering criteria specified in `config_hard_filter_cnv_loh.yaml` are listed below. Observe that regions with neutral copies but with high BAF signal is not filtered as they are candidates for copy neutral loss of heterozygosity.

* Filter cnvs with over 1.4 copies if the BAF between 0.3 and 0.7 as well as all duplication (copy number > 2.5)
* Filter cnvs found in more than 15% of the normal samples
* Filter cnvs not annotated in the INFO:Genes tag

## CNV report
A combined report in tsv format is generated based on the filtered vcf files by an **in-house script**. The script also looks for 1p19q deletions and add that to the report using the following options:

**1p19q options**

* del_1p19q_cn_limit: 1.4 - Both the 1p and 19q chromosome must have a region of the arm with under 1.4 copies
* del_1p19q_chr_arm_fraction: 0.3 - The fraction of the arm that must be deleted

**Result files**

* `results/dna/cnv/{sample}_{type}/{sample}_{type}.pathology.cnv_report.tsv`
* `results/dna/cnv/{sample}_{type}/{sample}_{type}.purecn.cnv_report.tsv`

## CNV visulization
In addition to the tsv report a html report is also generated giving an visualization of the underlying data as well as cnv information. This report is generated by an **in-house script** using the following configuration in the `config.yaml`:

```yaml
cnv_vcf:
  - annotation: cnv_loh_genes
    filter: cnv_hard_filter_loh
  - annotation: cnv_amp_genes
    filter: cnv_hard_filter_amp
template: config/cnv_report_template.html
```

**Result files**

* `results/dna/cnv/{sample}_{type}/{sample}_{type}.pathology.cnv.html`
* `results/dna/cnv/{sample}_{type}/{sample}_{type}.purecn.cnv.html`

## Small CNV deletions
CNVkit and GATK CNV sometimes miss small deletions where only a number of exons are involved. A specialized small CNV caller in the form of an **in-house script** is therefore run on the genes mentioned in [CNV gene annotation](#cnv-gene-annotation). The caller uses the log-odds values calculated by GATK for each region, where a region approximately corresponds to one region in the design file, i.e. exon, intron (if present) or CNV-probe. For each gene region, the caller uses a sliding window to find the biggest drop and subsequent rise in copy number in the region. It then determines if the drop in copy number is big enough and is significantly lower than the surrounding region. If a large part or the entire region is deleted it will not be called but will instead be called by the other tools.

**Options**

* regions_file: cnv_deletion_genes.tsv - Genes and surrounding regions to call small CNVs in
* window_size: 4 - Sliding window size of data points which governs the minimum size that can be found
* region_max_size: 30 - Max size of the region in the form of number of data points
* min_log_odds_diff: 0.3 - Min difference in copy numbers between deletion and the rest of the
* min_nr_stdev_diff: 2.5 - Min number of standard deviations difference

**Result file**

* `results/dna/cnv/{sample}_{type}/{sample}_{type}.deletions.tsv`

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
The purity of the samples is estimated in the lab by a pathologist and sometimes the estimation is incorrect. To obtain an alternative estimation we use **[PureCN](https://github.com/lima1/PureCN)** v2.2.0.

### PureCN coverage
PureCN calculates the coverage in 25kb windows of the BWA-mem aligned and merged bam files and normalizes these against a panel of normal (see [references](dna_references.md) on how the PoN was created).

**Options**

* intervals: `targets_twist-gms-st_hg19_25000_intervals.txt` - panel of normal

### PureCN purity estimation
Based on the coverage and a vcf file from Mutect2 PureCN make its own segmentation with the additional help of a normal db created from a panel of normal (see [references](dna_references.md) on how the PoN was created). The germline SNPs allele frequency in combination with the called copy numbers is then used to search for the optimal purity and ploidity combination. The optimal values are determined by looking for the best fit between the germline SNPs allele frequency and integer copy numbers.

**Options**

* genome: "hg19" - The genome fasta reference version
* interval_padding: 100 - Padding size of the 25kb regions
* segmentation_method: "internal" - Use PureCNs own segmentation (CNVkit segmentation can for example also be used)
* fun_segmentation: "PSCBS" - Segmentation function
* --model betabin - Recoomended for more than 10-15 normal samples
* --mapping-bias-file `mapping_bias_nextseq_27_hg19.rds` - Panel of normal
* normaldb: `normalDB_nextseq_27_hg19.rds` - Panel of normal

## Manta
**Manta** v1.6.0 is used to call larger INDELs and other structural variant events. However the results are only reported and not used in any clinical anaylsis.

**Options**

* --exome - use exome mode as this is not WGS data
* --callRegions design_bed_file - only make call in the designed region

**Resources**

* threads: 4
* time: "8:00:00"
* mem_mb: 24576
* mem_per_cpu: 6144

**Result file**

* `results/dna/cnv/{sample}_{type}/{sample}_{type}.manta_tumorSV.vcf.gz`
