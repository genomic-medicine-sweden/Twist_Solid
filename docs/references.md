# References, panel of normals and design files

## References overview
The following reference files, panel of normals and design files are need to run the Twist Solid Pipeline:

| Rule | Config name | File |
|-|-|-|
| reference | background | `background_panel_nextseq_noUmea_27_dp500_af015.tsv` |
| | artifacts | `artifact_panel_nextseq_36.tsv` |
| | fasta | `hg19.with.mt.fasta` |
| | fasta_rna | `GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa` |
| | dict | `hg19.with.mt.dict` |
| | fai | `hg19.with.mt.fai` |
| | design_bed | `pool1_pool2.sort.merged.padded20.cnv200.hg19.split_fusion_genes.reannotated.210608.bed` |
| | design_intervals | `pool1_pool2.sort.merged.padded20.cnv200.hg19.split_fusion_genes.MUC6_31_rm.exon_only.reannotated.210608.interval_list` |
| | design_intervals_gatk_cnv | `pool1_pool2_nochr_3c.sort.merged.padded20.cnv400.hg19.210311.met.annotated.bed.preprocessed.interval_list` |
| | design_bed_rna | `Twist_RNA_Design5.annotated.bed` |
|_ _| design_intervals_rna | `Twist_RNA_Design5.annotated.interval_list` |
| arriba | assembly | `hg19.with.mt.fasta` |
| | blacklist | `arriba/arriba_v2.3.0/database/blacklist_hg19_hs37d5_GRCh37_v2.3.0.tsv.gz` |
| | gtf | `hg19.refGene.gtf` |
| | extra | `arriba/arriba_v2.3.0/database/protein_domains_hg19_hs37d5_GRCh37_v2.3.0.gff3` |
|_ _| extra | `arriba/arriba_v2.3.0/database/known_fusions_hg19_hs37d5_GRCh37_v2.3.0.tsv.gz` |
| arriba_draw_fusion | cytobands | `arriba/arriba_v2.3.0/database/cytobands_hg19_hs37d5_GRCh37_v2.3.0.tsv` |
| | gtf | `hg19/gtf/hg19.refGene.gtf` |
|_ _| protein_domains | `arriba/arriba_v2.3.0/database/protein_domains_hg19_hs37d5_GRCh37_v2.3.0.gff3` |
| annotate_cnv | cnv_amp_genes | `cnv_amp_genes.bed` |
|_ _| cnv_loh_genes | `cnv_loh_genes.bed` |
| bcftools_annotate | annotation_db | `small_exac_common_3.hg19.vcf.gz` |
| bcftools_filter_include_region | exon | `pool1_pool2.sort.merged.padded20.cnv200.hg19.split_fusion_genes.MUC6_31_rm.exon_only.reannotated.210608.bed` |
| bcftools_filter_exclude_region | blacklist | `cnvkit_germline_blacklist_20221221.bed` |
| bcftools_id_snps | snps_bed | `ID_SNPs.bed` |
| bwa_mem | amb | `hg19.with.mt.amb` |
| | ann | `hg19.with.mt.ann` |
| | bwt | `hg19.with.mt.bwt` |
| | pac | `hg19.with.mt.pac` |
|_ _| sa | `hg19.with.mt.sa` |
| call_small_cnv_deletions | regions_file | `cnv_deletion_genes.tsv` |
| cnvkit_batch | normal_reference | `cnvkit_nextseq_36.cnn` |
| cnvkit_batch_hrd | normal_reference_hrd | `cnvkit_nextseq_27_HRD.cnn` |
| exon_skipping | design_bed | `Twist_RNA_Design5.annotated.bed` |
| fusioncatcher | genome_path | `human_v102/` |
| gatk_collect_allelic_counts | SNP_interval | `gnomad_SNP_0.001_target.annotated.interval_list` |
| gatk_denoise_read_counts | normal_reference | `gatk_cnv_nextseq_36.hdf5` |
| gatk_get_pileup_summaries | sites | `small_exac_common_3.hg19.vcf.gz` |
|_ _| variants | `small_exac_common_3.hg19.vcf.gz` |
| gene_fuse | genes | `GMS560_fusion_w_pool2.hg19.221117.csv` |
|_ _| fasta | `hg19.with.mt.fasta` |
| hotspot_annotation | hotspots | `Hotspots_combined_regions_nodups.csv` |
| hotspot_report | hotspot_mutations | `Hotspots_combined_regions_nodups.csv` |
| manta_config_t | extra | `pool1_pool2.sort.merged.padded20.cnv200.hg19.split_fusion_genes.210608.bed.gz` |
| msisensor_pro | PoN | `Msisensor_pro_reference_nextseq_36.list_baseline` |
| purecn | extra | `mapping_bias_nextseq_27_hg19.rds` |
| | normaldb | `normalDB_nextseq_27_hg19.rds` |
|_ _| intervals | `targets_twist-gms-st_hg19_25000_intervals.txt` |
| purecn_coverage | intervals | `targets_twist-gms-st_hg19_25000_intervals.txt` |
| report_fusions | annotation_bed | `Twist_RNA_fusionpartners.bed` |
| report_gene_fuse | filter_fusions | `filter_fusions_20221114.csv` |
| star | genome_index | `v2.7.10a_hg19/` |
|_ _| extra | `hg19.refGene.gtf` |
| star_fusion | genome_path | `GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/` |
| svdb_query | db_string | `all_TN_292_svdb_0.8_20220505.vcf` |
| vep | vep_cache | `VEP/` |

## Downloadable reference files

### Fasta reference genome
hg19 with mitochondria but without HLA and without decoys
```bash
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/*'
gunzip *.fa.gz
cat *.fa > hg19.with.mt.fasta
```

### BWA indexes
```bash
bwa index hg19.with.mt.fasta
```

### VEP v105
```bash
wget http://ftp.ensembl.org/pub/grch37/release-105/variation/vep/homo_sapiens_refseq_vep_105_GRCh37.tar.gz
```

### GNOMAD common SNPs
Used by GATK pileup (CNV calling) and GATK contamination
```bash
gsutil cp gs://gatk-best-practices/somatic-b37/small_exac_common_3.vcf* .
```
Adapt file to hg19 (add chr at all lines and in header)

### Mappability file
Used in CNVkit PoN creation
```bash
wget https://raw.githubusercontent.com/etal/cnvkit/master/data/access-5k-mappable.hg19.bed
```

### Arriba v2.3.0
```bash
singularity exec -B /path/to/references:/references docker://hydragenetics/arriba:2.3.0 download_references.sh hs37d5+RefSeq
wget https://github.com/suhrig/arriba/releases/download/v2.3.0/arriba_v2.3.0.tar.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.refGene.gtf.gz
```

### Star genome index
```bash
singularity exec docker://hydragenetics/star:2.7.10a STAR --runThreadN 8 --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles Human_genome.fasta
```

### Fusioncather v102
```bash
wget http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.aa
wget http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.ab
wget http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.ac
wget http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.ad
cat human_v102.tar.gz.* | tar xz
ln -s human_v102 current
```

### Star-fusion
```bash
wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play.tar.gz
```

## Panel of normals

Files needed by Twist Solid that are generated by the reference pipeline are listed below. Many of the panel of normals are available to download from the Uppsala Owncloud solution but should preferably be generated on in-house data.  
<br />
**Note**  
The input needed to create these files depend on the result files of the main pipeline itself, meaning that the samples need to be analysed with dummy panel of normals first.  
<br />
**Result files**

* `cnvkit.PoN.cnn`
* `gatk_cnv_panel_of_normal.hdf5`
* `Msisensor_pro_reference.list_baseline`
* `background_panel.tsv`
* `artifact_panel.tsv`
* `svdb_cnv.vcf`
* `normalDB_hg19.rds`
* `mapping_bias_nextseq_27_hg19.rds`

### Create samples and units
Files and samples used in the generation of the panel of normals are specified in samples.tsv and units.tsv. Required files and header names are listed down below.
Adapt out file specification (`workflow/rules/common_references.smk`) and comment out files that should not be generated.

### Run command
Run the pipeline in the same way as the [standard pipeline](running.md), using a reference specific profile and Snakefile:
```bash
snakemake --profile profiles/uppsala_ref/ -s workflow/Snakefile_references.smk
```
<br />
**Note**  
The `units.tsv` file needs to be adapted depending which panel of normals are created (see below) and should contain all the samples needed to create the panel of normals.

### CNVkit

**units.tsv**

| Header | Data | Description |
|-|-|-|
| bam | `bam_dna/{sample}_{type}.bam` | Merged bam files created as output of the Twist Solid pipeline from normal FFPE samples |

**Reference files**

* design bedfile
* fasta genome reference
* mappability file

### GATK CNV

**units.tsv**

| Header | Data | Description |
|-|-|-|
| bam | `bam_dna/{sample}_{type}.bam` | Merged bam files created as output of the Twist Solid pipeline from normal FFPE samples |

**Reference files**

* design bedfile
* fasta reference genome
* fasta reference genome dictionary

### MSISensor-pro

**units.tsv**

| Header | Data | Description |
|-|-|-|
| bam | `bam_dna/{sample}_{type}.bam` | Merged bam files created as output of the Twist Solid pipeline from normal FFPE samples |

**Reference files**

* design bedfile
* fasta reference genome

**Software settings**

| Options | Value | Description |
|-|-|-|
| extra | -c 50 | minimal coverage, recommended for WES: 20; WGS: 15 |

### SVDB

**units.tsv**

| Header | Data | Description |
|-|-|-|
| cnv_vcf | `results/dna/additional_files/cnv/{sample}_{type}/{sample}_{type}.pathology.svdb_query.vcf` | SVDB merged CNV vcf files created as output of the <br />Twist Solid pipeline from both normal and <br />tumor FFPE samples |

**Software settings**

| Options | Value | Description |
|-|-|-|
| extra | --overlap 0.8 | Overlap used to cluster variants (default 0.8) |

### Artifacts

**units.tsv**

| Header | Data | Description |
|-|-|-|
| vcf | `results/dna/additional_files/vcf/{sample}_{type}.annotated.vcf.gz` | Unfiltered and merged vcf files created as output of the Twist Solid pipeline <br />from normal FFPE samples |

### Background

**units.tsv**

| Header | Data | Description |
|-|-|-|
| vcf | `results/dna/additional_files/vcf/{sample}_{type}.annotated.vcf.gz` | Unfiltered and merged vcf files created as output of the Twist Solid pipeline <br />from normal FFPE samples |

**Software settings**

| Options | Value | Description |
|-|-|-|
| min_dp | 500 | Min read depth to be included (default: 500) |
| max_af | 0.015 | Max allele frequency to be included (default: 0.015) |

### PureCN
**OBS!** The best way to run PureCN is still to be determined. At present PureCN uses Mutect2 filtered vcf files (not hard filtered). This is not the same as the other PoNs that use ensembled vcf files.

#### Target interval file
Target interval file for hg19 with 25000 in target bin size also including of target regions.
```bash
singularity docker://hydragenetics/purecn:2.2.0 Rscript $PURECN/IntervalFile.R --fasta ${fasta_ref} --in-file ${design_bed} --out-file ${intervals_file} --export ${optimized_bed} --genome hg19 --average-off-target-width 25000 --off-target
```

#### PureCN Panel of normal

**units.tsv**

| Header | Data | Description |
|-|-|-|
| bam | `bam_dna/{sample}_{type}.bam` | Merged bam files created as output of the Twist Solid pipeline from <br />normal FFPE samples |
| vcf | `results/dna/additional_files/vcf/gatk_mutect2_{sample}_{type}.vcf.gz` | Mutect2 softfiltered vcf files created as output of the Twist Solid pipeline <br />from normal FFPE samples |

**Software settings**

| Options | Value | Description |
|-|-|-|
| intervals | Target interval file | File created by the command described above |

## Pipeline specific files
These are design files and other pipeline specific only available to download from the Uppsala Owncloud solution.

| File type | File | Description |
|-|-|-|
| Design files | `pool1_pool2.sort.merged.padded20.cnv200.hg19.split_fusion_genes`<br />`.reannotated.210608.bed` | Design bed |
| | `pool1_pool2.sort.merged.padded20.cnv200.hg19.split_fusion_genes`<br />`.MUC6_31_rm.exon_only.reannotated.210608.bed` | Design bed file containing only exons |
| | `pool1_pool2.sort.merged.padded20.cnv200.hg19.split_fusion_genes`<br />`.MUC6_31_rm.exon_only.reannotated.210608.interval_list` | Design interval file containing only exons |
| | `pool1_pool2_nochr_3c.sort.merged.padded20.cnv400.hg19.210311.met`<br />`.annotated.bed.preprocessed.interval_list` | Design interval file used in GATK CNV PoN creation |
| | `Twist_RNA_Design5.annotated.bed` | RNA design bed |
| | `Twist_RNA_Design5.annotated.interval_list` | RNA design interval file |
|_ _| `ID_SNPs.bed` | List of RNA ID SNPs |
| Hotspots | `Hotspots_combined_regions_nodups.csv` | Positions, transcript information, etc on clinically relevant regions |
| GeneFuse | `GMS560_fusion_w_pool2.hg19.221117.csv` | Genes and its exonic positions included in fusion calling |
| | `filter_fusions_20221114.csv` | Filtering criteria for false positive prone fusion partners |
| CNVkit | `cnvkit_germline_blacklist_20221221.bed` | List of regions excluded from the germline vcf file |
| GATK CNV | `gnomad_SNP_0.001_target.annotated.interval_list` | Bed file with CNV backbone SNPs which are selected from <br />GnomAD with over 0.1% global population frequency |
| Small CNV deletions | `cnv_deletion_genes.tsv` | File defining gene and its surrounding regions used for <br />small CNV deletion. Same deletion genes as in the <br />CNV deletion reports |
| Report RNA fusions | `Twist_RNA_fusionpartners.bed` | Bed file used for annotation of fusion partner exons |
