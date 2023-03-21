# References, panel of normals and design files

## References overview
The following reference files, panel of normals and design files are need to run the Twist Solid Pipeline:

| Rule | Config name | File |
|-|-|-|
| reference | <div id="background_db"> background </div> | `background_panel_nextseq_noUmea_27_dp500_af015.tsv` |
| | <div id="artifact_db"> artifacts </div> | `artifact_panel_nextseq_36.tsv` |
| | <div id="reference_fasta">fasta</div> | `hg19.with.mt.fasta` |
| | fasta_rna | `GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa` |
| | dict | `hg19.with.mt.dict` |
| | fai | `hg19.with.mt.fai` |
| | <div id="design_bed">design_bed</div> | `pool1_pool2.sort.merged.padded20.cnv200.hg19.split_fusion_genes.reannotated.210608.bed` |
| | design_intervals | `pool1_pool2.sort.merged.padded20.cnv200.hg19.split_fusion_genes.MUC6_31_rm.exon_only.reannotated.210608.interval_list` |
| | design_intervals_gatk_cnv | `pool1_pool2_nochr_3c.sort.merged.padded20.cnv400.hg19.210311.met.annotated.bed.preprocessed.interval_list` |
| | <div id="design_bed_rna">design_bed_rna</div> | `Twist_RNA_Design5.annotated.bed` |
|_ _| design_intervals_rna | `Twist_RNA_Design5.annotated.interval_list` |
| <div id="arriba_reference">arriba</div> | assembly | `hg19.with.mt.fasta` |
| | <div id="arriba_blacklist">blacklist</div> | `arriba/arriba_v2.3.0/database/blacklist_hg19_hs37d5_GRCh37_v2.3.0.tsv.gz` |
| | <div id="arriba_gtf">gtf</div> | `hg19.refGene.gtf` |
| | <div id="arriba_extra_gff3_">extra</div> | `arriba/arriba_v2.3.0/database/protein_domains_hg19_hs37d5_GRCh37_v2.3.0.gff3` |
|_ _| <div id="arriba_extra_tsv">extra | `arriba/arriba_v2.3.0/database/known_fusions_hg19_hs37d5_GRCh37_v2.3.0.tsv.gz` |
| <div id="arriba_draw_fusion_cytobands">arriba_draw_fusion</div> | cytobands | `arriba/arriba_v2.3.0/database/cytobands_hg19_hs37d5_GRCh37_v2.3.0.tsv` |
| | <div id="arriba_draw_fusion_gtf">gtf</div> | `hg19/gtf/hg19.refGene.gtf` |
|_ _| <div id="arriba_draw_fusion_protein_domains">protein_domains</div> | `arriba/arriba_v2.3.0/database/protein_domains_hg19_hs37d5_GRCh37_v2.3.0.gff3` |
| <div id="cnv_amp_genes">annotate_cnv</div> | cnv_amp_genes | `cnv_amp_genes.bed` |
|_ _| <div id="cnv_loh_genes">cnv_loh_genes</div> | `cnv_loh_genes.bed` |
| bcftools_annotate | annotation_db | `small_exac_common_3.hg19.vcf.gz` |
| <div id="bcftools_filter">bcftools_filter_include_region</div> | exon | `pool1_pool2.sort.merged.padded20.cnv200.hg19.split_fusion_genes.MUC6_31_rm.exon_only.reannotated.210608.bed` |
| <div id="bcftools_filter_exclude_region">bcftools_filter_exclude_region</div> | blacklist | `cnvkit_germline_blacklist_20221221.bed` |
| <div id="bcftools_id_snps">bcftools_id_snps</div> | snps_bed | `ID_SNPs.bed` |
| <div id="bwa_me_ref">bwa_mem | amb </div> | `hg19.with.mt.amb` |
| | ann | `hg19.with.mt.ann` |
| | bwt | `hg19.with.mt.bwt` |
| | pac | `hg19.with.mt.pac` |
|_ _| sa | `hg19.with.mt.sa` |
| <div id="call_small_cnv_deletions">call_small_cnv_deletions</div> | regions_file | `cnv_deletion_genes.tsv` |
| cnvkit_batch | normal_reference | `cnvkit_nextseq_36.cnn` |
| <div id="normal_reference_hrd">cnvkit_batch_hrd</div> | normal_reference_hrd | `cnvkit_nextseq_27_HRD.cnn` |
| <div id="exon_skipping">exon_skipping</div> | design_bed | `Twist_RNA_Design5.annotated.bed` |
| <div id="fusioncatcher">fusioncatcher</div> | genome_path | `human_v102/` |
| gatk_collect_allelic_counts | SNP_interval | `gnomad_SNP_0.001_target.annotated.interval_list` |
| gatk_denoise_read_counts | normal_reference | `gatk_cnv_nextseq_36.hdf5` |
| gatk_get_pileup_summaries | sites | `small_exac_common_3.hg19.vcf.gz` |
|_ _| variants | `small_exac_common_3.hg19.vcf.gz` |
| <div id="genefuse_transcripts">gene_fuse | genes | `GMS560_fusion_w_pool2.hg19.221117.csv` |
|_ _| <div id="genefuse_fasta">fasta</div> | `hg19.with.mt.fasta` |
| <div id="hotspot_file">hotspot_annotation</div> | hotspots | `Hotspots_combined_regions_nodups.csv` |
| <div id="hotspot_report">hotspot_report</div> | hotspot_mutations | `Hotspots_combined_regions_nodups.csv` |
| <div id="manta_design_bed">manta_config_t</div> | extra | `pool1_pool2.sort.merged.padded20.cnv200.hg19.split_fusion_genes.210608.bed.gz` |
| <div id="msisensor_pro_pon">msisensor_pro</div> | PoN | `Msisensor_pro_reference_nextseq_36.list_baseline` |
| <div id="purecn_estimation_mapping_pon">purecn</div> | extra | `mapping_bias_nextseq_27_hg19.rds` |
| | <div id="purecn_estimation_normaldb">normaldb</div> | `normalDB_nextseq_27_hg19.rds` |
|_ _| <div id="purecn_estimation_intervals">intervals</div> | `targets_twist-gms-st_hg19_25000_intervals.txt` |
| <div id="purecn_coverage_intervals">purecn_coverage</div> | intervals | `targets_twist-gms-st_hg19_25000_intervals.txt` |
| <div id="report_fusions">report_fusions</div> | annotation_bed | `Twist_RNA_fusionpartners.bed` |
| <div id="genefuse_filter_fusions">report_gene_fuse</div> | filter_fusions | `filter_fusions_20221114.csv` |
| <div id="star_genome_index">star</div> | genome_index | `v2.7.10a_hg19/` |
|_ _| <div id="star_genome_extra">extra</div> | `hg19.refGene.gtf` |
| <div id="star_fusion">star_fusion</div> | genome_path | `GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/` |
| svdb_query | db_string | `all_TN_292_svdb_0.8_20220505.vcf` |
| <div id="vep_cache">vep</div> | vep_cache | `VEP/` |

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

**Result files**

Files needed by Twist Solid that are generated by the pipeline itself are listed below. These are available to download from the Uppsala Owncloud solution but should preferably be generated on in-house data.

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

### CNVkit

**units.tsv**

* column header: bam
* column data: path to merged bam files

**References**

* design bedfile
* fasta genome reference
* mappability file

### GATK CNV

**units.tsv**

* column header: bam
* column data: path to merged bam files

**References**

* design bedfile
* fasta reference genome
* fasta reference genome dictionary

### MSISensor-pro

**units.tsv**

* column header: bam
* column data: path to merged bam files

**References**

* design bedfile
* fasta reference genome

**Options**

-c 50 - minimal coverage, recommended for WES: 20; WGS: 15

### SVDB

**units.tsv**

* column header: cnv_vcf
* column data: svdb merged cnv vcf files

**Options**

--overlap - Overlap used to cluster variants (default 0.8), Use "extra" parameter to set this in config

### Artifacts

**units.tsv**

* column header: vcf
* column data: unfiltered merged vcf files

### Background

**units.tsv**

* column header: vcf
* column data: unfiltered merged vcf files

**Options**

* min_dp - Min read depth to be included (default: 500)
* max_af - Max allele frequency to be included (default: 0.015)

### PureCN
**OBS!** The best way to run PureCN is still to be determined. At present PureCN uses Mutect2 filtered vcf files (not hard filtered). This is not the same as the other PoNs that use ensembled vcf files.

#### Target interval file
Target interval file for hg19 with 25000 in target bin size also including of target regions.
```bash
singularity docker://hydragenetics/purecn:2.2.0 Rscript $PURECN/IntervalFile.R --fasta ${fasta_ref} --in-file ${design_bed} --out-file ${intervals_file} --export ${optimized_bed} --genome hg19 --average-off-target-width 25000 --off-target
```

#### PureCN Panel of normal

**units.tsv**

* column header: bam
* column data: path to merged bam files
* column header: vcf
* column data: unfiltered GATK vcf files

**Reference**

* Target interval file

## Pipeline specific files
These are design files and other pipeline specific only available to download from the Uppsala Owncloud solution.

* Design files
    - `pool1_pool2.sort.merged.padded20.cnv200.hg19.split_fusion_genes.reannotated.210608.bed` - design bed
    - `pool1_pool2.sort.merged.padded20.cnv200.hg19.split_fusion_genes.MUC6_31_rm.exon_only.reannotated.210608.bed` - design bed file containing only exons
    - `pool1_pool2.sort.merged.padded20.cnv200.hg19.split_fusion_genes.MUC6_31_rm.exon_only.reannotated.210608.interval_list` - design interval file containing only exons
    - `pool1_pool2_nochr_3c.sort.merged.padded20.cnv400.hg19.210311.met.annotated.bed.preprocessed.interval_list` - design interval file used in GATK CNV PoN creation (must be identical)
    - `Twist_RNA_Design5.annotated.bed` - RNA design bed
    - `Twist_RNA_Design5.annotated.interval_list` - RNA design interval file
    - `ID_SNPs.bed` - List of RNA ID SNPs
* Hotspots (annotation of vcfs and mutation qc report)
    - `Hotspots_combined_regions_nodups.csv` - Positions, transcript information, etc on clinically relevant regions
* GeneFuse
    - `GMS560_fusion_w_pool2.hg19.221117.csv` - Genes and its exonic positions included in fusion calling
    - `filter_fusions_20221114.csv` - Filtering criteria for false positive prone fusion partners
* CNVkit
    - `cnvkit_germline_blacklist_20221221.bed` - list of regions excluded from the germline vcf file
* GATK CNV
    - `gnomad_SNP_0.001_target.annotated.interval_list` - Bed file with CNV backbone SNPs which are selected from GnomAD with over 0.1% global population frequency
* Small CNV deletions
    - `cnv_deletion_genes.tsv` - File defining gene and its surrounding regions used for small CNV deletion. Same deletion genes as in the CNV deletion reports
* Report RNA fusions
    - `Twist_RNA_fusionpartners.bed` - Bed file used for annotation of fusion partner exons
