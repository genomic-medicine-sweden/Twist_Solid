# References, panel of normals and design files

## References

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
```
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

* `results/cnvkit.PoN.cnn`
* `results/gatk_cnv_panel_of_normal.hdf5`
* `results/Msisensor_pro_reference.list_baseline`
* `results/background_panel.tsv`
* `results/artifact_panel.tsv`
* `results/svdb_cnv.vcf`
* `results/normalDB_hg19.rds`
* `results/mapping_bias_nextseq_27_hg19.rds`

### Create samples and units
Files and samples used in the generation of the panel of normals are specified in samples.tsv and units.tsv. For needed files and header names see below.
Adapt out file specification (`workflow/rules/common_references.smk`) and comment out files that should not be generated.

### Run command
Run the pipeline in the same way as the [standard pipeline](running.dm) but with a slightly modified run command:
```bash
snakemake --profile profiles/uppsala_ref/ -s workflow/Snakefile_references.smk
```

### CNVkit

**units.tsv**

column header: bam
column data: path to merged bam files

**References**

* design bedfile
* fasta genome reference
* mappability file

### GATK CNV

**units.tsv**

column header: bam
column data: path to merged bam files

**References**

* design bedfile
* fasta reference genome
* fasta reference genome dictionary

### MSISensor-pro

**units.tsv**

column header: bam
column data: path to merged bam files

**References**

* design bedfile
* fasta reference genome

**Options**

-c 50 - minimal coverage, recommended for WES: 20; WGS: 15

### SVDB

**units.tsv**

column header: cnv_vcf
column data: svdb merged cnv vcf files

**Options**

--overlap - Overlap used to cluster variants (default 0.8), Use "extra" parameter to set this in config

### Artifacts

**units.tsv**

column header: vcf
column data: unfiltered merged vcf files

### Background

**units.tsv**

column header: vcf
column data: unfiltered merged vcf files

**Options**

min_dp - Min read depth to be included (default: 500)
max_af - Max allele frequency to be included (default: 0.015)

## PureCN
**OBS!** PureCN uses Mutect2 filtered vcf files (not hard filtered). This is not the same as the other PoNs that use ensembled vcf files.

### Target interval file
Target interval file for hg19 with 25000 in target bin size also including of target regions.
```bash
singularity docker://hydragenetics/purecn:2.2.0 Rscript $PURECN/IntervalFile.R --fasta {fasta_ref} --in-file {design_bed} --out-file {intervals_file} --export {optimized_bed} --genome hg19 --average-off-target-width 25000 --off-target
```

### PureCN Panel of normal

**units.tsv**

column header: bam
column data: path to merged bam files
column header: vcf
column data: unfiltered GATK vcf files

**Reference**

Target interval file

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