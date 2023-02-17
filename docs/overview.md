# Overview of the pipeline
Here is a brief overview of the entire pipeline. For details see subsections and the (hydra-genetics)[https://github.com/hydra-genetics/hydra-genetics] documentation.

## DNA
1. **Input files**: fastq
2. **Trimming** using fastp
3. **Alignment** using BWA-mem
4. **Mark duplicates** using Picard
5. **SNV and INDEL**  
  5.1 Calling using Mutect2 and Vardict  
  5.2 Annotation using VEP and hydra-genetics  
  5.3 Filtering using bcftools and hydra-genetics  
6. **CNV**  
  6.1 Calling using CNVkit and GATK CNV  
  6.2 Merging using SVDB  
  6.3 Annotation using SVDB and hydra-genetics  
  6.4 Filtering using hydra-genetics  
  6.5 CNV html report using hydra-genetics
7. **Fusion** calling using GeneFuse
8. **Biomarkers**  
  8.1 TMB using hydra-genetics  
  8.2 MSI score using MSIsensor-Pro  
  8.3 HRD using CNVkit and ScarHRD  
9. **QC**  
  9.1 QC measures from Samtools, Picard, FastQC, GATK  
  9.2 MultiQC hmtl report  
  9.3 Hotspot coverage report  

## RNA
1. **Input files**: fastq
2. **Alignment** using Star
3. **Fusions**  
  3.1 Fusion calling using Arriba, StarFusion, FusionCatcher  
  3.2 Filtering and report using in-house script  
  3.3 Fusion images using Arriba  
4. **Exon skipping** using in-house script
5. **ID-SNP** calling using bcftools
6. **QC**  
  6.1 QC measures from Samtools, Picard, FastQC, Mosdepth  
  6.2 MultiQC hmtl report  
  6.3 House keeping gene coverage  
