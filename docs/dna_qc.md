# QC
See the [qc hydra-genetics module](https://snv_indels.readthedocs.io/en/latest/) documentation for more details on the softwares for the quality control.

**Result file**

* `results/dna/qc/multiqc_DNA.html`

## MultiQC
A MultiQC html report is generated using **MultiQC** v1.11. The report starts with a general statistics table showing the most important QC-values followed by additional QC data and diagrams. The qc data is generated using FastQC, samtools, picard, and GATK.

**Options**

* `config/multiqc_config_dna.yaml` - Config of the general statistics table
* `config/config.yaml` - Configuration of input files to MultiQC in the config file
```yaml
multiqc:
  container: "docker://hydragenetics/multiqc:1.11"
  reports:
    DNA:
      config: "config/multiqc_config_dna.yaml"
      included_unit_types: ["N", "T"]
      qc_files:
        - "qc/fastqc/{sample}_{type}_{flowcell}_{lane}_{barcode}_fastq1_fastqc.zip"
        - "qc/fastqc/{sample}_{type}_{flowcell}_{lane}_{barcode}_fastq2_fastqc.zip"
        - "qc/picard_collect_alignment_summary_metrics/{sample}_{type}.alignment_summary_metrics.txt"
        - "qc/picard_collect_duplication_metrics/{sample}_{type}.duplication_metrics.txt"
        - "qc/picard_collect_hs_metrics/{sample}_{type}.HsMetrics.txt"
        - "qc/picard_collect_insert_size_metrics/{sample}_{type}.insert_size_metrics.txt"
        - "qc/samtools_stats/{sample}_{type}.samtools-stats.txt"
        - "qc/gatk_calculate_contamination/{sample}_{type}.contamination.table"
```

## FastQC
**FastQC** v0.11.9 is run on the raw fastq-files.

## Samtools
**Samtools stats** v1.15 is run on BWA-mem aligned and merged bam files.

## Picard
**Picard** v2.25.0 is run on BWA-mem aligned and merged bam files collecting a number of metrics. The metrics calculated are listed below:

* **picard CollectAlignmentSummaryMetrics** - using a fasta reference genome file
* **picard CollectDuplicationMetrics**
* **picard CollectHsMetrics** - using a fasta reference genome file, a design bed file, and with the option COVERAGE_CAP=5000
* **picard CollectInsertSizeMetrics**

## GATK
Cross-sample contamination is estimated using **GATK** v4.1.9.0. The contamination is calculated by **gatk GetPileupSummaries** making pileups for input SNP positions on BWA-mem aligned and merged bam files followed by evaluating these pileups with **gatk CalculateContamination**. Contamination levels should very low so that already at 1% there is reason to be concerned.
