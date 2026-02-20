# Twist Solid

## GMS560 gene panel
The GMS560 gene panel is a combined gDNA target capture panel and a cDNA junction capture panel for comprehensive genomic profiling of solid tumors, designed and validated as a collaborative effort within the Genomic Medicine Sweden Solid (GMS) tumor work package. The design allows for clinical use as well as for research applications.

## Twist RNA panel
The RNA panel is a smaller panel designed to capture fusion transcripts and specific exon skipping events (e.g., MET exon 14). It complements the DNA panel by providing RNA-level evidence for fusions and splice variants.

## FFPE and ctDNA pipeline
The Twist Solid pipeline analyses both tumor-only FFPE and ctDNA samples and reports small variant, CNVs, fusions and biomarkers (TMB, MSI, HRD) as well as a QC report and a sample report.
The RNA analysis pipeline focuses on fusion calling using multiple callers (Arriba, Star-Fusion, FusionCatcher) and exon skipping detection.
The pipeline is build using the [hydra-genetics](https://github.com/hydra-genetics/hydra-genetics) framework. All new releases are run through real data sets and compared to previous versions to ensure that the expected result are generated.
