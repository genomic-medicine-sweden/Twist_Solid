# Fusion calling in DNA
See the [fusions hydra-genetics module](https://snv_indels.readthedocs.io/en/latest/) documentation for more details on the softwares for fusion calling.

**Result file**

* `results/dna/fusion/{sample}_{type}.gene_fuse_report.tsv`

## Fusions calling using GeneFuse
DNA fusion calling is performed by **[GeneFuse](https://github.com/OpenGene/GeneFuse)** v0.6.1 on fastq-files. It uses a gene transcript target file to limit the number of targets to analyse.

**References**

* Fasta reference genome
* Gene transcript file with genomic positions for all exons include in the analysis

**Resources**

* threads: 6
* time: "8:00:00"
* mem_mb: 36864
* mem_per_cpu: 6144

## Filtering and report
The output from GeneFuse is filtered and then reported into a fusion report using the in-house script [report_gene_fuse.py](https://github.com/genomic-medicine-sweden/Twist_Solid/blob/develop/workflow/scripts/report_gene_fuse.py) ([rule](https://github.com/genomic-medicine-sweden/Twist_Solid/blob/develop/workflow/rules/report_gene_fuse.smk)). The following filter criteria is used:

* Fusions must have at least 6 unique supporting reads.
* Very noisy fusion pairs found in almost all samples (defined in `filter_fusions_20230214.csv`) are removed:
    - NPM1::ALK
    - CLTC::NTRK3
    - MSH2_ALK
    - MSH2_HIP1
* Noisy fusion pairs found in some samples (defined in `filter_fusions_20230214.csv`) are filtered individually on the number of uniquely supporting reads:
    - LMNA::EZR 9
    - ABL1::STRN 7
    - EZR::ALK 8
    - RSPO2::BRAF 8
    - LMNA::HIP1 12
    - NPM1::BICC1 11
    - RSPO2::ERG 13

**Result file**

* `results/dna/fusion/{sample}_{type}.gene_fuse_report.tsv`
