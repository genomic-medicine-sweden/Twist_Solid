# Fusion calling in DNA
See the [fusions hydra-genetics module](https://hydra-genetics-fusions.readthedocs.io/en/latest/) documentation for more details on the softwares for fusion calling. Default hydra-genetics settings/resources are used if no configuration is specified.

<br />
![dag plot](images/fusions.png){: style="height:25%;width:25%"}

## Pipeline output files:

* `results/dna/fusion/{sample}_{type}.gene_fuse_report.tsv`
* `results/dna/fusion/{sample}_{type}.fuseq_wes.report.csv`

## Fusions calling using GeneFuse
DNA fusion calling is performed by **[GeneFuse](https://github.com/OpenGene/GeneFuse)** v0.6.1 on fastq-files. It uses a gene transcript target file to limit the number of targets to analyze.

### Configuration

**References**

* [Fasta reference](references.md#genefuse_fasta) genome
* [Gene transcript](references.md#genefuse_transcripts) file with genomic positions for all exons include in the analysis  

<br />
**Cluster resources**

| **Options** | **Value** |
|-------------|-|
| mem_mb | 36864 |
| mem_per_cpu | 6144 |
| threads | 6 |
| time | "8:00:00" |

## GeneFuse Filtering and report
The output from GeneFuse is filtered and then reported into a fusion report using the in-house script [report_gene_fuse.py](https://github.com/genomic-medicine-sweden/Twist_Solid/blob/develop/workflow/scripts/report_gene_fuse.py) ([rule and config](softwares.md#report_gene_fuse)). The following filter criteria is used:

* Fusions must have at least 6 unique supporting reads.
* Very noisy fusion pairs found in almost all samples (defined in [`filter_fusions_20230214.csv`](references.md#genefuse_filter_fusions)) are removed:
    - NPM1::ALK
    - CLTC::NTRK3
    - MSH2_ALK
    - MSH2_HIP1
* Noisy fusion pairs found in some samples (defined in [`filter_fusions_20230214.csv`](references.md#genefuse_filter_fusions)) are filtered individually on the number of uniquely supporting reads:
    - LMNA::EZR 9
    - ABL1::STRN 7
    - EZR::ALK 8
    - RSPO2::BRAF 8
    - LMNA::HIP1 12
    - NPM1::BICC1 11
    - RSPO2::ERG 13

### Result file

* `results/dna/fusion/{sample}_{type}.gene_fuse_report.tsv`

<br />

## Fusions calling using FuSeq_WES
DNA fusion calling is performed by **[FuSeq_WES](https://github.com/nghiavtr/FuSeq_WES)** v1.0.1 on bam-files. It uses a gene transcript target file to limit the number of targets to analyze.

### Configuration
**References**

* [Transcript json](references.md#fuseq_wes_json)
* [Transcript sqlite](references.md#fuseq_wes_sqlite)
* [Fusion database](references.md#fuseq_wes_fusion_db)
* [Paralog database](references.md#fuseq_wes_paralog_db)

<br />

**Software settings**

| **Options** | **Value** | **Description** |
|-------------|-|-|
| params | [`cnv_amplfication_genes.tsv`](references.md#call_small_cnv_amplifications) | Genes and surrounding regions to call small CNVs in |

**Cluster resources**

| **Options** | **Value** |
|-------------|-|
| mem_mb | 12288 |
| mem_per_cpu | 6144 |
| threads | 2 |
| time | "24:00:00" |

## FuSeq_WES Filtering and report
The output from FuSeq_WES is filtered and then reported into a fusion report using the in-house script [filter_report_fuseq_wes.py](https://github.com/genomic-medicine-sweden/Twist_Solid/blob/develop/workflow/scripts/filter_report_fuseq_wes.py) ([rule and config](softwares.md#report_fuseq_wes)). The following filter criteria is used:

### Configuration
**Software settings**

| **Options** | **Value** | **Description** |
|-------------|-|-|
| filter_on_fusiondb | True | Only keep fusions found in the fusion database |
| gene_white_list | [`fuseq_wes_gene_white_list.txt`](references.md#fuseq_wes_white_list) | Only keep fusions with at least one gene in the gene white list |
| gtf | [`hg19.refGene.gtf`](references.md#filter_report_fuseq_wes) | Transcript annotation |
| min_support | 30 | Minimal total number of supporting reads |
| transcript_black_list | [`fuseq_wes_transcript_black_list.txt`](references.md#fuseq_wes_transcript_black_list) | Transcripts that should not be used in annotation |

### Result file

* `results/dna/fusion/{sample}_{type}.fuseq_wes_report.tsv`

<br />
