# Fusions in RNA
See the [fusions hydra-genetics module](https://snv_indels.readthedocs.io/en/latest/) documentation for more details on the softwares for fusion calling.

Fusion calling is performed using three different fusion callers; [Arriba](https://github.com/suhrig/arriba), [Star-Fusion](https://github.com/STAR-Fusion/STAR-Fusion) and [fusioncatcher](https://github.com/ndaniel/fusioncatcher). Both Arriba and Star-Fusion uses the [Star](https://github.com/alexdobin/STAR) for alignment but with different settings while fusioncatcher uses its own aligner. After fusion calling the fusions are filtered depending on software and then merged into a fusion report.

**Result files**

* `results/rna/fusion/{sample}_{type}.fusion_report.tsv`
* `results/rna/fusion/{sample}_{type}.arriba.fusions.tsv`
* `results/rna/fusion/{sample}_{type}.arriba.fusions.pdf`
* `results/rna/fusion/{sample}_{type}.fusioncatcher.fusion_predictions.txt`
* `results/rna/fusion/{sample}_{type}.star-fusion.fusion_predictions.tsv`

## Arriba

### Alignment with STAR
Merged fastq files are aligned with [Star](https://github.com/alexdobin/STAR) v2.7.10a before fusion calling.

**References**

Star genome index (see [references](references.md#star-genome-index))

**Options** (Recommended options by Arriba)

* --quantMode GeneCounts
* --sjdbGTFfile `hg19.refGene.gtf` (see [references](references.md#arriba-230))
* --outSAMtype BAM SortedByCoordinate
* --chimSegmentMin 10
* --chimOutType WithinBAM SoftClip
* --chimJunctionOverhangMin 10
* --chimScoreMin 1
* --chimScoreDropMax 30
* --chimScoreJunctionNonGTAG 0
* --chimScoreSeparation 1
* --alignSJstitchMismatchNmax 5 -1 5 5
* --chimSegmentReadGapMax 3


### Fusion calling with Arriba
Star aligned bam-files are used for fusion calling with [Arriba](https://github.com/suhrig/arriba) v2.3.0.

**References**

* assembly: fasta reference genome
* blacklist: `blacklist_hg19_hs37d5_GRCh37_v2.3.0.tsv.gz` - (see [references](references.md#arriba-230))
* gtf: `hg19.refGene.gtf` - (see [references](references.md#arriba-230))
* -p `protein_domains_hg19_hs37d5_GRCh37_v2.3.0.gff3` - (see [references](references.md#arriba-230))
* -k `known_fusions_hg19_hs37d5_GRCh37_v2.3.0.tsv.gz` - (see [references](references.md#arriba-230))

**Result file**

* `results/rna/fusion/{sample}_{type}.arriba.fusions.tsv`

### Fusion images
Arriba produces a pdf file containing a figure for every fusion called with a schematic presentation of the exons involved, breakpoints, coverage and directions of the fusion partners in the fusion.

**References**

cytobands: `cytobands_hg19_hs37d5_GRCh37_v2.3.0.tsv` - (see [references](references.md#arriba-230))
gtf: `hg19.refGene.gtf` - (see [references](references.md#arriba-230))
protein_domains: `protein_domains_hg19_hs37d5_GRCh37_v2.3.0.gff3` - (see [references](references.md#arriba-230))

**Result file**

* `results/rna/fusion/{sample}_{type}.arriba.fusions.pdf`

## Star-Fusion
[Star-Fusion](https://github.com/STAR-Fusion/STAR-Fusion) v1.10.1

**Result file**

* `results/rna/fusion/{sample}_{type}.star-fusion.fusion_predictions.tsv`

## Fusioncatcher
[fusioncatcher](https://github.com/ndaniel/fusioncatcher) v1.33 and genome v102

**Result file**

* `results/rna/fusion/{sample}_{type}.fusioncatcher.fusion_predictions.txt`

## Fusion report

**Result file**

* `results/rna/fusion/{sample}_{type}.fusion_report.tsv`
