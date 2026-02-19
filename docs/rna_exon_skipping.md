# Exon skipping
Exon skipping cannot be called by the fusion callers as they only search for fusion in-between genes at a certain distance. To be able to call the clinically relevant MET exon 14 skipping and EGFRvIII the in-house script [exon_skipping.py](https://github.com/genomic-medicine-sweden/Twist_Solid/blob/develop/workflow/scripts/exon_skipping.py) ([rule and config](softwares.md#house_keeping_gene_coverage)) is used that uses the junction file from the Star aligner output. The junction file contains reads split information that in combination with exon start and end positions can be used to find reads skipping exons. The current analysis is limited to the MET and EGFR genes and is reported in a tsv file. To be reported the number of skipping reads must be at least 10% of the total number of reads in that junction. 

## Configuration
**References**

* [Annotated RNA design bed](references.md#exon_skipping)

## Result file

* `results/rna/{sample}_{type}/fusion/{sample}_{type}.exon_skipping.tsv`

<br />

## CTAT-Splicing
[CTAT-Splicing](https://github.com/TrinityCTAT/CTAT-SPLICING) is used to find cancer relevant introns with abnormal splicing, like MET exon 14 skipping and EGFR exon deletions. It uses the junction file from the Star aligner output and a reference genome library.

### Configuration
**References**
* [CTAT genome library](references.md#star_fusion)

**Software settings**
* [ctat_splicing_call](softwares.md#ctat_splicing_call)
* [ctat_splicing_filter](softwares.md#ctat_splicing_filter)

### Result file
* `fusions/ctat_splicing_call/{sample}_{type}.cancer.introns`
* `fusions/ctat_splicing_call/{sample}_{type}.ctat-splicing.igv.html`
* `fusions/ctat_splicing_filter/{sample}_{type}.cancer.introns.filtered.tsv`
