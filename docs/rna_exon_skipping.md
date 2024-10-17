# Exon skipping
Exon skipping cannot be called by the fusion callers as they only search for fusion in-between genes at a certain distance. To be able to call the clinically relevant MET exon 14 skipping and EGFRvIII the in-house script [exon_skipping.py](https://github.com/genomic-medicine-sweden/Twist_Solid/blob/develop/workflow/scripts/exon_skipping.py) ([rule and config](softwares.md#house_keeping_gene_coverage)) is used that uses the junction file from the Star aligner output. The junction file contains reads split information that in combination with exon start and end positions can be used to find reads skipping exons. The current analysis is limited to the MET and EGFR genes and is reported in a tsv file. To be reported the number of skipping reads must be at least 10% of the total number of reads in that junction. 

## Configuration
**References**

* [Annotated RNA design bed](references.md#exon_skipping)

## Result file

* `results/rna/{sample}_{type}/fusion/{sample}_{type}.exon_skipping.tsv`

<br />
