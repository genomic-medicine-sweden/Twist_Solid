# Exon skipping
Exon skipping cannot be called by the fusion callers as they only search for fusion in-between genes at a certain distance. To be able to call the clinically relevant MET exon 14 skipping and EGFRvIII an **in-house script** is used that uses the junction file from the Star aligner output. The junction file contains reads split information that in combination with exon start and end positions can be used to find reads skipping exons. The current analysis is limited to the MET and EGFR genes and is reported in a tsv file.

**References**

* annotated RNA design bed

**Result files**

* `results/rna/fusion/{sample}_{type}.exon_skipping.tsv`
