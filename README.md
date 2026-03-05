<p align="center">
<a href="https://twist-solid.readthedocs.io">https://twist-solid.readthedocs.io</a>
</p>

# Twist_Solid
Pipeline for Solid tumours


## Code style validation
[![Lint](https://github.com/genomic-medicine-sweden/Twist_Solid/actions/workflows/lint.yaml/badge.svg)](https://github.com/genomic-medicine-sweden/Twist_Solid/actions/workflows/lint.yaml)
[![pycodestyle](https://github.com/genomic-medicine-sweden/Twist_Solid/actions/workflows/pycodestyle.yaml/badge.svg)](https://github.com/genomic-medicine-sweden/Twist_Solid/actions/workflows/pycodestyle.yaml)

## Code testing
[![pytest](https://github.com/genomic-medicine-sweden/Twist_Solid/actions/workflows/pytest.yaml/badge.svg)](https://github.com/genomic-medicine-sweden/Twist_Solid/actions/workflows/pytest.yaml)
[![Snakemake dry run](https://github.com/genomic-medicine-sweden/Twist_Solid/actions/workflows/snakemake-dry-run.yaml/badge.svg)](https://github.com/genomic-medicine-sweden/Twist_Solid/actions/workflows/snakemake-dry-run.yaml)

## Somalier Best Match Report
The pipeline includes a Somalier-based relatedness check that identifies the best genetic match for each sample. Specifically, it reports the best matching RNA sample for each DNA sample, and vice versa. This is used to verify sample identity and detect potential mixups. The report is delivered to `results/qc/somalier/sample_mixup_check_somalier.tsv`.
