# Biomarkers
See the [biomarkers hydra-genetics module](https://snv_indels.readthedocs.io/en/latest/) documentation for more details on the softwares for the respective biomarkers.

**Result files**

* `results/dna/tmb/{sample}_{type}.TMB.txt`
* `results/dna/msi/{sample}_{type}.msisensor_pro.score.tsv`
* `results/dna/hrd/{sample}_{type}.purecn.scarhrd_cnvkit_score.txt`
* `results/dna/hrd/{sample}_{type}.pathology.scarhrd_cnvkit_score.txt`

## Tumor mutational burden (TMB)
TMB is a measure of the frequency of somatic mutations and is usually measured as mutations per megabase. The size of design of the exons is approximately 1.55Mb. However, by validating the TMB for GMS560 against Foundation One and TSO500 TMB the effective design size is adjusted to 0.98Mb. The TMB is calculated using an **in-house script** which counts the number of nsSNVs and divide by the adjusted design size. Variants must fulfill the following criteria to be counted:

* filter_nr_observations: 1 - Max seen once in panel of normal samples
* dp_limit: 100 - Minimum read depth of 100
* vd_limit: 10 - Minimum 10 observations of variant allele
* af_lower_limit: 0.05 - Minimum 5% allele frequency
* af_upper_limit: 0.45 - Maximum 45% allele frequency
* gnomad_limit: 0.0001 - Germline filter of 0.01% population frequency
* db1000g_limit: 0.0001 - Germline filter of 0.01% population frequency
* background_sd_limit: 5 - At least 5 standard deviation above background
* nssnv_tmb_correction: 1.02 - Variant size times correction factor (correction factor = 1 / adjusted design size)

The main result is the TMB calculated using nsSNV only. However, TMB calculated using both nsSNVs and sSNVs are also provided as well as all the variants passing all filters.

**Result file**

* `results/dna/tmb/{sample}_{type}.TMB.txt`

## Microsatellite instability (MSI)
To determine MSS or MSI status of the samples the percentage of sites that have microsatellite instability are calculated using **[MSIsensor-pro]([https://github.com/xjtu-omics/msisensor-pro])** v1.1.a. When the more than 10% of the sites are instable the sample is determined to have MSI status. The program uses a panel of normal to determine the normal level of instability in the used sites.

**Reference**

* Panel of normal for MSIsensor-pro (see [references](dna_references.md) on how the PoN was created)

**Result file**

* `results/dna/msi/{sample}_{type}.msisensor_pro.score.tsv`

## Homologous recombination deficiency (HRD) - in development
**OBS! The Homologous recombination deficiency score is still under development**  
A homologous recombination deficiency score is calculated using **scarhrd** v20200825 using cnvkit segmentation files as input. The cnvkit panel of normal for HRD is created from a design file where the extra CNV-probes were removed as coverage in these regions tended to be more affected in low quality samples. The segmentation is sensitive to the estimated purity. Therefore, a score based on both the pathology and purecn estimated tumor content is reported. The cutoff for HRD is still to be determined but is somewhere around 50 which is slightly higher than the Myriad HRD score cutoff of 42.

**Reference for cnvkit**

* Panel of normals created by cnvkit with extra CNV-probes removed (see [references](dna_references.md) on how the PoN was created)

**Options**

* reference_name: "grch37"
* seqz: FALSE

**Result files**

* `results/dna/hrd/{sample}_{type}.purecn.scarhrd_cnvkit_score.txt`
* `results/dna/hrd/{sample}_{type}.pathology.scarhrd_cnvkit_score.txt`
