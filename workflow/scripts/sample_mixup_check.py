
import os

vcf_dna_filenames = snakemake.input.id_snp_vcf_dna
vcf_rna_filenames = snakemake.input.id_snp_vcf_rna
report = open(snakemake.ouput.mixup_report, "w")

vcf_dict_dna = {}
dna_samples = {}
vcf_dict_rna = {}
rna_samples = {}


def read_vcf(vcf_filename, vcf_dict, samples):
    vcf = open(vcf_filename)
    sample_type = os.path.split(vcf_filename)[1].split(".")[0]
    vcf_dict[sample_type] = []
    samples[sample_type] = {}
    header = True
    for line in vcf:
        if header:
            if line[:6] == "#CHROM":
                header = False
            continue
        columns = line.strip().split("\t")
        chrom = columns[0][3:]
        pos = columns[1]
        key = chrom + "_" + pos
        FORMAT = columns[8].split(":")
        DATA = columns[9].split(":")
        GT_index = 0
        i = 0
        for f in FORMAT:
            if f == "GT":
                GT_index = i
            i += 1
        GT = DATA[GT_index]
        vcf_dict[sample_type].append(GT)
    vcf.close()


for vcf_filename in vcf_dna_filenames:
    read_vcf(vcf_filename, vcf_dict_dna, dna_samples)
for vcf_filename in vcf_rna_filenames:
    read_vcf(vcf_filename, vcf_dict_rna, rna_samples)

for rna_sample in rna_samples:
    for dna_sample in dna_samples:
        rna_samples[rna_sample][dna_sample] = 0
        i = 0
        for GT_rna in vcf_dict_rna[rna_sample]:
            GT_dna = vcf_dict_dna[dna_sample][i]
            if GT_rna == GT_dna:
                rna_samples[rna_sample][dna_sample] += 1
            i += 1

report.write("RNA_sample\tDNA_sample\tnr_matches\t%_match\n")
for rna_sample in rna_samples:
    best_dna_sample = ""
    best_gt_match = 0
    for dna_sample in rna_samples[rna_sample]:
        if rna_samples[rna_sample][dna_sample] > best_gt_match:
            best_dna_sample = dna_sample
            best_gt_match = rna_samples[rna_sample][dna_sample]
    report.write(f"{rna_sample}\t{best_dna_sample}\t{best_gt_match}\t{round(best_gt_match * 100 / 44.0, 1)}%\n")
