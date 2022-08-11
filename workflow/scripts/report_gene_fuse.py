
fusions = open(snakemake.input.fusions)
report = open(snakemake.output.report, "w")
min_unique_reads = int(open(snakemake.params.min_unique_reads))

report.write("Gene1\tGene2\tUnique_reads\n")

FP_gene_pairs = ["NPM1_ALK", "CLTC_NTRK3", "MASH2_ALK"]
Noisy_genes = {"RSPO2": 8, "ABL1": 8, "BRAF": 8, "EZR": 8}

for line in fusions:
    lline = line.strip().split("\t")
    if line[:8] != "#Fusion:":
        continue
    gene1 = line.split(" ")[1].split("_")[0]
    gene2 = line.split("___")[1].split("_")[0]
    unique_reads = int(line.split("unique:")[1].split(")")[0])
    min_u_r = min_unique_reads
    if gene1 in Noisy_gene_pairs:
        min_u_r = Noisy_gene_pairs[gene1]
    if gene2 in Noisy_gene_pairs:
        min_u_r = Noisy_gene_pairs[gene2]
    if unique_reads < min_u_r:
        continue
    key = gene1 + "_" + gene2
    if key in FP_gene_pairs:
        continue
    report.write(f"{gene1}\t{gene2}\t{unique_reads}\n")
report.close()
