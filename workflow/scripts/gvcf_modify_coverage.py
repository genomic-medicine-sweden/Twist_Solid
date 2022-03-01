
import gzip

in_coverage = gzip.open(snakemake.input.coverage, "rt")
in_gvcf = gzip.open(snakemake.input.gvcf, "rt")
out_gvcf = gzip.open(snakemake.output.gvcf, "wt")

coverage_data = in_coverage.read().split("\n")
gvcf_data = in_gvcf.read().split("\n")

coverage_data_i = -1
coverage_data_chrom = "chr1"
coverage_data_endpos = 0
coverage_data_cov = 0
header = True
for line in gvcf_data:
    if header:
        out_gvcf.write(line + "\n")
        if line[:6] == "#CHROM" :
            header = False
        continue
    columns = line.split("\t")
    if len(columns) < 10 :
        continue
    chrom = columns[0]
    pos = int(columns[1])
    format = columns[8].split(":")
    data = columns[9].split(":")
    DP_i = 0
    for f in format:
        if f == "DP":
            break
        DP_i += 1
    while not(chrom == coverage_data_chrom and pos <= coverage_data_endpos):
        coverage_data_i += 1
        coverage_data_chrom = coverage_data[coverage_data_i].split("\t")[0]
        coverage_data_endpos = int(coverage_data[coverage_data_i].split("\t")[2])
        coverage_data_cov = coverage_data[coverage_data_i].split("\t")[3]
    data[DP_i] = str(coverage_data_cov)
    out_gvcf.write(columns[0])
    i = 1
    for column in columns[1:]:
        if i == 9:
            out_gvcf.write("\t" + data[0])
            for d in data[1:]:
                out_gvcf.write(":" + d)
        else:
            out_gvcf.write("\t" + column)
        i += 1
    out_gvcf.write("\n")
out_gvcf.close()
