
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
        if line[:6] == "#CHROM" :
            out_gvcf.write("##FORMAT=<ID=DP_mosdepth,Number=1,Type=Integer,Description=\"Read depth reported by mosdepth\">\n")
            header = False
        out_gvcf.write(line + "\n")
        continue
    columns = line.split("\t")
    if len(columns) < 10 :
        continue
    chrom = columns[0]
    pos = int(columns[1])
    format = "%s:DP_mosdepth" % (columns[8])
    while not(chrom == coverage_data_chrom and pos <= coverage_data_endpos):
        coverage_data_i += 1
        coverage_data_chrom = coverage_data[coverage_data_i].split("\t")[0]
        coverage_data_endpos = int(coverage_data[coverage_data_i].split("\t")[2])
        coverage_data_cov = coverage_data[coverage_data_i].split("\t")[3]
    data = "%s:%s" % (columns[9], str(coverage_data_cov))
    out_gvcf.write(columns[0])
    for column in columns[1:]:
        out_gvcf.write("\t" + column)
    out_gvcf.write("\n")
out_gvcf.close()
