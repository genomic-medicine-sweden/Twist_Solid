
import subprocess
import time

genes = ["GAPDH", "GUSB", "OAZ1", "POLR2A"]

bam_file = snakemake.input.bam
bedfilename = snakemake.input.bed
outfile = open(snakemake.output.result, "w")
outfile.write("Gene\tAvg_coverage\n")


for gene in genes:
    regions = []
    bedfile = open(bedfilename)
    for line in bedfile:
        if line.find(gene) != -1:
            lline = line.strip().split("\t")
            regions.append([lline[0], lline[1], lline[2]])
    bedfile.close()
    coverage_sum = 0
    coverage_nr_pos = 0
    coverage_list = []
    for region in regions:
        region = region[0] + ":" + region[1] + "-" + region[2]
        samtools_coverage = subprocess.check_output(f"samtools depth -d 5000000 -a -r {region} {bam_file}",
                                                    shell=True).decode('ascii')
        for line in samtools_coverage.split("\n"):
            if len(line.strip().split("\t")) < 2:
                continue
            coverage = int(line.strip().split("\t")[2])
            coverage_sum += coverage
            coverage_nr_pos += 1
            coverage_list.append(coverage)
    outfile.write(gene + "\t" + str(round(coverage_sum/float(coverage_nr_pos), 1)) + "\n")
outfile.close()
