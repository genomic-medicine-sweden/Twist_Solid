
import subprocess
import sys
import gzip

bed = open(snakemake.input.hotspots)
vcf = open(snakemake.input.vcf)
bam_file = snakemake.input.bam
background_panel_filename = snakemake.input.background_panel
background_run = open(snakemake.input.background_run)
gvcf = snakemake.input.gvcf
min_coverage = snakemake.params.min_coverage
outfile = open(snakemake.output.low_coverage, "w")
outfile2 = open(snakemake.output.hotspot_info, "w")

header1 = "#Chr\tStart_hg19\tEnd_hg19\tGene\tCDS_mut_syntax\tAA_mut_syntax\tReport\tcomment\tExon\tAccession_number"
header1 += "\tCoverage\tPosition"
header2 = header1
header2 += "\tpanel_median\tpanel_sd\trun_median\talt_AF\tSD_from_median"
header2 += "\tDP\tRef_DP\tAlt_DP\tAF\tAA_change\tCDS_change\n"

outfile.write(header1 + "\n")
outfile2.write(header2)


'''Find positions to report and gene regions to analyse'''
inv_pos = {}
gene_regions = []
gene_region_dict = {}
hotspot_dict = {}
header = True
prev_gene = ""
prev_chrom = ""
first_gene = True
gene_start_pos = 0
gene_end_pos = 0
for line in bed:
    if header:
        header = False
        continue
    lline = line.strip().split("\t")
    Report = lline[6]
    if Report == "region":
        continue
    chrom = str(int(lline[0].split(".")[0].split("_")[1]))
    start_pos = lline[1]
    end_pos = lline[2]
    gene = lline[3]
    pos = int(start_pos)
    while pos <= int(end_pos):
        inv_pos[chrom + "_" + str(pos)] = lline
        if Report == "hotspot" or Report == "region_all":
            hotspot_dict[chrom + "_" + str(pos)] = ""
        pos += 1
    if first_gene:
        prev_gene = gene
        prev_chrom = chrom
        gene_start_pos = start_pos
        gene_end_pos = end_pos
        first_gene = False
        continue
    if prev_gene != gene:
        gene_regions.append([prev_chrom, int(gene_start_pos), int(gene_end_pos)])
        prev_gene = gene
        prev_chrom = chrom
        gene_start_pos = start_pos
        gene_end_pos = end_pos
    else:
        gene_end_pos = end_pos
gene_regions.append([prev_chrom, int(gene_start_pos), int(gene_end_pos)])


'''find all calls in the vcf overlapping hotspots'''
vcf_dict = {}
header = True
for line in vcf:
    if header:
        if line[:6] == "#CHROM":
            header = False
        continue
    lline = line.strip().split("\t")
    chrom = lline[0][3:]
    pos = lline[1]
    key = chrom + "_" + pos
    INFO = lline[7].split(";")
    FORMAT = lline[8].split(":")
    DATA = lline[9].split(":")
    if INFO[:3] == "AA=":
        continue
    AD_index = 0
    DP_index = 0
    RD_index = 0
    i = 0
    for f in FORMAT:
        if f == "AD":
            AD_index = i
        if f == "DP":
            DP_index = i
        if f == "RD":
            RD_index = i
        i += 1
    AD = DATA[AD_index].split(",")
    Ref_DP = 0
    Alt_DP = 0
    if len(AD) == 2:
        Ref_DP = AD[0]
        Alt_DP = AD[1]
    else:
        Ref_DP = DATA[RD_index]
        Alt_DP = DATA[AD_index]
    DP = DATA[DP_index]
    AF_index = 0
    i = 0
    for info in INFO:
        if info[:3] == "AF=":
            AF_index = i
        i += 1
    AF = INFO[AF_index][3:]
    VEP = INFO[-1]
    AA_change = VEP.split(":p.")
    if len(AA_change) == 2:
        AA_change = AA_change[1].split("|")[0]
        if AA_change[-3:] == "%3D":
            AA_change = AA_change[:-3]
    else:
        AA_change = ""
    CDS_change = VEP.split(":c.")
    if len(CDS_change) == 2:
        CDS_change = CDS_change[1].split("|")[0]
    else:
        CDS_change = ""
    if key in inv_pos:
        vcf_dict[key] = [DP, Ref_DP, Alt_DP, AF, AA_change, CDS_change]


'''find all positions in the gvcfs overlapping hotspots (not including indels)'''
gvcf_panel_dict = {}
gvcf_run_dict = {}
gvcf_sample_dict = {}
if background_panel_filename != "":
    background_panel = open(background_panel_filename)
    next(background_panel)
    for line in background_panel:
        columns = line.strip().split()
        chrom = columns[0]
        pos = columns[1]
        key = chrom + "_" + pos
        if key in hotspot_dict:
            median = float(columns[2])
            sd = float(columns[3])
            gvcf_panel_dict[key] = [median, sd]
next(background_run)
for line in background_run:
    columns = line.strip().split()
    chrom = columns[0]
    pos = columns[1]
    key = chrom + "_" + pos
    if key in hotspot_dict:
        median = float(columns[2])
        gvcf_run_dict[key] = median
with gzip.open(gvcf, 'rt') as infile:
    file_content = infile.read().split("\n")
    header = True
    for line in file_content:
        if header:
            if line[:6] == "#CHROM":
                header = False
            continue
        columns = line.strip().split("\t")
        if len(columns) <= 1:
            continue
        chrom = columns[0][3:]
        pos = columns[1]
        key = chrom + "_" + pos
        if key in hotspot_dict:
            format = columns[8].split(":")
            data = columns[9].split(":")
            AD_id = 0
            DP_mosdepth_id = 0
            i = 0
            for f in format:
                if f == "AD":
                    AD_id = i
                if f == "DP_mosdepth":
                    DP_mosdepth_id = i
                i +=1
            AD_info = data[AD_id].split(",")
            DP_mosdepth = data[DP_mosdepth_id]
            ref_AD = int(AD_info[0])
            alt_AD = 0
            for AD in AD_info[1:]:
                alt_AD += int(AD)
            DP = ref_AD + alt_AD
            alt_AF = 0.0
            if DP > 0:
                alt_AF = alt_AD / float(DP)
            gvcf_sample_dict[key] = [alt_AF, DP, DP_mosdepth]


'''Report all interesting positions and also all positions (All but region) with coverage < 200'''
depth_dict = {}
for region in gene_regions:
    start_pos = region[1]
    end_pos = region[2]
    chrom = region[0]
    while start_pos <= end_pos:
        pos = str(start_pos)
        key = chrom + "_" + pos
        if key in inv_pos:
            coverage = 0
            if key in gvcf_sample_dict:
                coverage = int(gvcf_sample_dict[key][2])
            if coverage < min_coverage:
                for info in inv_pos[key]:
                    outfile.write(info + "\t")
                outfile.write(str(coverage) + "\t" + pos + "\n")
            for info in inv_pos[key]:
                outfile2.write(info + "\t")
            outfile2.write(str(coverage) + "\t" + pos)
            panel_median = 1000
            panel_sd = 1000
            run_median = 1000
            alt_AF = 0.0
            pos_sd = 1000
            if key in gvcf_panel_dict:
                panel_median = gvcf_panel_dict[key][0]
                panel_sd = gvcf_panel_dict[key][1]
            if key in gvcf_run_dict:
                run_median = gvcf_run_dict[key]
            if key in gvcf_sample_dict:
                alt_AF = gvcf_sample_dict[key][0]
                if panel_sd > 0.0:
                    pos_sd = (alt_AF - panel_median) / panel_sd
            found_hotspot = False
            if key in hotspot_dict and ((alt_AF >= 0.005 and pos_sd >= 3.0) or key in vcf_dict):
                found_hotspot = True
                outfile2.write(
                    "\t" + "{:.4f}".format(panel_median) + "\t" + "{:.4f}".format(panel_sd) + "\t" + "{:.4f}".format(run_median) +
                    "\t" + "{:.4f}".format(alt_AF) + "\t" + "{:.2f}".format(pos_sd)
                )
            if key in vcf_dict:
                if not found_hotspot:
                    outfile2.write("\t\t\t\t\t")
                for info in vcf_dict[key]:
                    outfile2.write("\t" + str(info))
            outfile2.write("\n")
        start_pos += 1
outfile.close()
outfile2.close()
