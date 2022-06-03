
import sys
import subprocess

input_bed = open(snakemake.input.bed)
input_arriba = open(snakemake.input.arriba)
input_starfusion = open(snakemake.input.star_fusion)
input_fusioncatcher = open(snakemake.input.fusioncathcer)
input_bam = snakemake.input.bam
output_fusions = open(snakemake.output.fusions, "w")
output_coverage_file_name = snakemake.output.coverage


housekeeping_genes = ["GAPDH", "GUSB", "OAZ1", "POLR2A"]
artefact_genes = {"MAML2" : ["FRMPD3", "NCOA6", "ATXN3", "SRP14", "KMT2D", "CHD1", "NFAT5", "FOXP2", "NUMBL", "GLG1", "VEZF1", "AAK1", "NCOR2"]}

output_fusions.write("Caller\tgene1\tgene2\texon1\texon2\tconfidence\tpredicted_effect\tbreakpoint1\tbreakpoint2\tcoverage1\tcoverage2\tsplit_reads\tSpanning_pairs\tBreakpoint1_covarage/SplitReads\tBreakpoint2_covarage/SplitReads\n")

#Only keep fusions with one gene that are in the design
design_genes = {}
for line in input_bed :
    lline = line.strip().split("\t")
    chrom = lline[0]
    start = int(lline[1])
    end = int(lline[2])
    exon = lline[3]
    gene = lline[3].split("_")[0]
    if len(lline) > 4 :
        pool = lline[4]
    if gene in design_genes :
        design_genes[gene].append([chrom, start, end, exon])
    else :
        design_genes[gene] = [[chrom, start, end, exon]]

#Arriba fusions
header = True
for line in input_arriba :
    if header :
        header = False
        continue
    lline = line.strip().split("\t")
    gene1 = lline[0]
    gene2 = lline[1]
    #Only keep fusions with one gene that are in the design
    if (gene1 not in design_genes and gene2 not in design_genes) :
        continue
    confidence = lline[14]
    breakpoint1 = lline[4]
    breakpoint2 = lline[5]
    split_reads1 = lline[9]
    split_reads2 = lline[10]
    total_split_reads = int(split_reads1) + int(split_reads2)
    discordant_mates = lline[11]
    coverage1 = lline[12]
    coverage2 = lline[13]
    predicted_effect = lline[15]
    #Compare fusion coverage with coverage in breakpoints
    chrom1 = "chr" + breakpoint1.split(":")[0]
    pos1 = breakpoint1.split(":")[1]
    chrom2 = "chr" + breakpoint2.split(":")[0]
    pos2 = breakpoint2.split(":")[1]
    cov1 = 0
    cov2 = 0
    subprocess.call("samtools depth -d 500000 -a -r " + chrom1 + ":" + pos1 + "-" + pos1 + " " + input_bam + " > " + output_coverage_file_name, shell=True)
    output_coverage = open(output_coverage_file_name)
    for line in output_coverage :
        cov1 = int(line.strip().split("\t")[2])
    output_coverage.close()
    subprocess.call("samtools depth -d 500000 -a -r " + chrom2 + ":" + pos2 + "-" + pos2 + " " + input_bam + " > " + output_coverage_file_name, shell=True)
    output_coverage = open(output_coverage_file_name)
    for line in output_coverage :
        cov2 = int(line.strip().split("\t")[2])
    output_coverage.close()
    q1 = round((cov1 / (float(split_reads1) + float(split_reads2) + 0.01)),1)
    q2 = round((cov2 / (float(split_reads1) + float(split_reads2) + 0.01)),1)
    #Get exon name if it is in design
    exon1 = ""
    exon2 = ""
    if gene1 in design_genes :
        for region in design_genes[gene1] :
            if int(pos1) >= region[1] and int(pos1) <= region[2] :
                exon1 = region[3]
    if gene2 in design_genes :
        for region in design_genes[gene2] :
            if int(pos2) >= region[1] and int(pos2) <= region[2] :
                exon2 = region[3]
    output_fusions.write("Arriba\t" + gene1 + "\t" + gene2 + "\t" + exon1 + "\t" + exon2 + "\t" + confidence + "\t" + predicted_effect + "\t" + breakpoint1 + "\t" + breakpoint2 + "\t" + coverage1 + "\t" + coverage2 + "\t" + str(total_split_reads) + "\t" + discordant_mates + "\t" + str(q1) + "\t" + str(q2) + "\n")


#Star-fusions
header = True
for line in input_starfusion :
    if header :
        header = False
        continue
    lline = line.strip().split("\t")
    gene1 = lline[0].split("--")[0]
    gene2 = lline[0].split("--")[1]
    #Only keep fusions with one gene that are in the design
    if (gene1 not in design_genes and gene2 not in design_genes) :
        continue
    Junction_read_count = lline[1]
    Spanning_Frag_count = lline[2]
    predicted_effect = lline[21]
    #Flag fusions with junction_read_count < 15 and Spanning_Frag_count < 2
    confidence = ""
    if int(Junction_read_count) < 15 :
        confidence = "Low support"
    #Remove Fusions with very weak read support
    if int(Junction_read_count) <= 6 and predicted_effect != "INFRAME":
        continue
    if int(Junction_read_count) <= 2 :
        continue
    #Higher demand of read support for genes with frequent FP, house keeping genes
    if (((gene1 in artefact_genes and gene2 in artefact_genes[gene1]) or (gene2 in artefact_genes and gene1 in artefact_genes[gene2])) or
        gene1 in housekeeping_genes or gene2 in housekeeping_genes) :
        if int(Junction_read_count) < 20 :
            continue
    breakpoint1 = lline[7]
    breakpoint2 = lline[9]
    FFPM = lline[11]
    DBs = lline[16]
    #Compare fusion coverage with coverage in breakpoints
    chrom1 = breakpoint1.split(":")[0]
    pos1 = breakpoint1.split(":")[1]
    chrom2 = breakpoint2.split(":")[0]
    pos2 = breakpoint2.split(":")[1]
    cov1 = 0
    cov2 = 0
    subprocess.call("samtools depth -d 500000 -a -r " + chrom1 + ":" + pos1 + "-" + pos1 + " " + input_bam + " > " + output_coverage_file_name, shell=True)
    output_coverage = open(output_coverage_file_name)
    for line in output_coverage :
        cov1 = int(line.strip().split("\t")[2])
    output_coverage.close()
    subprocess.call("samtools depth -d 500000 -a -r " + chrom2 + ":" + pos2 + "-" + pos2 + " " + input_bam + " > " + output_coverage_file_name, shell=True)
    output_coverage = open(output_coverage_file_name)
    for line in output_coverage :
        cov2 = int(line.strip().split("\t")[2])
    output_coverage.close()
    q1 = round((cov1 / (float(Junction_read_count) + 0.01)),1)
    q2 = round((cov2 / (float(Junction_read_count) + 0.01)),1)
    #Get exon name if it is in design
    exon1 = ""
    exon2 = ""
    if gene1 in design_genes :
        for region in design_genes[gene1] :
            if int(pos1) >= region[1] and int(pos1) <= region[2] :
                exon1 = region[3]
    if gene2 in design_genes:
        for region in design_genes[gene2] :
            if int(pos2) >= region[1] and int(pos2) <= region[2] :
                exon2 = region[3]
    output_fusions.write("StarFusion\t" + gene1 + "\t" + gene2 + "\t" + exon1 + "\t" + exon2 + "\t" + confidence + "\t" + predicted_effect + "\t" + breakpoint1 + "\t" + breakpoint2 + "\t\t\t" + Junction_read_count + "\t" + Spanning_Frag_count + "\t" + str(q1) + "\t" + str(q2) + "\n")



#FusionCatcher
header = True
for line in input_fusioncatcher :
    if header :
        header = False
        continue
    lline = line.strip().split("\t")
    gene1 = lline[0]
    gene2 = lline[1]
    #Only keep fusions with one gene that are in the design
    if (gene1 not in design_genes and gene2 not in design_genes) :
        continue
    fp_filters = lline[2].split(",")
    DBs = lline[2]
    common_mapping = lline[3]
    Spanning_pairs = lline[4]
    Spanning_reads_unique = lline[5]
    Fusion_finding_method = lline[7]
    breakpoint1 = lline[8]
    breakpoint2 = lline[9]
    predicted_effect = lline[15]
    #Flag fusions with Spanning_reads_unique < 5
    confidence = ""
    if int(Spanning_reads_unique) < 15:
        confidence = "Low support"
    #Filter fusions with very low support
    if int(Spanning_reads_unique) <= 6 and predicted_effect != "in-frame":
        continue
    if int(Spanning_reads_unique) <= 3 :
        continue
    #Higher demand of read support for genes with frequent FP, house keeping genes, and pool2 genes without fusion to pool1 gene
    if (((gene1 in artefact_genes and gene2 in artefact_genes[gene1]) or (gene2 in artefact_genes and gene1 in artefact_genes[gene2])) or
        gene1 in housekeeping_genes or gene2 in housekeeping_genes) :
        if int(Spanning_reads_unique) < 20 :
            continue
    #Flag fusions annotated that are fusions with very high probability
    fp_db = ["banned", "bodymap2", "cacg", "1000genomes", "conjoing", "cortex", "distance1000bp", "ensembl_fully_overlapping", "ensembl_same_strand_overlapping", "gtex", "hpa", "mt", "paralogs", "refseq_fully_overlapping", "refseq_same_strand_overlapping", "rrna", "similar_reads", "similar_symbols", "ucsc_fully_overlapping", "ucsc_same_strand_overlapping"]
    fp_found = ""
    for fp in fp_db :
        if fp in fp_filters :
            fp_found = "FP"
    #Compare fusion coverage with coverage in breakpoints
    pos1 = "0"
    pos2 = "0"
    if len(breakpoint1.split(":")) == 3 and len(breakpoint2.split(":")) == 3 :
        chrom1 = "chr" + breakpoint1.split(":")[0]
        pos1 = breakpoint1.split(":")[1]
        chrom2 = "chr" + breakpoint2.split(":")[0]
        pos2 = breakpoint2.split(":")[1]
        cov1 = 0
        cov2 = 0
        subprocess.call("samtools depth -d 500000 -a -r " + chrom1 + ":" + pos1 + "-" + pos1 + " " + input_bam + " > " + output_coverage_file_name, shell=True)
        output_coverage = open(output_coverage_file_name)
        for line in output_coverage :
            cov1 = int(line.strip().split("\t")[2])
        output_coverage.close()
        subprocess.call("samtools depth -d 500000 -a -r " + chrom2 + ":" + pos2 + "-" + pos2 + " " + input_bam + " > " + output_coverage_file_name, shell=True)
        output_coverage = open(output_coverage_file_name)
        for line in output_coverage :
            cov2 = int(line.strip().split("\t")[2])
        output_coverage.close()
        q1 = round((cov1 / (float(Spanning_pairs) + 0.01)),1)
        q2 = round((cov2 / (float(Spanning_pairs) + 0.01)),1)
    else :
        q1 = "NA"
        q2 = "NA"
    #Get exon name if it is in design
    exon1 = ""
    exon2 = ""
    if gene1 in design_genes :
        for region in design_genes[gene1] :
            if int(pos1) >= region[1] and int(pos1) <= region[2] :
                exon1 = region[3]
    if gene2 in design_genes :
        for region in design_genes[gene2] :
            if int(pos2) >= region[1] and int(pos2) <= region[2] :
                exon2 = region[3]
    output_fusions.write("FusionCatcher\t" + gene1 + "\t" + gene2 + "\t" + exon1 + "\t" + exon2 + "\t" + confidence + "\t" + predicted_effect + "\t" + breakpoint1 + "\t" + breakpoint2 + "\t\t\t" + Spanning_pairs + "\t" + Spanning_reads_unique + "\t" + str(q1) + "\t" + str(q2) + "\n")
