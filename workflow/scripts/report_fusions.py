
import sys
import subprocess

input_bed = open(snakemake.input.bed)
input_bed_extra_annotation = open(snakemake.input.bed_extra_annotation)
input_arriba = open(snakemake.input.arriba)
input_starfusion = open(snakemake.input.star_fusion)
input_fusioncatcher = open(snakemake.input.fusioncatcher)
input_bam = snakemake.input.bam
output_fusions = open(snakemake.output.fusions, "w")
star_fusion_flag_low_support = snakemake.params.star_fusion_flag_low_support
star_fusion_low_support = snakemake.params.star_fusion_low_support
star_fusion_low_support_inframe = snakemake.params.star_fusion_low_support_inframe
fusioncatcher_flag_low_support = snakemake.params.fusioncatcher_flag_low_support
fusioncatcher_low_support = snakemake.params.fusioncatcher_low_support
fusioncatcher_low_support_inframe = snakemake.params.fusioncatcher_low_support_inframe
fp_fusions_filename = snakemake.params.fp_fusions


housekeeping_genes = {}
artefact_gene_dict = {}
if fp_fusions_filename != "":
    with open(fp_fusions_filename) as fp_fusions:
        for line in fp_fusions:
            if line.startswith("#"):
                continue
            columns = line.strip().split("\t")
            gene1 = columns[0]
            gene2 = columns[1]
            read_limit_SF = int(columns[2])
            read_limit_FC = int(columns[3])
            if gene2 == "housekeeping":
                housekeeping_genes[gene1] = [read_limit_SF, read_limit_FC]
            if gene1 not in artefact_gene_dict:
                artefact_gene_dict[gene1] = {}
            artefact_gene_dict[gene1][gene2] = [read_limit_SF, read_limit_FC]

output_fusions.write("caller\tgene1\tgene2\texon1\texon2\tconfidence\tFC-callers\tpredicted_effect\tbreakpoint1\tbreakpoint2\t")
output_fusions.write("coverage1\tcoverage2\tsplit_reads\tspanning_pairs\ttotal_supporting_reads\n")

# Only keep fusions with one gene that are in the design
design_genes = {}
for line in input_bed:
    lline = line.strip().split("\t")
    chrom = lline[0]
    start = int(lline[1])
    end = int(lline[2])
    exon = lline[3]
    gene = lline[3].split("_")[0]
    if gene in design_genes:
        design_genes[gene].append([chrom, start, end, exon])
    else:
        design_genes[gene] = [[chrom, start, end, exon]]

# Extra annotation of partner genes
annotation_genes = {}
for line in input_bed_extra_annotation:
    lline = line.strip().split("\t")
    chrom = lline[0]
    start = int(lline[1])
    end = int(lline[2])
    exon = lline[3]
    gene = lline[3].split("_")[0]
    if gene in design_genes:
        design_genes[gene].append([chrom, start, end, exon])
    else:
        design_genes[gene] = [[chrom, start, end, exon]]

# Arriba fusions
header = True
for line in input_arriba:
    if header:
        header = False
        continue
    lline = line.strip().split("\t")
    gene1 = lline[0]
    gene2 = lline[1]
    # Only keep fusions with one gene that are in the design
    if (gene1 not in design_genes and gene2 not in design_genes):
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
    # Compare fusion coverage with coverage in breakpoints
    chrom1 = "chr" + breakpoint1.split(":")[0]
    pos1 = breakpoint1.split(":")[1]
    chrom2 = "chr" + breakpoint2.split(":")[0]
    pos2 = breakpoint2.split(":")[1]
    # Get exon name if it is in design
    exon1 = ""
    exon2 = ""
    if gene1 in design_genes:
        for region in design_genes[gene1]:
            if int(pos1) >= region[1] and int(pos1) <= region[2]:
                exon1 = region[3]
    elif gene1 in annotation_genes:
        for region in annotation_genes[gene1]:
            if int(pos1) >= region[1] and int(pos1) <= region[2]:
                exon1 = region[3]
    if gene2 in design_genes:
        for region in design_genes[gene2]:
            if int(pos2) >= region[1] and int(pos2) <= region[2]:
                exon2 = region[3]
    elif gene2 in annotation_genes:
        for region in annotation_genes[gene2]:
            if int(pos2) >= region[1] and int(pos2) <= region[2]:
                exon2 = region[3]
    total_supporting_reads = int(total_split_reads) + int(discordant_mates)
    output_fusions.write(f"Arriba\t{gene1}\t{gene2}\t{exon1}\t{exon2}\t{confidence}\t\t{predicted_effect}\t{breakpoint1}")
    output_fusions.write(f"\t{breakpoint2}\t{coverage1}\t{coverage2}\t{total_split_reads}\t{discordant_mates}")
    output_fusions.write(f"\t{total_supporting_reads}\n")


# Star-fusions
header = True
for line in input_starfusion:
    if header:
        header = False
        continue
    lline = line.strip().split("\t")
    gene1 = lline[0].split("--")[0]
    gene2 = lline[0].split("--")[1]
    # Only keep fusions with one gene that are in the design
    if (gene1 not in design_genes and gene2 not in design_genes):
        continue
    Junction_read_count = lline[1]
    Spanning_Frag_count = lline[2]
    predicted_effect = lline[21]
    # Flag fusions with junction_read_count < 15 and Spanning_Frag_count < 2
    confidence = ""
    if int(Junction_read_count) < star_fusion_flag_low_support:
        confidence = "Low support"
    # Remove Fusions with very weak read support
    if int(Junction_read_count) <= star_fusion_low_support_inframe and predicted_effect != "INFRAME":
        continue
    if int(Junction_read_count) <= star_fusion_low_support:
        continue
    # Higher demand of read support for genes with frequent FP, house keeping genes
    if (gene1 in artefact_gene_dict and gene2 in artefact_gene_dict[gene1]):
        if int(Junction_read_count) < artefact_gene_dict[gene1][gene2][0]:
            continue
    if (gene2 in artefact_gene_dict and gene1 in artefact_gene_dict[gene2]):
        if int(Junction_read_count) < artefact_gene_dict[gene2][gene1][0]:
            continue
    if gene1 in housekeeping_genes:
        if int(Junction_read_count) < housekeeping_genes[gene1]["housekeeping"][0]:
            continue
    if gene2 in housekeeping_genes:
        if int(Junction_read_count) < housekeeping_genes[gene2]["housekeeping"][0]:
            continue
    breakpoint1 = lline[7]
    breakpoint2 = lline[9]
    FFPM = lline[11]
    DBs = lline[16]
    # Compare fusion coverage with coverage in breakpoints
    chrom1 = breakpoint1.split(":")[0]
    pos1 = breakpoint1.split(":")[1]
    chrom2 = breakpoint2.split(":")[0]
    pos2 = breakpoint2.split(":")[1]
    # Get exon name if it is in design
    exon1 = ""
    exon2 = ""
    if gene1 in design_genes:
        for region in design_genes[gene1]:
            if int(pos1) >= region[1] and int(pos1) <= region[2]:
                exon1 = region[3]
    if gene2 in design_genes:
        for region in design_genes[gene2]:
            if int(pos2) >= region[1] and int(pos2) <= region[2]:
                exon2 = region[3]
    total_supporting_reads = int(Junction_read_count) + int(Spanning_Frag_count)
    output_fusions.write(f"StarFusion\t{gene1}\t{gene2}\t{exon1}\t{exon2}\t{confidence}\t\t{predicted_effect}\t{breakpoint1}")
    output_fusions.write(f"\t{breakpoint2}\t\t\t{Junction_read_count}\t{Spanning_Frag_count}\t{total_supporting_reads}\n")


# FusionCatcher
header = True
for line in input_fusioncatcher:
    if header:
        header = False
        continue
    lline = line.strip().split("\t")
    gene1 = lline[0]
    gene2 = lline[1]
    # Only keep fusions with one gene that are in the design
    if (gene1 not in design_genes and gene2 not in design_genes):
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
    # Flag fusions with Spanning_reads_unique < 5
    confidence = ""
    if int(Spanning_reads_unique) < fusioncatcher_flag_low_support:
        confidence = "Low support"
    # Filter fusions with very low support
    if int(Spanning_reads_unique) <= fusioncatcher_low_support_inframe and predicted_effect != "in-frame":
        continue
    if int(Spanning_reads_unique) <= fusioncatcher_low_support:
        continue
    # Higher demand of read support for genes with frequent FP, house keeping genes, and pool2 genes without fusion to pool1 gene
    if (gene1 in artefact_gene_dict and gene2 in artefact_gene_dict[gene1]):
        if int(Spanning_reads_unique) < artefact_gene_dict[gene1][gene2][1]:
            continue
    if (gene2 in artefact_gene_dict and gene1 in artefact_gene_dict[gene2]):
        if int(Spanning_reads_unique) < artefact_gene_dict[gene2][gene1][1]:
            continue
    if gene1 in housekeeping_genes:
        if int(Spanning_reads_unique) < housekeeping_genes[gene1]["housekeeping"][1]:
            continue
    if gene2 in housekeeping_genes:
        if int(Spanning_reads_unique) < housekeeping_genes[gene2]["housekeeping"][1]:
            continue
    # Flag fusions annotated that are fusions with very high probability
    fp_db = [
        "banned", "bodymap2", "cacg", "1000genomes", "conjoing", "cortex", "distance1000bp", "ensembl_fully_overlapping",
        "ensembl_same_strand_overlapping", "gtex", "hpa", "mt", "paralogs", "refseq_fully_overlapping",
        "refseq_same_strand_overlapping", "rrna", "similar_reads", "similar_symbols", "ucsc_fully_overlapping",
        "ucsc_same_strand_overlapping"
    ]
    fp_found = ""
    for fp in fp_db:
        if fp in fp_filters:
            fp_found = "FP"
    # Compare fusion coverage with coverage in breakpoints
    pos1 = "0"
    pos2 = "0"
    if len(breakpoint1.split(":")) == 3 and len(breakpoint2.split(":")) == 3:
        chrom1 = "chr" + breakpoint1.split(":")[0]
        pos1 = breakpoint1.split(":")[1]
        chrom2 = "chr" + breakpoint2.split(":")[0]
        pos2 = breakpoint2.split(":")[1]
    # Get exon name if it is in design
    exon1 = ""
    exon2 = ""
    if gene1 in design_genes:
        for region in design_genes[gene1]:
            if int(pos1) >= region[1] and int(pos1) <= region[2]:
                exon1 = region[3]
    if gene2 in design_genes:
        for region in design_genes[gene2]:
            if int(pos2) >= region[1] and int(pos2) <= region[2]:
                exon2 = region[3]
    total_supporting_reads = int(Spanning_pairs) + int(Spanning_reads_unique)
    output_fusions.write(f"FusionCatcher\t{gene1}\t{gene2}\t{exon1}\t{exon2}\t{confidence}\t{Fusion_finding_method}")
    output_fusions.write(f"\t{predicted_effect}\t{breakpoint1}\t{breakpoint2}\t\t\t{Spanning_pairs}")
    output_fusions.write(f"\t{Spanning_reads_unique}\t{total_supporting_reads}\n")
