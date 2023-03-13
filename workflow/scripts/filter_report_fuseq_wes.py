
import logging

log = logging.getLogger()


def get_breakpoints(breakpoint_file, sample):
    fusion_breakpoint_dict = {}
    next(breakpoint_file)
    for fusion in breakpoint_file:
        columns = fusion.strip().split("\t")
        fusion_name = columns[18]
        chrom1 = columns[3].split("__")[1]
        break_point1_1 = int(columns[3].split("__")[2]) + int(columns[5])
        break_point1_2 = break_point1_1 + int(columns[7])
        break_point1 = f"{chrom1}:{break_point1_1}-{break_point1_2}"
        chrom2 = columns[8].split("__")[1]
        break_point2_1 = int(columns[8].split("__")[2]) + int(columns[10])
        break_point2_2 = break_point2_1 + int(columns[12])
        break_point2 = f"{chrom2}:{break_point2_1}-{break_point2_2}"
        fusion_breakpoint_dict[fusion_name] = [break_point1, break_point2]
    return fusion_breakpoint_dict


def get_report_genes(gene_white_list):
    report_genes = []
    for gene in gene_white_list:
        report_genes.append(gene.strip())
    return report_genes


def filter_fusion(sample, fusion_breakpoint_dict, report_genes, fusion_file, min_support, filter_on_fusiondb):
    nr_report_genes = len(report_genes)
    filtered_fusions = []
    next(fusion_file)
    for fusion in fusion_file:
        columns = fusion.strip().split("\t")
        fusion_name = columns[0]
        gene1 = fusion_name.split("--")[0]
        gene2 = fusion_name.split("--")[1]
        reverse_fusion_name = f"{gene2}--{gene1}"
        if nr_report_genes > 0 and not (gene1 in report_genes or gene2 in report_genes):
            continue
        support = int(columns[8])
        if support < min_support:
            continue
        SR_support = int(columns[9])
        MR_support = int(columns[10])
        fusiondb = int(columns[11])
        paralog = columns[12]
        if fusiondb == 1 or filter_on_fusiondb is False:
            break_points = ["", ""]
            break_point1 = ""
            break_point2 = ""
            if SR_support > 0:
                if fusion_name in fusion_breakpoint_dict:
                    break_points = fusion_breakpoint_dict[fusion_name]
                if reverse_fusion_name in fusion_breakpoint_dict:
                    break_points = fusion_breakpoint_dict[reverse_fusion_name]
                break_point1 = break_points[0]
                break_point2 = break_points[1]
            filtered_fusions.append([fusion_name, break_point1, "", break_point2, "", paralog, SR_support, MR_support, support])
    return filtered_fusions


def annotate_fusion(filtered_fusions, input_gtf):
    gene_dict = {}
    chr_pos_dict = {}
    transcript_exon_max = {}
    i = 0
    for fusion in filtered_fusions:
        if not (fusion[1] == "" or fusion[3] == ""):
            gene1 = fusion[0].split("--")[0]
            gene2 = fusion[0].split("--")[1]
            gene_dict[gene1] = ""
            gene_dict[gene2] = ""
            bp1_chrom = fusion[1].split(":")[0]
            bp1_pos = (int(fusion[1].split(":")[1].split("-")[0]) + int(fusion[1].split(":")[1].split("-")[1])) / 2
            bp2_chrom = fusion[3].split(":")[0]
            bp2_pos = (int(fusion[3].split(":")[1].split("-")[0]) + int(fusion[3].split(":")[1].split("-")[1])) / 2
            if bp1_chrom in chr_pos_dict:
                chr_pos_dict[bp1_chrom].append([bp1_pos, i, 2, 100000, 0, "", ""])
            else:
                chr_pos_dict[bp1_chrom] = [[bp1_pos, i, 2, 100000, 0, "", ""]]
            if bp2_chrom in chr_pos_dict:
                chr_pos_dict[bp2_chrom].append([bp2_pos, i, 4, 100000, 0, "", ""])
            else:
                chr_pos_dict[bp2_chrom] = [[bp2_pos, i, 4, 100000, 0, "", ""]]
        i += 1
    for gtf in input_gtf:
        columns = gtf.strip().split("\t")
        type = columns[2]
        if type != "CDS":
            continue
        gene_name = columns[8].split("gene_id \"")[1].split("\";")[0]
        if gene_name in gene_dict:
            chrom = columns[0]
            pos = (int(columns[3]) + int(columns[4])) / 2
            if chrom in chr_pos_dict:
                for bp in chr_pos_dict[chrom]:
                    direction = columns[6]
                    if direction == "-":
                        distance = bp[0] - pos
                    else:
                        distance = pos - bp[0]
                    if distance < 0:
                        distance = 100000
                    transcript_id = columns[8].split("transcript_id \"")[1].split("\";")[0]
                    if transcript_id == "NM_001353765":
                        continue
                    exon_number = int(columns[8].split("exon_number \"")[1].split("\";")[0])
                    if distance < bp[3]:
                        bp[3] = distance
                        bp[4] = exon_number
                        bp[5] = transcript_id
                        bp[6] = direction
                    transcript_exon_max[transcript_id] = exon_number
    for chrom in chr_pos_dict:
        for bp in chr_pos_dict[chrom]:
            transcript_id = bp[5]
            if bp[6] == "-":
                exon_number = transcript_exon_max[transcript_id] - bp[4] + 1
            else:
                exon_number = bp[4]
            filtered_fusions[bp[1]][bp[2]] = f"exon {exon_number} in {transcript_id}"
    return filtered_fusions


def write_fusions(annotated_filtered_fusions, out_file):
    out_file.write("fusion\tbreak_point1\texon1\tbreak_point2\texon2\tparalog\tsplit_reads\tmate_pairs\ttotal_support\n")
    for data in annotated_filtered_fusions:
        first = True
        for d in data:
            if first:
                out_file.write(f"{d}")
                first = False
            else:
                out_file.write(f"\t{d}")
        out_file.write(f"\n")
    out_file.close()


if __name__ == "__main__":
    log = snakemake.log_fmt_shell(stdout=False, stderr=True)

    sample = snakemake.input.breakpoint.split("/")[-2].split(".")[0]
    fusion_breakpoint_dict = get_breakpoints(open(snakemake.input.breakpoint), sample)
    report_genes = get_report_genes(open(snakemake.params.gene_white_list))
    filtered_fusions = filter_fusion(
        sample,
        fusion_breakpoint_dict,
        report_genes,
        open(snakemake.input.fusions),
        snakemake.params.min_support,
        snakemake.params.filter_on_fusiondb,
    )
    annotated_filtered_fusions = annotate_fusion(filtered_fusions, open(snakemake.params.gtf))
    write_fusions(annotated_filtered_fusions, open(snakemake.output.fusions, "w"))
