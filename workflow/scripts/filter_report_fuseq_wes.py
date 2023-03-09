
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


def filter_report_fusion(sample, fusion_breakpoint_dict, report_genes, fusion_file, min_support, filter_on_fusiondb, out_file):
    out_file.write("fusion\tbreak_point1\tbreak_point2\tparalog\tsplit_reads\tmate_pairs\ttotal_support\n")
    nr_report_genes = len(report_genes)
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
            if fusion_name in fusion_breakpoint_dict:
                break_points = fusion_breakpoint_dict[fusion_name]
            if reverse_fusion_name in fusion_breakpoint_dict:
                break_points = fusion_breakpoint_dict[reverse_fusion_name]
            break_point1 = break_points[0]
            break_point2 = break_points[1]
            out_file.write(
                f"{fusion_name}\t{break_point1}\t{break_point2}\t{paralog}\t"
                f"{SR_support}\t{MR_support}\t{support}\n"
            )
    out_file.close()


if __name__ == "__main__":
    log = snakemake.log_fmt_shell(stdout=False, stderr=True)

    sample = snakemake.input.breakpoint.split("/")[-2].split(".")[0]
    fusion_breakpoint_dict = get_breakpoints(open(snakemake.input.breakpoint), sample)
    report_genes = get_report_genes(open(snakemake.params.gene_white_list))
    filter_report_fusion(
        sample,
        fusion_breakpoint_dict,
        report_genes,
        open(snakemake.input.fusions),
        snakemake.params.min_support,
        snakemake.params.filter_on_fusiondb,
        open(snakemake.output.fusions, "w")
    )
