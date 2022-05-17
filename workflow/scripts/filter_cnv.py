
from pysam import VariantFile

def variant_in_genelist(chrom, start, end, gene_dict):
    keep_variant = False
    genes = ""
    if chrom in gene_dict:
        for gene_region in gene_dict[chrom]:
            g_start = gene_region[0]
            g_end = gene_region[1]
            if ((start >= g_start and start <= g_end) or
                (end >= g_start and end <= g_end) or
                (start < g_start and end > g_end)):
                keep_variant = True
                if genes == "":
                    genes = gene_region[2]
                else:
                    genes = "%s,%s" % (genes, gene_region[2])
    return keep_variant, genes


def filter_variants(in_vcf, out_vcf, filter_bed_file):
    gene_dict = {}
    for line in filter_bed_file:
        columns = line.strip().split("\t")
        chrom = columns[0]
        start = int(columns[1])
        end = int(columns[2])
        gene = columns[3]
        if chrom not in gene_dict:
            gene_dict[chrom] = [[start, end, gene]]
        else:
            gene_dict[chrom].append([start, end, gene])

    vcf_out = open(out_vcf, "a")
    vcf_in = open(in_vcf)
    header = True
    info_ids = {}
    for line in vcf_in:
        if header:
            if line[:6] == "#CHROM":
                vcf_out.write("##INFO=<ID=Genes,Number=1,Type=String,Description=\"Gene names\">\n")
                vcf_out.write(line)
                header = False
            elif line[:11] == "##INFO=<ID=" and line.find("pipe separated list of all details in the") != -1:
                header_id = "##INFO=<ID=CNV_" + line.split(",")[0].split("_")[-1]
                header_id_new = "CNV_" + line.split(",")[0].split("_")[-1]
                header_id_org = line.split(",")[0].split("ID=")[1]
                info_ids[header_id_org] = header_id_new
                for hi in line.split(",")[1:]:
                    header_id += "," + hi
                vcf_out.write(header_id)
            else:
                vcf_out.write(line)
            continue
        columns = line.strip().split("\t")
        chrom = columns[0]
        start = int(columns[1])
        INFO = columns[7]
        end = int(INFO.split("END=")[1].split(";")[0])
        INFO_mod = INFO.split(";")[0]
        for info in INFO.split(";")[1:]:
            if info.split("=")[0] in info_ids:
                info_mod = info_ids[info.split("=")[0]] + "=" + info.split("=")[1]
                INFO_mod = "%s;%s" % (INFO_mod, info_mod)
            else:
                INFO_mod = "%s;%s" % (INFO_mod, info)

        keep_variant, genes = variant_in_genelist(chrom, start, end, gene_dict)
        if keep_variant:
            INFO_mod = "Genes=%s;%s" % (genes, INFO_mod)
            columns[7] = INFO_mod
            vcf_out.write(columns[0])
            for column in columns[1:]:
                vcf_out.write("\t" + column)
            vcf_out.write("\n")
    vcf_out.close()


if __name__ == "__main__":
    in_vcf = snakemake.input.vcf
    out_vcf = snakemake.output.vcf
    filter_bed_file = open(snakemake.params.filter_config)

    filter_variants(in_vcf, out_vcf, filter_bed_file)
