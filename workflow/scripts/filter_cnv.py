
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
        chrom = colums[0]
        start = int(colums[1])
        end = int(colums[2])
        gene = colums[3]
        if chrom not in gene_dict:
            gene_dict[chrom] = [[start, end, gene]]
        else:
            gene_dict[chrom].append([start, end, gene])

    vcf_in = VariantFile(in_vcf)
    new_header = in_vcf.header
    new_header.info.add("Gene", "1", "String", "Gene name")
    out_vcf = VariantFile(out_vcf_filename, 'w', header=new_header)
    out_vcf.write(vcf_in.header)

    for variant in vcf_in:
        chrom = variant.contig
        start = int(variant.pos)
        end = int(variant.info.split("END=")[1].split(";")[0])
        keep_variant, genes = variant_in_genelist(chrom, start, end, gene_dict)
        if keep_variant:
            variant.info.__setitem__('Gene',genes)
            out_vcf.write(variant)

    out_vcf.close()


if __name__ == "__main__":
    in_vcf = snakemake.input.vcf
    out_vcf = open(snakemake.output.vcf, "w")
    filter_gene_file = snakemake.params.filter_config

    filter_variants(in_vcf, out_vcf, filter_gene_file)
