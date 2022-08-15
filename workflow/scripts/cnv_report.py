import logging
from pysam import VariantFile
from hydra_genetics.utils.io import utils

log = logging.getLogger()


def create_tsv_report(input_vcfs, input_org_vcfs, output_txt):
    gene_all_dict = {}
    for input_org_vcf in input_org_vcfs:
        log.info(f"Opening vcf file: {input_org_vcf}")
        variants = VariantFile(input_org_vcf)
        samples = list(variants.header.samples)
        if len(samples) > 1:
            raise Exeception(f"Unable to process vcf with more then one sample: {samples}")
        else:
            samples = samples[0]
        for variant in variants:
            genes = utils.get_annotation_data_info(variant, "Genes")
            log.debug(f"Processing variant: {variant}")
            if isinstance(genes, tuple):
                genes = ",".join(genes)
            chr = variant.chrom
            start = variant.pos
            end = variant.pos + int(utils.get_annotation_data_info(variant, "SVLEN")) - 1
            callers = utils.get_annotation_data_info(variant, "CALLER")
            cn = utils.get_annotation_data_info(variant, "CORR_CN")
            AF = utils.get_annotation_data_info(variant, "Twist_AF")
            if AF is None:
                AF = 0.0
            if genes is not None:
                for gene in genes.split(","):
                    if gene not in gene_all_dict:
                        gene_all_dict[gene] = []
                        gene_all_dict[gene].append([chr, start, end, callers, cn, AF])
                    else:
                        duplicate = False
                        for variant in gene_all_dict[gene]:
                            if chr == variant[0] and start == variant[1] and end == variant[2] and callers == variant[3]:
                                duplicate = True
                                break
                        if not duplicate:
                            gene_all_dict[gene].append([chr, start, end, callers, cn, AF])

    first_vcf = True
    for input_vcf in input_vcfs:
        gene_variant_dict = {}
        log.info(f"Opening vcf file: {input_vcf}")
        variants = VariantFile(input_vcf)
        samples = list(variants.header.samples)
        if len(samples) > 1:
            raise Exeception(f"Unable to process vcf with more then one sample: {samples}")
        else:
            samples = samples[0]
        log.info(f"Opening output tsv file: {output_txt}")
        counter = 0
        output_mode = "a"
        if first_vcf:
            output_mode = "w"
        with open(output_txt, output_mode) as writer:
            if first_vcf:
                writer.write("sample\tgene(s)\tchrom\tregion\tcallers\tcopy_number")
            for variant in variants:
                genes = utils.get_annotation_data_info(variant, "Genes")
                log.debug(f"Processing variant: {variant}")
                if isinstance(genes, tuple):
                    genes = ",".join(genes)
                chr = variant.chrom
                start = variant.pos
                end = variant.pos + int(utils.get_annotation_data_info(variant, "SVLEN")) - 1
                callers = utils.get_annotation_data_info(variant, "CALLER")
                cn = utils.get_annotation_data_info(variant, "CORR_CN")
                AF = utils.get_annotation_data_info(variant, "Twist_AF")
                if AF is None:
                    AF = 0.0
                writer.write(f"\n{samples}\t{genes}\t{chr}\t{start}-{end}\t{callers}\t{AF:.2f}\t{cn:.2f}")
                counter += 1
                for gene in genes.split(","):
                    if gene not in gene_variant_dict:
                        gene_variant_dict[gene] = []
                    gene_variant_dict[gene].append([chr, start, end, callers, cn,AF])
            for gene in gene_variant_dict:
                if len(gene_variant_dict[gene]) == 1:
                    caller = gene_variant_dict[gene][0][3]
                    for cnv in gene_all_dict[gene]:
                        if cnv[3] != caller:
                            chr = cnv[0]
                            start = cnv[1]
                            end = cnv[2]
                            callers = cnv[3]
                            cn = cnv[4]
                            AF = cnv[5]
                            if (
                                (start >= gene_variant_dict[gene][0][1] and start <= gene_variant_dict[gene][0][2]) or
                                (end >= gene_variant_dict[gene][0][1] and end <= gene_variant_dict[gene][0][2]) or
                                (gene_variant_dict[gene][0][1] >= start and gene_variant_dict[gene][0][1] <= end) or
                                (gene_variant_dict[gene][0][2] >= end and gene_variant_dict[gene][0][2] <= start)
                            ):
                                writer.write(f"\n{samples}\t{gene}\t{chr}\t{start}-{end}\t{callers}\t{AF:.2f}\t{cn:.2f}")
        log.info(f"Processed {counter} variants")
        first_vcf = False


if __name__ == "__main__":
    in_vcfs = snakemake.input.vcfs
    in_org_vcfs = snakemake.input.org_vcfs
    out_tsv = snakemake.output.tsv
    create_tsv_report(in_vcfs, in_org_vcfs, out_tsv)
