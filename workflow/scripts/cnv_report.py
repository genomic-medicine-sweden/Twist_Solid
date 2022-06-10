import logging
from pysam import VariantFile
from hydra_genetics.utils.io import utils

log = logging.getLogger()


def create_tsv_report(input_vcf, output_txt):
    log.info(f"Opening vcf file: {input_vcf}")
    variants = VariantFile(input_vcf)
    samples = list(variants.header.samples)
    if len(samples) > 1:
        raise Exeception(f"Unable to process vcf with more then one sample: {samples}")
    else:
        samples = samples[0]
    log.info(f"Opening output tsv file: {output_txt}")
    counter = 0
    with open(output_txt, "w") as writer:
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
            writer.write(f"\n{samples}\t{genes}\t{chr}\t{start}-{end}\t{callers}\t{cn:.2f}")
            counter += 1
    log.info(f"Processed {counter} variants")


if __name__ == "__main__":
    in_vcf = snakemake.input.vcf
    out_tsv = snakemake.output.tsv
    create_tsv_report(in_vcf, out_tsv)
