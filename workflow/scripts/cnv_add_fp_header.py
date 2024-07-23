from pysam import VariantFile


def add_fp_header(in_vcf_filename, out_vcf_filename):
    variants = VariantFile(in_vcf_filename)
    header = variants.header
    header.add_meta(
        'INFO', items=[('ID', "FP_FLAG"), ('Number', "1"), ('Type', 'String'), ('Description', 'CNV false positive flag')]
    )
    out_vcf = VariantFile(out_vcf_filename, "w", header=header)
    for variant in variants.fetch():
        out_vcf.write(variant)


if __name__ == "__main__":
    in_vcf_filename = snakemake.input.vcf
    out_vcf_filename = snakemake.output.vcf
    add_fp_header(in_vcf_filename, out_vcf_filename)
