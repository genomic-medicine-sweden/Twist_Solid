import logging
from pysam import VariantFile
from hydra_genetics.utils.io import utils

log = logging.getLogger()


def create_tsv_report(
    input_vcfs, input_org_vcfs, input_del, input_amp, amp_cn_limit,
    output_txt, out_additional_only, del_1p19q_cn, del_1p19q_chr_arm_fraction, TC
):
    gene_all_dict = {}
    nr_writes = 0
    log.info(f"Opening output tsv file: {output_txt}")
    sample_name = ""
    with open(output_txt, "w") as writer:
        writer.write("sample\tgene(s)\tchrom\tregion\tcallers\tfreq_in_db\tcopy_number")
        out_additional_only.write("sample\tgene(s)\tchrom\tregion\tcallers\tfreq_in_db\tcopy_number")
        for input_org_vcf in input_org_vcfs:
            del_1p19q = {
                "1p_cnvkit": 0, "19q_cnvkit": 0, "1p_gatkcnv": 0, "19q_gatkcnv": 0,
                "1p": [0, 125000000, 125000000], "19q": [26500000, 59128983, 32628983],
            }
            log.info(f"Opening vcf file: {input_org_vcf}")
            variants = VariantFile(input_org_vcf)
            samples = list(variants.header.samples)
            if len(samples) > 1:
                raise Exception(f"Unable to process vcf with more then one sample: {samples}")
            else:
                samples = samples[0]
                sample_name = samples
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
                if cn < del_1p19q_cn and chr == "chr1" and start >= del_1p19q["1p"][0] and start <= del_1p19q["1p"][1]:
                    if callers == "cnvkit":
                        del_1p19q["1p_cnvkit"] += end - start + 1
                    elif callers == "gatk":
                        del_1p19q["1p_gatkcnv"] += end - start + 1
                if cn < del_1p19q_cn and chr == "chr19" and start >= del_1p19q["19q"][0] and start <= del_1p19q["19q"][1]:
                    if callers == "cnvkit":
                        del_1p19q["19q_cnvkit"] += end - start + 1
                    elif callers == "gatk":
                        del_1p19q["19q_gatkcnv"] += end - start + 1
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
            if (del_1p19q["1p_cnvkit"] / del_1p19q["1p"][2] > del_1p19q_chr_arm_fraction and
                    del_1p19q["19q_cnvkit"] / del_1p19q["19q"][2] > del_1p19q_chr_arm_fraction):
                if nr_writes < 2:
                    writer.write(f"\n{samples}\t1p19q\tNA\tNA\tcnvkit\tNA\tNA")
                    out_additional_only.write(f"\n{samples}\t1p19q\tNA\tNA\tcnvkit\tNA\tNA")
                    nr_writes += 1
            if (del_1p19q["1p_gatkcnv"] / del_1p19q["1p"][2] > del_1p19q_chr_arm_fraction and
                    del_1p19q["19q_gatkcnv"] / del_1p19q["19q"][2] > del_1p19q_chr_arm_fraction):
                if nr_writes < 2:
                    writer.write(f"\n{samples}\t1p19q\tNA\tNA\tgatk_cnv\tNA\tNA")
                    out_additional_only.write(f"\n{samples}\t1p19q\tNA\tNA\tgatk_cnv\tNA\tNA")
                    nr_writes += 1

        for input_vcf in input_vcfs:
            gene_variant_dict = {}
            log.info(f"Opening vcf file: {input_vcf}")
            variants = VariantFile(input_vcf)
            samples = list(variants.header.samples)
            if len(samples) > 1:
                raise Exception(f"Unable to process vcf with more then one sample: {samples}")
            else:
                samples = samples[0]
            counter = 0
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
                    gene_variant_dict[gene].append([chr, start, end, callers, cn, AF])
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

        deletions = open(input_del)
        header_list = next(deletions).split("\t")
        for deletion in deletions:
            columns = {k: v for k, v in zip(header_list, deletion.strip().split("\t"))}
            gene = columns['Gene(s)']
            chr = columns['Chromosome']
            start = columns['Gene_start']
            end = columns['Gene_end']
            callers = "small_deletion"
            AF = "NA"
            log_odds_ratio = float(columns['Median_L2R_deletion'])
            cn = 2*pow(2, float(log_odds_ratio))
            ccn = cn
            if TC > 0.0:
                ccn = round(2 + (cn - 2) * (1/float(TC)), 2)
            writer.write(f"\n{sample_name}\t{gene}\t{chr}\t{start}-{end}\t{callers}\t{AF}\t{ccn:.2f}")
            out_additional_only.write(f"\n{sample_name}\t{gene}\t{chr}\t{start}-{end}\t{callers}\t{AF}\t{ccn:.2f}")
        amplifications = open(input_amp)
        header_list = next(amplifications).split("\t")
        for amplification in amplifications:
            columns = {k: v for k, v in zip(header_list, amplification.strip().split("\t"))}
            gene = columns['Gene(s)']
            chr = columns['Chromosome']
            start = columns['Gene_start']
            end = columns['Gene_end']
            callers = "small_amplification"
            AF = "NA"
            log_odds_ratio = float(columns['Median_L2R_amplification'])
            cn = 2*pow(2, float(log_odds_ratio))
            ccn = cn
            if TC > 0.0:
                ccn = round(2 + (cn - 2) * (1/float(TC)), 2)
            if ccn > amp_cn_limit:
                writer.write(f"\n{sample_name}\t{gene}\t{chr}\t{start}-{end}\t{callers}\t{AF}\t{ccn:.2f}")
                out_additional_only.write(f"\n{sample_name}\t{gene}\t{chr}\t{start}-{end}\t{callers}\t{AF}\t{ccn:.2f}")


if __name__ == "__main__":
    in_vcfs = snakemake.input.vcfs
    in_org_vcfs = snakemake.input.org_vcfs
    in_del = snakemake.input.deletions
    in_amp = snakemake.input.amplifications
    amp_cn_limit = snakemake.params.call_small_amplifications_cn_limit
    out_tsv = snakemake.output.tsv
    del_1p19q_cn = snakemake.params.del_1p19q_cn_limit
    del_1p19q_chr_arm_fraction = snakemake.params.del_1p19q_chr_arm_fraction
    TC = float(snakemake.params.tc)
    with open(snakemake.output.tsv_additional_only, "w") as out_additional_only:
        create_tsv_report(
            in_vcfs,
            in_org_vcfs,
            in_del,
            in_amp,
            amp_cn_limit,
            out_tsv,
            out_additional_only,
            del_1p19q_cn,
            del_1p19q_chr_arm_fraction,
            TC,
        )
