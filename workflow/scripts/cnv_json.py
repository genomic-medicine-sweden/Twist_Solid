from collections import defaultdict
import cyvcf2
import json


def parse_cns(cnvkit_segments_filename):
    cns_dict = defaultdict(list)
    with open(cnvkit_segments_filename) as f:
        # skip header
        next(f)
        for line in f:
            (
                chrom,
                start,
                end,
                gene,
                log2,
                depth,
                probes,
                weight,
                ci_lo,
                ci_hi,
            ) = line.strip().split()
            start = int(start)
            end = int(end)
            log2 = float(log2)

            cns_dict[chrom].append(
                dict(
                    start=start,
                    end=end,
                    log2=log2,
                )
            )

    return cns_dict


def parse_gatk_segments(gatk_segments_filename):
    segment_dict = defaultdict(list)
    with open(gatk_segments_filename) as f:
        for line in f:
            if line.startswith("@") or line.startswith("CONTIG"):
                continue
            (
                chrom,
                start,
                end,
                n_points,
                log2,
            ) = line.strip().split()
            start = int(start)
            end = int(end)
            log2 = float(log2)

            segment_dict[chrom].append(
                dict(
                    start=start,
                    end=end,
                    log2=log2,
                )
            )
    return segment_dict


def get_vaf(vcf_filename):
    vafs = defaultdict(list)
    vcf = cyvcf2.VCF(vcf_filename)
    for variant in vcf:
        vafs[variant.CHROM].append(dict(pos=variant.POS, vaf=variant.INFO.get("AF", None)))
    return vafs


def parse_cnr(cnvkit_ratios_filename):
    cnr_dict = defaultdict(list)
    with open(cnvkit_ratios_filename) as f:
        # skip header
        next(f)
        for line in f:
            chrom, start, end, gene, depth, log2, weight = line.strip().split()
            start = int(start)
            end = int(end)
            log2 = float(log2)

            cnr_dict[chrom].append(
                dict(
                    start=start,
                    end=end,
                    log2=log2,
                )
            )

    return cnr_dict


def parse_gatk_ratios(gatk_ratios_filename):
    ratio_dict = defaultdict(list)
    with open(gatk_ratios_filename) as f:
        for line in f:
            # skip header
            if line.startswith("@") or line.startswith("CONTIG"):
                continue
            (
                chrom,
                start,
                end,
                log2,
            ) = line.strip().split()
            start = int(start)
            end = int(end)
            log2 = float(log2)

            ratio_dict[chrom].append(
                dict(
                    start=start,
                    end=end,
                    log2=log2,
                )
            )

    return ratio_dict


def parse_bed(bed_filename):
    bed = defaultdict(list)
    with open(bed_filename) as f:
        for line in f:
            chrom, start, end, gene = line.strip().split()
            bed[chrom].append(
                dict(
                    start=int(start),
                    end=int(end),
                    gene=gene,
                )
            )
    return bed


def parse_fai(fai_filename):
    chroms = dict()
    with open(fai_filename) as f:
        for line in f:
            line = line.strip().split()
            chroms[line[0]] = int(line[1])
    return chroms


def to_json(cnvkit_segments, cnvkit_ratios, gatk_segments, gatk_ratios, chroms, amp, loh, vaf, skip=None):
    cnvkit_list = []
    for chrom, length in chroms.items():
        if skip is not None and chrom in skip:
            continue
        cnvkit_list.append(
            dict(
                chromosome=chrom,
                label=chrom,
                length=length,
                cnvkit_segments=cnvkit_segments.get(chrom, []),
                cnvkit_ratios=cnvkit_ratios.get(chrom, []),
                gatk_segments=gatk_segments.get(chrom, []),
                gatk_ratios=gatk_ratios.get(chrom, []),
                genes=amp.get(chrom, []) + loh.get(chrom, []),
                vaf=vaf.get(chrom, []),
            )
        )

    return json.dumps(cnvkit_list)


def main():
    cnvkit_segments_filename = snakemake.input.cns
    cnvkit_ratios_filename = snakemake.input.cnr
    gatk_segments_filename = snakemake.input.gatk_segments
    gatk_ratios_filename = snakemake.input.gatk_ratios
    vcf_filename = snakemake.input.vcf
    fai_filename = snakemake.input.fai
    amp_filename = snakemake.input.amp_bed
    loh_filename = snakemake.input.loh_bed
    json_filename = snakemake.output.json

    skip_chromosomes = snakemake.params.skip_chromosomes

    cnvkit_segments = parse_cns(cnvkit_segments_filename)
    cnvkit_ratios = parse_cnr(cnvkit_ratios_filename)
    gatk_segments = parse_gatk_segments(gatk_segments_filename)
    gatk_ratios = parse_gatk_ratios(gatk_ratios_filename)
    chroms = parse_fai(fai_filename)
    amp = {}
    if len(amp_filename) > 0:
        amp = parse_bed(amp_filename)
    loh = {}
    if len(loh_filename) > 0:
        loh = parse_bed(loh_filename)
    vaf = get_vaf(vcf_filename)

    cnvkit_json = to_json(
        cnvkit_segments,
        cnvkit_ratios,
        gatk_segments,
        gatk_ratios,
        chroms,
        amp,
        loh,
        vaf,
        skip=skip_chromosomes,
    )
    with open(json_filename, "w") as f:
        f.write(cnvkit_json)


if __name__ == "__main__":
    main()
