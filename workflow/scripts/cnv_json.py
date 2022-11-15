from collections import defaultdict
import cyvcf2
import functools
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


def parse_svdb_vcf(vcf_filename):
    vcf = cyvcf2.VCF(vcf_filename)
    cnvs = defaultdict(list)

    for variant in vcf:
        if len(variant.INFO.get("Genes", "")) == 0:
            continue
        genes = variant.INFO["Genes"].split(",")
        if len(genes) == 0:
            continue

        variant_type = variant.ALT

        assert len(variant_type) == 1
        variant_type = variant_type[0].strip("<>")

        caller = variant.INFO.get("CALLER")
        if "gatk" in caller:
            caller = "gatk"

        for g in genes:
            for v in cnvs[g]:
                if variant.CHROM == v["chromosome"] and variant.INFO.get("SVLEN") == v["length"] and caller == v["caller"]:
                    break
            else:
                cnvs[g].append(
                    dict(
                        chromosome=variant.CHROM,
                        start=variant.POS,
                        type=variant_type,
                        length=variant.INFO.get("SVLEN"),
                        caller=caller,
                        cn=variant.INFO.get("CORR_CN"),
                    )
                )

    return cnvs


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


def to_json(
    cnvkit_segments,
    cnvkit_ratios,
    gatk_segments,
    gatk_ratios,
    chroms,
    amp,
    loh,
    vaf,
    cnvs,
    skip=None,
):
    cnvkit_list = []
    for chrom, length in chroms.items():
        if skip is not None and chrom in skip:
            continue
        cnvkit_list.append(
            dict(
                chromosome=chrom,
                label=chrom,
                length=length,
                cnvs=cnvs.get(chrom, []),
                cnvkit_segments=cnvkit_segments.get(chrom, []),
                cnvkit_ratios=cnvkit_ratios.get(chrom, []),
                gatk_segments=gatk_segments.get(chrom, []),
                gatk_ratios=gatk_ratios.get(chrom, []),
                genes=amp.get(chrom, []) + loh.get(chrom, []),
                vaf=vaf.get(chrom, []),
            )
        )

    return json.dumps(cnvkit_list)


def merge_chrom_dicts(d1, d2):
    res = {}
    for k in set(list(d1.keys()) + list(d2.keys())):
        res[k] = d1.get(k, []) + d2.get(k, [])
    return res


def merge_cnvs(unfiltered_cnvs, filtered_cnvs, n_callers=2):
    cnv_union = defaultdict(list)

    for (cnv_dict, unfiltered_cnv_dict) in zip(filtered_cnvs, unfiltered_cnvs):
        for g, cnv_list in cnv_dict.items():
            callers = set([c["caller"] for c in cnv_list])
            for cnv in cnv_list:
                cnv_union[cnv["chromosome"]].append(
                    dict(
                        gene=g,
                        start=cnv["start"],
                        length=cnv["length"],
                        caller=cnv["caller"],
                        type=cnv["type"],
                        cn=cnv["cn"],
                    )
                )
            if len(callers) < n_callers:
                # Results for at least one caller is missing
                for raw_cnv in unfiltered_cnv_dict[g]:
                    if not raw_cnv in cnv_list:
                        cnv_union[raw_cnv["chromosome"]].append(
                            dict(
                                gene=g,
                                start=raw_cnv["start"],
                                length=raw_cnv["length"],
                                caller=raw_cnv["caller"],
                                type=raw_cnv["type"],
                                cn=raw_cnv["cn"],
                            )
                        )

    return cnv_union


def main():
    cnvkit_segments_filename = snakemake.input.cns
    cnvkit_ratios_filename = snakemake.input.cnr
    gatk_segments_filename = snakemake.input.gatk_segments
    gatk_ratios_filename = snakemake.input.gatk_ratios
    svdb_vcfs = snakemake.input.svdb_vcfs
    svdb_filtered_vcfs = snakemake.input.svdb_filtered_vcfs
    vcf_filename = snakemake.input.germline_vcf
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

    unfiltered_cnvs = list(map(parse_svdb_vcf, svdb_vcfs))
    filtered_cnvs = list(map(parse_svdb_vcf, svdb_filtered_vcfs))
    reported_cnvs = merge_cnvs(unfiltered_cnvs, filtered_cnvs)

    cnvkit_json = to_json(
        cnvkit_segments,
        cnvkit_ratios,
        gatk_segments,
        gatk_ratios,
        chroms,
        amp,
        loh,
        vaf,
        reported_cnvs,
        skip=skip_chromosomes,
    )
    with open(json_filename, "w") as f:
        f.write(cnvkit_json)


if __name__ == "__main__":
    main()
