from collections import defaultdict
import cyvcf2
import json


def parse_cns(cns_filename):
    cns_dict = defaultdict(list)
    with open(cns_filename) as f:
        # skip header
        next(f)
        for line in f:
            chrom, start, end, gene, log2, depth, probes, weight, ci_lo, ci_hi = line.strip().split()
            start = int(start)
            end = int(end)
            genes = list(set(map(lambda x: x.split("_")[0], gene.strip().split(","))))
            log2 = float(log2)
            depth = float(depth)
            probes = int(probes)
            weight = float(weight)
            ci_lo = float(ci_lo)
            ci_hi = float(ci_hi)

            cns_dict[chrom].append(dict(
                start=start,
                end=end,
                genes=genes,
                depth=depth,
                log2=log2,
                weight=weight,
                ci_lo=ci_lo,
                ci_hi=ci_hi
            ))

    return cns_dict


def get_vaf(vcf_filename):
    vafs = defaultdict(list)
    vcf = cyvcf2.VCF(vcf_filename)
    for variant in vcf:
        vafs[variant.CHROM].append(dict(
            pos=variant.POS,
            vaf=variant.INFO.get("AF", None)
        ))
    return vafs


def parse_cnr(cnr_filename):
    cnr_dict = defaultdict(list)
    with open(cnr_filename) as f:
        # skip header
        next(f)
        for line in f:
            chrom, start, end, gene, depth, log2, weight = line.strip().split()
            start = int(start)
            end = int(end)
            depth = float(depth)
            log2 = float(log2)
            weight = float(weight)

            cnr_dict[chrom].append(dict(
                gene=gene,
                start=start,
                end=end,
                depth=depth,
                log2=log2,
                weight=weight
            ))

    return cnr_dict


def parse_bed(bed_filename):
    bed = defaultdict(list)
    with open(bed_filename) as f:
        for line in f:
            chrom, start, end, gene = line.strip().split()
            bed[chrom].append(dict(
                start=int(start),
                end=int(end),
                gene=gene,
            ))
    return bed


def parse_fai(fai_filename):
    chroms = dict()
    with open(fai_filename) as f:
        for line in f:
            line = line.strip().split()
            chroms[line[0]] = int(line[1])
    return chroms


def to_json(cns, cnr, chroms, amp, loh, vaf, skip=None):
    cnvkit_list = []
    for chrom, length in chroms.items():
        if not skip is None and chrom in skip:
            continue
        cnvkit_list.append(dict(
            chromosome=chrom,
            label=chrom,
            length=length,
            segments=cns.get(chrom, []),
            regions=cnr.get(chrom, []),
            genes=amp.get(chrom, []) + loh.get(chrom, []),
            vaf=vaf.get(chrom, []),
        ))

    return json.dumps(cnvkit_list)


def main():
    cns_filename = snakemake.input.cns
    cnr_filename = snakemake.input.cnr
    vcf_filename = snakemake.input.vcf
    fai_filename = snakemake.input.fai
    amp_filename = snakemake.input.amp_bed
    loh_filename = snakemake.input.loh_bed
    json_filename = snakemake.output.json

    skip_chromosomes = snakemake.params.skip_chromosomes

    cns = parse_cns(cns_filename)
    cnr = parse_cnr(cnr_filename)
    chroms = parse_fai(fai_filename)
    amp = {}
    if len(amp_filename) > 0:
        amp = parse_bed(amp_filename)
    loh = {}
    if len(loh_filename) > 0:
        loh = parse_bed(loh_filename)
    vaf = get_vaf(vcf_filename)

    cnvkit_json = to_json(cns, cnr, chroms, amp, loh, vaf, skip=skip_chromosomes)
    with open(json_filename, "w") as f:
        f.write(cnvkit_json)


if __name__ == "__main__":
    main()
