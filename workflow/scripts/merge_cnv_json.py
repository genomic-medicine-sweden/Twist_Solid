import cyvcf2
import json


def parse_fai(filename, skip=None):
    with open(filename) as f:
        for line in f:
            chrom, length = line.strip().split()[:2]
            if not skip is None and chrom in skip:
                continue
            yield chrom, int(length)


def parse_annotation_bed(filename, skip=None):
    with open(filename) as f:
        for line in f:
            chrom, start, end, name = line.strip().split()[:4]
            if not skip is None and chrom in skip:
                continue
            yield chrom, int(start), int(end), name


def get_vaf(vcf_filename, skip=None):
    vcf = cyvcf2.VCF(vcf_filename)
    for variant in vcf:
        if not skip is None and variant.CHROM in skip:
            continue
        yield variant.CHROM, variant.POS, variant.INFO.get("AF", None)


def merge_cnv_dicts(dicts, vaf, annotations, chromosomes):
    callers = list(map(lambda x: x["caller"], dicts))
    cnvs = {}
    for chrom, chrom_length in chromosomes:
        cnvs[chrom] = dict(
            chromosome=chrom,
            length=chrom_length,
            vaf=[],
            annotations=[],
            callers={c: dict(name=c, ratios=[], segments=[]) for c in callers},
        )

    for a in annotations:
        cnvs[a[0]]["annotations"].append(
            dict(
                start=a[1],
                end=a[2],
                name=a[3],
            )
        )

    for v in vaf:
        cnvs[v[0]]["vaf"].append(
            dict(
                pos=v[1],
                vaf=v[2],
            )
        )

    for d in dicts:
        for r in d["ratios"]:
            cnvs[r["chromosome"]]["callers"][d["caller"]]["ratios"].append(
                dict(
                    start=r["start"],
                    end=r["end"],
                    log2=r["log2"],
                )
            )
        for s in d["segments"]:
            cnvs[s["chromosome"]]["callers"][d["caller"]]["segments"].append(
                dict(
                    start=s["start"],
                    end=s["end"],
                    log2=s["log2"],
                )
            )

    for v in cnvs.values():
        v["callers"] = list(v["callers"].values())

    return list(cnvs.values())


def main():
    annotation_bed = snakemake.input["annotation_bed"]
    fasta_index_file = snakemake.input["fai"]
    germline_vcf = snakemake.input["germline_vcf"]
    json_files = snakemake.input["json"]

    output_file = snakemake.output["json"]

    skip_chromosomes = snakemake.params["skip_chromosomes"]

    cnv_dicts = []
    for fname in json_files:
        with open(fname) as f:
            cnv_dicts.append(json.load(f))

    fai = parse_fai(fasta_index_file, skip_chromosomes)
    vaf = get_vaf(germline_vcf)
    annotations = []
    if len(annotation_bed) > 0:
        annotations = parse_annotation_bed(annotation_bed, skip_chromosomes)

    cnvs = merge_cnv_dicts(cnv_dicts, vaf, annotations, fai)

    with open(output_file, "w") as f:
        print(json.dumps(cnvs), file=f)


if __name__ == "__main__":
    main()
