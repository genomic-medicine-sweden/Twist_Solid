import collections
import csv
import functools
import json
import sys


# The functions `parse_*_ratios` functions take a filename of a file containing
# copy number log2-ratios across the genome for a specific CNV caller. The
# functions `parse_*_segments` takes a filename of a file containing log2-ratio
# segments across the genome for a specific caller. The return value from both
# of these functions should be a list of dictionaries, where the dictionaries
# look like
#
# {
#     "chromosome": str,
#     "start": int,
#     "end": int,
#     "log2": float,
# }
#


PARSERS = collections.defaultdict(dict)


def cnv_parser(file_format, header=True, skip=0, comment="#"):
    """
    Decorator for parsers of CNV result files. The first argument of
    the wrapped function should be a path to a file, and this argument
    is replaced with the contents of that file. How the content is represented
    depends on the file format:

    - tsv, csv: a generator over lines, each line being a list of values
    - vcf:

    arguments:
        file_format     a file path
        skip            the number of lines to skip before starting to read
                        the file
    """

    def decorator_cnv_parser(func):
        caller_filetype = func.__name__.split("_")[1:]
        caller = "_".join(caller_filetype[:-1])
        filetype = caller_filetype[-1]

        def line_generator(file, delim):
            for _ in range(skip):
                next(file)
            found_header = False
            for line in csv.reader(file, delimiter=delim):
                if line[0].strip()[0] == comment:
                    continue
                if header and not found_header:
                    found_header = True
                    continue
                yield line

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            f = None
            try:
                f = open(args[0], "r")
                if file_format == "tsv":
                    lines = line_generator(f, delim="\t")
                elif file_format == "csv":
                    lines = line_generator(f, delim=",")
                else:
                    raise IOError(f"invalid filetype: {file_format}")
                args = [lines] + list(args[1:])
                return func(*args, **kwargs)
            finally:
                if f is not None:
                    f.close()

        PARSERS[caller][filetype] = wrapper

    return decorator_cnv_parser


@cnv_parser("tsv", header=True)
def parse_cnvkit_ratios(file):
    ratios = []
    for line in file:
        ratios.append(
            dict(
                chromosome=line[0],
                start=int(line[1]),
                end=int(line[2]),
                log2=float(line[5]),
            )
        )
    return ratios


@cnv_parser("tsv", header=True)
def parse_cnvkit_segments(file):
    segments = []
    for line in file:
        segments.append(
            dict(
                chromosome=line[0],
                start=int(line[1]),
                end=int(line[2]),
                log2=float(line[4]),
            )
        )
    return segments


@cnv_parser("tsv", header=True, comment="@")
def parse_gatk_ratios(file):
    ratios = []
    for line in file:
        ratios.append(
            dict(
                chromosome=line[0],
                start=int(line[1]),
                end=int(line[2]),
                log2=float(line[3]),
            )
        )
    return ratios


@cnv_parser("tsv", header=True, comment="@")
def parse_gatk_segments(file):
    segments = []
    for line in file:
        segments.append(
            dict(
                chromosome=line[0],
                start=int(line[1]),
                end=int(line[2]),
                log2=float(line[4]),
            )
        )
    return segments


def to_json(caller, ratios, segments):
    json_dict = dict(
        caller=caller,
        ratios=ratios,
        segments=segments,
    )
    return json.dumps(json_dict)


def main():
    caller = snakemake.wildcards["caller"]
    ratio_filename = snakemake.input["ratios"]
    segment_filename = snakemake.input["segments"]

    output_filename = snakemake.output["json"]

    if caller not in PARSERS:
        print(f"error: no parser for {caller} implemented", file=sys.stderr)
        sys.exit(1)

    ratios = PARSERS[caller]["ratios"](ratio_filename)
    segments = PARSERS[caller]["segments"](segment_filename)

    with open(output_filename, "w") as f:
        print(to_json(caller, ratios, segments), file=f)


if __name__ == "__main__":
    main()
