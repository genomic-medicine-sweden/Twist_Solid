
import statistics


def read_cnv_data(cnv_file_name, sample_name, region):
    probe_data = []
    gene_probe_index = []
    cnv_file = open(cnv_file_name)
    header = True
    i = 0
    for line in cnv_file:
        if header:
            if line.find("CONTIG") != -1:
                header = False
            continue
        columns = line.strip().split("\t")
        chrom = columns[0]
        start = int(columns[1])
        end = int(columns[2])
        log2_copy_ratio = float(columns[3])
        if chrom == region[0] and ((start >= region[1] and start <= region[2]) or (end >= region[1] and end <= region[2])):
            probe_data.append(log2_copy_ratio)
            if chrom == region[0] and ((start >= region[3] and start <= region[4]) or (end >= region[3] and end <= region[4])):
                gene_probe_index.append(i)
            i += 1
    return probe_data, gene_probe_index


def find_max_probe_diff(probe_data):
    probe_len = len(probe_data)
    i = 0
    j = i + 1 + window_size
    probe_windows1 = []
    probe_windows2 = []
    while j + window_size + 1 < probe_len:
        probe_windows1.append([sum(probe_data[i:i + window_size]) / window_size, i, i + window_size])
        probe_windows2.append([sum(probe_data[j:j + window_size]) / window_size, j, j + window_size])
        i += 1
        j += 1
    k = 0
    max_probe_window_diff = 0.0
    min_probe_window_diff = 0.0
    k_max_index = 0
    k_min_index = 0
    for probe_window in probe_windows1:
        probe_window_diff = probe_window[0] - probe_windows2[k][0]
        if probe_window_diff > max_probe_window_diff:
            max_probe_window_diff = probe_window_diff
            k_max_index = k
        if probe_window_diff < min_probe_window_diff:
            min_probe_window_diff = probe_window_diff
            k_min_index = k
        k += 1
    return k_max_index, k_min_index


def filter_deletions(max_probe_diff_index, min_probe_diff_index, probe_data, gene_probe_index, region, deletions):
    deletions.write(
        "Gene(s)\tLog2_ratio_diff\tMedian_L2R_deletion\tMedian_L2R_surrounding\tNumber_of_data_points\tNumber_of_stdev\n"
    )
    # Filter amplifications
    if max_probe_diff_index >= min_probe_diff_index:
        return
    # Filter large deletions (captured by other tools)
    if min_probe_diff_index - max_probe_diff_index >= region_max_size:
        return
    # Filter deletions without start or end breakpoint
    probe_len = len(probe_data)
    if max_probe_diff_index == 0 or min_probe_diff_index + window_size + 1 == probe_len:
        return
    # Filter deletions not in the actual gene of interest
    in_gene = False
    k = max_probe_diff_index + window_size + 1
    while k < min_probe_diff_index + window_size:
        if k in gene_probe_index:
            in_gene = True
            break
        k += 1
    if not in_gene:
        return
    # Filter too short deletions
    low_probes = probe_data[max_probe_diff_index + window_size + 1:min_probe_diff_index + window_size]
    if len(low_probes) < window_size:
        return
    # Filter deletions with too few probes outside outside the deletion
    high_probes = probe_data[:max_probe_diff_index + window_size] + probe_data[min_probe_diff_index + window_size + 1:]
    if len(high_probes) < 4:
        return
    # Calculate high probes window averages
    high_probes_window_averages = []
    k = 0
    while k + window_size < max_probe_diff_index + window_size:
        high_probes_window_averages.append(sum(probe_data[k:k+window_size]) / window_size)
        k += 1
    k = min_probe_diff_index + window_size + 2
    while k + window_size < len(probe_data):
        high_probes_window_averages.append(sum(probe_data[k:k+window_size]) / window_size)
        k += 1
    # Calculate high and low medians and stdev for high probes window averages
    median_low = statistics.median(low_probes)
    median_high = statistics.median(high_probes_window_averages)
    stdev_high = statistics.stdev(high_probes_window_averages)
    # Filter based on difference in log odds ration between high and low median (min_log_odds_diff)
    if median_high - median_low < min_log_odds_diff:
        return
    # Filter based on the number of standard deviation between the high and low median (min_nr_stdev_diff)
    if median_low > median_high - stdev_high * min_nr_stdev_diff:
        return
    median_diff = median_low - median_high
    nr_std_diff = median_diff / stdev_high
    deletions.write(f"{region[5]}\t{median_diff}\t{median_low}\t{median_high}\t{len(low_probes)}\t{nr_std_diff}\n")


def read_regions_data(regions_file):
    next(regions_file)
    for line in regions_file:
        columns = line.strip().split("\t")
        gene_name = columns[0]
        chrom = columns[1]
        start_cnv_region = int(columns[2])
        stop_cnv_region = int(columns[3])
        gene_start = int(columns[4])
        gene_stop = int(columns[5])
        regions.append([chrom, start_cnv_region, stop_cnv_region, gene_start, gene_stop, gene_name])
    return regions


def call_small_cnv_deletions(
    cnv_file_name, regions_file, deletions, window_size, region_max_size, min_nr_stdev_diff, min_log_odds_diff,
):
    regions = read_regions_data(regions_file)
    sample_name = cnv_file_name.split("/")[1].split("_")[0]
    for region in regions:
        probe_data, gene_probe_index = read_cnv_data(cnv_file_name, sample_name, region)
        max_probe_diff_index, min_probe_diff_index = find_max_probe_diff(probe_data)
        filter_deletions(max_probe_diff_index, min_probe_diff_index, probe_data, gene_probe_index, region, deletions)


if __name__ == "__main__":
    log = snakemake.log_fmt_shell(stdout=False, stderr=True)

    call_small_cnv_deletions(
        snakemake.input.cnv_data,
        open(snakemake.input.regions_file),
        open(snakemake.output.deletions, "w"),
        snakemake.params.window_size,
        snakemake.params.region_max_size,
        snakemake.params.min_nr_stdev_diff,
        snakemake.params.min_log_odds_diff,
    )
