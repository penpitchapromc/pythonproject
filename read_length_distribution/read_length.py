
"""Calculate read length statistics from sequencing reads.

This module calculates read length distributions, N50 and N90 read lengths,
the percentage of reads above a specified length, and cumulative sequencing
yield grouped by barcode.
"""

def get_read_lengths(records):
    """Return the sequence length of each FASTQ record."""
    lengths = []
    for record in records:
        lengths.append(len(record["seq"]))
    return lengths

def read_length_distribution(lengths, bin_size):
    """Group read lengths into ranges of a specified bin size."""
    distribution = {}

    if bin_size <= 0:
        raise ValueError("bin_size must be greater than 0")

    for current_length in lengths:
        range_start = (current_length // bin_size) * bin_size
        distribution.setdefault(range_start, 0)
        distribution[range_start] += 1

    return distribution

def n50_read_length(lengths):
    """Return the N50 read length."""

    if not lengths:
        return 0

    sorted_lengths = sorted(lengths, reverse=True)
    half_len = sum(lengths) / 2

    cumulative_length = 0

    for length in sorted_lengths:
        cumulative_length += length

        if cumulative_length >= half_len:
            return length

def n90_read_length(lengths):
    """Return the N90 read length."""

    if not lengths:
        return 0
    
    sorted_lengths = sorted(lengths, reverse=True)
    target_yield = sum(lengths) * 0.9

    cumulative_length = 0

    for current_length in sorted_lengths:
        cumulative_length += current_length

        if cumulative_length >= target_yield:
            return current_length

def percentage_above_length(lengths, cutoff):
    """Return the percentage of reads at or above a specified length."""
    num_above = 0
    for length in lengths:
        if length >= cutoff:
            num_above += 1
    try:
        percentage_above = (num_above / len(lengths)) * 100
    except ZeroDivisionError:
        return 0
    else:
        return percentage_above

def cumulative_yield_per_barcode(records):
    """Return the total sequencing yield in bases for each barcode."""
    barcode_yield = {}

    for record in records:
        barcode = record["barcode"]
        read_length = len(record["seq"])
        barcode_yield.setdefault(barcode, 0)
        barcode_yield[barcode] += read_length
    return barcode_yield
