"""quality_stats.py - Module 2: Quality score statistics.

Nanopore FASTQ quality scores are ASCII characters that encode Phred scores
(Phred+33): phred_scores() turns each character into a number with ord(char) - 33.

Calculates
    - mean quality score per read
    - median quality score per read
    - mean quality score per barcode
    - percentage of reads with Q >= 10, 15, 20

Usage (from the repository folder):
    from barcode_qc import read_fastq
    from quality_score_statistics.quality_stats import collect_qualities, quality_stats

    data = collect_qualities(read_fastq("data/sample.fastq"))
    stats = quality_stats(data)
"""

import math


# ---------- small helpers ----------

def mean(values):
    """Average of a list of numbers (0 for an empty list)."""
    if len(values) == 0:
        return 0
    return sum(values) / len(values)


def median(values):
    """Middle value of a list of numbers (0 for an empty list)."""
    if len(values) == 0:
        return 0
    values = sorted(values)
    n = len(values)
    middle = n // 2
    if n % 2 == 1:
        return values[middle]
    return (values[middle - 1] + values[middle]) / 2


def percent(part, whole):
    """part as a percentage of whole (0 if whole is 0)."""
    if whole == 0:
        return 0
    return 100 * part / whole

Q_THRESHOLDS = [10, 15, 20]


def mean_quality(quals):
    """Mean Phred quality of ONE read.

    Phred scores are logarithmic, so we cannot average them directly.
    1. turn each score into an error probability: 10 ** (-Q / 10)
    2. average the error probabilities
    3. turn the average back into a Phred score: -10 * log10(average)
    """
    if len(quals) == 0:
        return 0.0
    total_error = 0.0
    for q in quals:
        total_error = total_error + 10 ** (-q / 10)
    mean_error = total_error / len(quals)
    return -10 * math.log10(mean_error)


def median_quality(quals):
    """Median Phred quality of ONE read."""
    return median(quals)


def summarize_quality(read_means, read_medians, thresholds=Q_THRESHOLDS):
    """Quality statistics for one group of reads (one barcode or all)."""
    result = {
        "mean_q": mean(read_means),
        "median_q": median(read_means),
        "mean_of_read_medians": mean(read_medians),
    }
    for t in thresholds:
        count = 0
        for q in read_means:
            if q >= t:
                count = count + 1
        result[f"pct_q{t}"] = percent(count, len(read_means))
    return result


def phred_scores(quality_string):
    """Turn a FASTQ quality line into numbers (Phred+33).

    '+' -> 43 - 33 = Q10,  '5' -> 53 - 33 = Q20,  '?' -> 63 - 33 = Q30
    """
    scores = []
    for char in quality_string:
        scores.append(ord(char) - 33)
    return scores


def collect_qualities(records):
    """Mean Q and median Q of every read, grouped by barcode.

    records: FASTQ records as dictionaries with "barcode" and "qual",
             e.g. from barcode_qc.read_fastq() (the same records Module 3 uses)
    Returns {barcode: {"quals": [...], "median_quals": [...]}}
    """
    data = {}
    for record in records:
        barcode = record["barcode"]
        scores = phred_scores(record["qual"])
        if barcode not in data:
            data[barcode] = {"quals": [], "median_quals": []}
        data[barcode]["quals"].append(mean_quality(scores))
        data[barcode]["median_quals"].append(median_quality(scores))
    return data


def quality_stats(data, thresholds=Q_THRESHOLDS):
    """Module 2: quality statistics per barcode, plus "all".

    data = {barcode: {"quals": [...], "median_quals": [...]}}  (from collect_qualities)
    """
    stats = {}
    all_means = []
    all_medians = []
    for barcode in sorted(data):
        means = data[barcode]["quals"]
        medians = data[barcode]["median_quals"]
        stats[barcode] = summarize_quality(means, medians, thresholds)
        all_means.extend(means)
        all_medians.extend(medians)
    stats["all"] = summarize_quality(all_means, all_medians, thresholds)
    return stats
