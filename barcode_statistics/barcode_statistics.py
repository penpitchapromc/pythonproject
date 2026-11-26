#!/usr/bin/env python3
"""Calculate sequencing summary statistics grouped by barcode.

The script accepts FASTQ or FASTQ.GZ files, or directories containing them,
and writes one summary row per barcode to a CSV file.
"""

import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path

from read_fastq.read_fastq import read_fastq

FASTQ_SUFFIXES = (".fastq", ".fastq.gz", ".fq", ".fq.gz")


def _read_quality_scores(quality):
    """Convert FASTQ Phred+33 ASCII quality characters to numeric scores."""
    scores = [ord(character) - 33 for character in quality]
    if any(score < 0 for score in scores):
        raise ValueError("Quality string contains a character below Phred+33 range")
    return scores


def _nx_read_length(lengths, fraction):
    """Return the N-fraction read length (for example, N50 or N90)."""
    if not lengths:
        return 0
    target_bases = sum(lengths) * fraction
    cumulative_bases = 0
    for length in sorted(lengths, reverse=True): # Sort the read lengths in descending order
        cumulative_bases += length
        if cumulative_bases >= target_bases: # Return the read length at which the cumulative base count reaches or exceeds the target fraction of total bases
            return length
    return 0


def calculate_barcode_stats(fastq_files):
    """Return summary dictionaries for one or more FASTQ files.

    Quality is summarized in two ways:
    * ``mean_base_quality`` and ``Q*_base_pct`` use all bases as the unit.
    * ``mean_read_quality`` and ``reads_Q*_pct`` use each read's mean Q-score.
    """
    if isinstance(fastq_files, (str, Path)):
        fastq_files = [fastq_files]

    # defaultdict calls the lambda and creates that barcode’s stats.
    # The lambda makes a new dictionary each time, so different barcodes don’t share their lists or counters.
    barcode_stats = defaultdict(
        lambda: {
            "read_lengths": [],
            "read_quality_means": [],
            "total_bases": 0,
            "quality_sum": 0,
            "q10_bases": 0,
            "q15_bases": 0,
            "q20_bases": 0,
            "reads_q10": 0,
            "reads_q15": 0,
            "reads_q20": 0,
        }
    )

    total_reads = 0
    for fastq_file in fastq_files:
        for record in read_fastq(fastq_file):
            total_reads += 1
            barcode = record["barcode"]
            read_length = len(record["seq"])
            q_scores = _read_quality_scores(record["qual"])

            stats = barcode_stats[barcode]
            stats["read_lengths"].append(read_length)
            stats["total_bases"] += read_length
            stats["quality_sum"] += sum(q_scores)
            stats["q10_bases"] += sum(score >= 10 for score in q_scores)
            stats["q15_bases"] += sum(score >= 15 for score in q_scores)
            stats["q20_bases"] += sum(score >= 20 for score in q_scores)

            mean_read_quality = statistics.mean(q_scores) if q_scores else 0.0
            stats["read_quality_means"].append(mean_read_quality)
            stats["reads_q10"] += mean_read_quality >= 10
            stats["reads_q15"] += mean_read_quality >= 15
            stats["reads_q20"] += mean_read_quality >= 20

    results = []
    for barcode, stats in sorted(barcode_stats.items()):
        read_lengths = stats["read_lengths"]
        read_quality_means = stats["read_quality_means"]
        read_count = len(read_lengths)
        total_bases = stats["total_bases"]
        quality_bases = total_bases

        results.append(
            {
                "barcode": barcode,
                "read_count": read_count,
                "read_proportion_pct": read_count / total_reads * 100,
                "total_bases": total_bases,
                "mean_read_length": statistics.mean(read_lengths),
                "median_read_length": statistics.median(read_lengths),
                "min_read_length": min(read_lengths),
                "max_read_length": max(read_lengths),
                "N50_read_length": _nx_read_length(read_lengths, 0.50),
                "N90_read_length": _nx_read_length(read_lengths, 0.90),
                "mean_base_quality": (
                    stats["quality_sum"] / quality_bases if quality_bases else 0.0
                ),
                "mean_read_quality": statistics.mean(read_quality_means),
                "median_read_quality": statistics.median(read_quality_means),
                "Q10_base_pct": (
                    stats["q10_bases"] / quality_bases * 100
                    if quality_bases
                    else 0.0
                ),
                "Q15_base_pct": (
                    stats["q15_bases"] / quality_bases * 100
                    if quality_bases
                    else 0.0
                ),
                "Q20_base_pct": (
                    stats["q20_bases"] / quality_bases * 100
                    if quality_bases
                    else 0.0
                ),
                "reads_Q10_pct": stats["reads_q10"] / read_count * 100,
                "reads_Q15_pct": stats["reads_q15"] / read_count * 100,
                "reads_Q20_pct": stats["reads_q20"] / read_count * 100,
            }
        )

    return results


def discover_fastq_files(inputs):
    """Resolve input files and directories into a sorted list of FASTQ files."""
    files = []
    for value in inputs:
        path = Path(value)
        if path.is_dir():
            files.extend(
                child
                for child in path.rglob("*")
                if child.is_file()
                and any(child.name.lower().endswith(suffix) for suffix in FASTQ_SUFFIXES)
            )
        elif path.is_file():
            files.append(path)
        else:
            raise FileNotFoundError(f"Input does not exist: {path}")

    files = sorted(set(files))
    if not files:
        raise FileNotFoundError("No FASTQ or FASTQ.GZ files found in the input")
    return files


def write_csv(rows, output_path):
    """Write summary rows to a CSV file."""
    if not rows:
        raise ValueError("No reads were found; no summary CSV was written")
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Calculate sequencing statistics grouped by barcode."
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help="FASTQ/FASTQ.GZ files or directories containing FASTQ files",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="barcode_statistics.csv",
        help="output CSV path (default: barcode_statistics.csv)",
    )
    args = parser.parse_args(argv)

    fastq_files = discover_fastq_files(args.inputs)
    rows = calculate_barcode_stats(fastq_files)
    write_csv(rows, args.output)
    read_count = sum(row["read_count"] for row in rows)
    print(
        f"Processed {read_count} reads from {len(fastq_files)} FASTQ file(s); "
        f"wrote {args.output}"
    )


if __name__ == "__main__":
    main()
