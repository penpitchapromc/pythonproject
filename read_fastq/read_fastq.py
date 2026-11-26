#!/usr/bin/env python3
""" Read records from a FASTQ file.
and yield them as dictionaries with name, seq, qual, and barcode fields.
Reads without a barcode=... header tag are labeled unclassified.
"""


import argparse
import gzip
import re
from pathlib import Path
from typing import Counter


BARCODE_PATTERN = re.compile(r"barcode=(\S+)")

FASTQ_SUFFIXES = (".fastq", ".fastq.gz", ".fq", ".fq.gz")


def check_record(header, sequence, plus, quality):
    # Check the structure and sequence/quality lengths of a FASTQ record.
    if not header.startswith("@"):
        raise ValueError("FASTQ header does not start with '@': " + header[:30])
    if not plus.startswith("+"):
        raise ValueError("FASTQ separator line does not start with '+': " + header[:30])
    if len(sequence) != len(quality):
        raise ValueError(
            "Sequence and quality lengths differ for "
            f"{header[:30]} ({len(sequence)} != {len(quality)})"
        )



def open_fastq(path):
    # Open a plain-text or gzipped FASTQ file as ASCII text."""
    path = Path(path)
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="ascii") # "rt" = read text
    return path.open("rt", encoding="ascii")



def read_fastq(path):
    """Yield records from a standard four-line FASTQ file.

    Each yielded dictionary has ``name``, ``seq``, ``qual``, and ``barcode``
    fields. Reads without a ``barcode=...`` header tag are labeled
    ``unclassified``.
    """
    with open_fastq(path) as file:
        record_number = 0
        while True:
            header_line = file.readline() # Read the header line of the FASTQ record
            if header_line == "":
                break
            record_number += 1
            header = header_line.rstrip("\r\n") # Remove trailing newline characters

            sequence_line = file.readline() # Read the sequence line of the FASTQ record
            plus_line = file.readline() # Read the plus line of the FASTQ record
            quality_line = file.readline() # Read the quality line of the FASTQ record
            if not sequence_line or not plus_line or not quality_line:
                raise ValueError(
                    f"Truncated FASTQ record {record_number} in {path}"
                )

            sequence = sequence_line.rstrip("\r\n")
            plus = plus_line.rstrip("\r\n")
            quality = quality_line.rstrip("\r\n")

            check_record(header, sequence, plus, quality) # Validate the FASTQ record structure and lengths

            # Extract barcode
            barcode_match = BARCODE_PATTERN.search(header)
            barcode = barcode_match.group(1) if barcode_match else "unclassified" # If barcode is not found, set to "unclassified"

            yield {
                "name": header.split()[0][1:],
                "seq": sequence,
                "qual": quality,
                "barcode": barcode,
            }


def main():
    import argparse
    from collections import Counter

    barcode_counts = Counter()

    parser = argparse.ArgumentParser(
        description="Read records from one or more FASTQ files."
    )
    parser.add_argument(
        "fastq_files",
        nargs="+",
        help="Path(s) to FASTQ or FASTQ.GZ files.",
    )
    args = parser.parse_args()

    for fastq_file in args.fastq_files:
        print(f"\nFile: {fastq_file}")
        for record in read_fastq(fastq_file):
            print(record["name"])
            print(record["barcode"])
            print(len(record["seq"]))

            barcode = record["barcode"] or "unclassified"

            barcode_counts[barcode] += 1

    print(barcode_counts)

if __name__ == "__main__":
    main()
