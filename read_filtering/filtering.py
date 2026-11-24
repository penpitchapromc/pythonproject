"""Filter sequencing reads using user-defined thresholds.

This module filters FASTQ records based on read length and quality.
"""

def filter_min_length(records, min_length):
    """Return FASTQ records with read length at or above the minimum threshold."""
    filtered_records = []

    for record in records:
        read_length = len(record["seq"])

        if read_length >= min_length:
            filtered_records.append(record)

    return filtered_records

def filter_max_length(records, max_length):
    """Return FASTQ records with read length at or below the maximum threshold."""
    filtered_records = []

    for record in records:
        read_length = len(record["seq"])

        if read_length <= max_length:
            filtered_records.append(record)

    return filtered_records

def average_qscore(qual):
    """Calculate the average Phred+33 quality score of a read."""
    scores = []

    for letter in qual:
        phred = ord(letter) - 33
        scores.append(phred)
    
    try:
        mean_phred = sum(scores) / len(scores)
    except ZeroDivisionError:
        return 0
    else:
        return mean_phred

def filter_min_qscore(records, min_qscore):
    """Return FASTQ records with average Q-score at or above the minimum threshold."""
    filtered_records = []

    for record in records:
        mean_qscore = average_qscore(record["qual"])

        if mean_qscore >= min_qscore:
            filtered_records.append(record)

    return filtered_records

def filter_reads(records, min_length, max_length, min_qscore):
    """Return FASTQ records that pass all length and average Q-score thresholds."""
    filtered_records = []

    for record in records:
        read_length = len(record["seq"])
        mean_qscore = average_qscore(record["qual"])

        if min_length <= read_length <= max_length and mean_qscore >= min_qscore:
            filtered_records.append(record)
    return filtered_records

