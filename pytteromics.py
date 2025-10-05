TRANSCR_RULES = {
    "A": "A",
    "T": "U",
    "U": "U",  # in case of virus RNA(-) -> RNA(+) transcription
    "G": "G",
    "C": "C",
    "a": "a",
    "t": "u",
    "u": "u",  # in case of virus RNA(-) -> RNA(+) transcription
    "g": "g",
    "c": "c",
}

DNA_COMPL_RULES = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
    "a": "t",
    "t": "a",
    "g": "c",
    "c": "g",
}

RNA_COMPL_RULES = {
    "A": "U",
    "U": "A",
    "G": "C",
    "C": "G",
    "a": "u",
    "u": "a",
    "g": "c",
    "c": "g",
}


def filter_fastq(seqs, gc_bounds = (0, 100), length_bounds = (0, 2**32), quality_threshold = 0):
    """
    Filters FASTQ reads by GC-content, read length, and average quality score.

    Arguments:
    seqs: dict — dictionary with reads in format {name: (sequence, quality_string)}
    gc_bounds: tuple or int — GC% range; if single number, treated as upper bound
    length_bounds: tuple or int — read length range; if single number, treated as lower bound
    quality_threshold: float — minimal average Q-score to keep the read

    Returns:
    dict — filtered reads {name: (sequence, quality_string)} that passed all conditions.

    Notes:
    - Sequences failing `is_nucleic_acid()` check are skipped.
    """
    if type(gc_bounds) in (int, float):
        gc_bounds = (0, gc_bounds)
    if type(length_bounds) in (int, float):
        length_bounds = (length_bounds, 2**32)

    gc_low, gc_high = gc_bounds
    length_low, length_high = length_bounds

    filtered = {}

    for name, (sequence, quality_str) in seqs.items():
        if not is_nucleic_acid(sequence):
            continue
        gc = gc_content(sequence)
        seq_len = len(sequence)
        avg_q = qscore(quality_str)

        gc_pass = (gc_low <= gc <= gc_high)
        length_pass = (length_low <= seq_len <= length_high)
        quality_pass = (avg_q >= quality_threshold)

        if gc_pass and length_pass and quality_pass:
            filtered[name] = (sequence, quality_str)

    return filtered