from modules.dna_rna_tools import is_nucleic_acid


def gc_content(sequence: str) -> float:
    """
    Calculates GC content for a nucleotide sequence.

    Arguments:
    sequence: str — input DNA or RNA sequence

    Returns:
    float — GC percentage (0–100). 
    Returns 0 for empty sequences.
    """
    sequence = str(sequence).upper()
    sequence_length = len(sequence)
    if sequence_length == 0:
        return 0
    
    g_count = sequence.count("G")
    c_count = sequence.count("C")

    gc_percent = ((g_count + c_count) / sequence_length) * 100
    return gc_percent


def qscore(quality_str: str) -> float:
    """
    Calculates the mean phred33 quality score for a given quality string.

    Arguments:
    quality_str: str — ASCII-encoded quality string

    Returns:
    float — mean Q-score across all bases.
    Returns 0.0 for empty strings.
    """
    str_length = len(str(quality_str))
    if str_length == 0:
        return 0.0
    q_scores = [ord(symbol) - 33 for symbol in quality_str]
    return sum(q_scores) / str_length