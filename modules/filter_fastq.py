from modules.dna_rna_tools import is_nucleic_acid


def calc_gc_percent(sequence: str) -> float:
    sequence = str(sequence).upper()
    sequence_len = len(sequence)
    if sequence_len == 0:
        return 0.0
    g_count = sequence.count("G")
    c_count = sequence.count("C")
    return ((g_count + c_count) / sequence_len) * 100


def is_gc_passing(
    sequence: str, 
    gc_bounds: float | tuple[float, float] = (0, 100)
) -> bool:
    if type(gc_bounds) in (int, float):
        gc_bounds = (0, gc_bounds)
    min_gc, max_gc = gc_bounds
    gc_percent = calc_gc_percent(sequence)
    return min_gc <= gc_percent <= max_gc


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