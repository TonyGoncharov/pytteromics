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


def calc_q_score(quality: str) -> float:
    quality = str(quality)
    quality_len = len(quality)
    if quality_len == 0:
        return 0.0
    q_scores = [ord(char) - 33 for char in quality]
    return sum(q_scores) / quality_len


def is_qscore_passing(
    quality: str, 
    quality_threshold: int | float = 0
) -> bool:
    q_score = calc_q_score(quality)
    return q_score >= quality_threshold