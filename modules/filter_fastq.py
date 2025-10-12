import os
from modules.dna_rna_tools import is_nucleic_acid


def calc_gc_percent(sequence: str) -> float:
    sequence = str(sequence).upper().strip()
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
    quality = str(quality).strip()
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


def calc_len_seq(sequence: str) -> int:
    sequence = str(sequence).strip()
    return len(sequence)


def is_len_passing(
    sequence: str,
    length_bounds: int | tuple[int, int] = (0, 2**32)
) -> bool:
    if type(length_bounds) in (int, float):
        length_bounds = (length_bounds, 2**32)
    min_len, max_len = length_bounds
    return min_len <= calc_len_seq(sequence) <= max_len


def validate_input(file_path):
    if not os.path.isfile(file_path):
        raise ValueError("Error: at least 2 arguments expected!")
    return file_path


def validate_output(file_path, output_mode):
    if not os.path.exists("filtered"):
        os.mkdir("filtered")
    output_fastq = os.path.join("filtered", file_path)

    if output_mode == "append":
        output_mode = "a"
    elif output_mode == "rewrite":
        output_mode = "w"
    return output_fastq, output_mode