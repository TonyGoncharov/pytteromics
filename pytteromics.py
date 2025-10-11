from modules.filter_fastq import is_len_passing, is_gc_passing, is_qscore_passing
from modules.dna_rna_tools import is_nucleic_acid, validate_input, validate_mode, validate_function


def filter_fastq(
    seqs: dict[str, tuple[str, str]],
    gc_bounds: float | tuple[float, float] = (0, 100),
    length_bounds: int | tuple[int, int] = (0, 2**32),
    quality_threshold: int | float = 0
) -> dict[str, tuple[str, str]]: 
    filtered_seqs = {}
    for name, (sequence, quality) in seqs.items():
        if (
            is_len_passing(sequence, length_bounds)
            and is_gc_passing(sequence, gc_bounds)
            and is_qscore_passing(quality, quality_threshold)
        ):
            filtered_seqs[name] = (sequence, quality)
    return filtered_seqs


def run_dna_rna_tools(*args, behavior="mild"):
    sequences = validate_input(args)
    behavior = validate_mode(behavior)
    tool = validate_function(args[-1])

    result = []

    for seq in sequences:
        if is_nucleic_acid(seq):
            result.append(tool(seq))
        else:
            if behavior == "harsh":
                raise ValueError(f"Error: {seq} is not a DNA or RNA sequence!")
            result.append(False)

    if len(result) == 1:
        return result[0]
    return result