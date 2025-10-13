from modules.filter_fastq import is_len_passing, is_gc_passing, is_qscore_passing, validate_input, validate_output
from modules.dna_rna_tools import is_nucleic_acid, prepare_sequences, prepare_tool, transcribe, reverse, complement, reverse_complement


def filter_fastq(
    input_fastq: str,
    output_fastq: str = "output_fastq.fastq",
    gc_bounds: float | tuple[float, float] = (0, 100),
    length_bounds: int | tuple[int, int] = (0, 2**32),
    quality_threshold: int | float = 0,
    output_mode: str = "append"
) -> None:
    """
    Filters reads in a FASTQ file by length, GC content, and quality score.

    Args:
        input_fastq: Path to the input FASTQ file.
        output_fastq: Path to the output FASTQ file.
        gc_bounds: Maximum or (min, max) GC content percentage.
        length_bounds: Minimum or (min, max) read length.
        quality_threshold: Minimum average quality score.
        output_mode: File writing mode ("append" or "rewrite").

    Returns:
        None
    """
    input_fastq = validate_input(input_fastq)
    output_fastq, output_mode = validate_output(output_fastq, output_mode)
    with open(input_fastq) as input_fastq, open(output_fastq, output_mode) as output_fastq:
        while True:
            header = input_fastq.readline()
            if header == "":
                break
            sequence = input_fastq.readline()
            comment = input_fastq.readline()
            quality = input_fastq.readline()

            if (
                is_len_passing(sequence, length_bounds)
                and is_gc_passing(sequence, gc_bounds)
                and is_qscore_passing(quality, quality_threshold)
            ):
                output_fastq.write(header + sequence + comment + quality)


def run_dna_rna_tools(*args: str) -> str | bool | list[str | bool]:
    """
    Applies a selected DNA/RNA processing tool to one or more input sequences.

    Args:
        *args: Input arguments, including one or more sequences 
        and selected tool name.

    Returns:
        A single result (str or bool) if one sequence is provided,
        otherwise a list of results for all sequences.
    """
    sequences = prepare_sequences(args)
    tool = prepare_tool(args)
    result = []
    for seq in sequences:
        if is_nucleic_acid(seq):
            result.append(tool(seq))
        else:
            result.append(False)
    if len(result) == 1:
        return result[0]
    return result
