def gc_content(sequence) -> float:
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