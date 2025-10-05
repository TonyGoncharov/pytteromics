def is_nucleic_acid(sequence) -> bool:
    """
    Checks if the input sequence is DNA or RNA.

    Arguments:
    sequence: str — input sequence to check

    Returns:
    bool — True if sequence contains only 
    DNA (A, T, G, C) or RNA (A, U, G, C) bases.
    
    Returns False for empty, mixed, or invalid sequences.
    """
    if len(sequence) == 0:
        return False
    seq_set = set(sequence)
    dna = set("ATGCatgc")
    rna = set("AUGCaugc")
    return (seq_set <= dna) or (seq_set <= rna)