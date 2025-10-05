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


def transcribe(sequence):
    compl_seq = []
    for i in range(0, len(sequence)):
        compl_seq.append(TRANSCR_RULES[sequence[i]])
    return "".join(compl_seq)


def reverse(sequence):
    return sequence[::-1]