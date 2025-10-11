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


def is_nucleic_acid(sequence: str) -> bool:
    if len(sequence) == 0:
        return False
    seq_set = set(sequence)
    dna = set("ATGCatgc")
    rna = set("AUGCaugc")
    return (seq_set <= dna) or (seq_set <= rna)


def transcribe(sequence: str) -> str:
    compl_seq = []
    for i in range(0, len(sequence)):
        compl_seq.append(TRANSCR_RULES[sequence[i]])
    return "".join(compl_seq)


def reverse(sequence: str) -> str:
    return sequence[::-1]


def complement(sequence: str) -> str:
    compl_seq = []
    if "U" in sequence.upper():
        for i in range(0, len(sequence)):
            compl_seq.append(RNA_COMPL_RULES[sequence[i]])
        return "".join(compl_seq)
    else:
        for i in range(0, len(sequence)):
            compl_seq.append(DNA_COMPL_RULES[sequence[i]])
        return "".join(compl_seq)
    

def reverse_complement(sequence: str) -> str:
    compl_seq = complement(sequence)
    return reverse(compl_seq)
    

def prepare_sequences(args: tuple[str, ...]) -> tuple[str, ...]:
    if len(args) < 2:
        raise ValueError("Error: at least 2 arguments expected!")
    return args[:-1]


def prepare_tool(args: tuple[str, ...]) -> callable:
    tool = args[-1]
    TOOLS = {
    "is_nucleic_acid": is_nucleic_acid,
    "transcribe": transcribe,
    "reverse": reverse,
    "complement": complement,
    "reverse_complement": reverse_complement,
}
    if tool not in TOOLS:
        raise ValueError("Error: " f"{tool} is unsupported tool!")
    return TOOLS[tool]