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