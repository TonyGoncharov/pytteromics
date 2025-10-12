import os


def convert_multiline_fasta_to_oneline(
        input_fasta: str = 'data/example_multiline_fasta.fasta',
        output_fasta: str = 'data/oneline_fasta.fasta'
) -> str:
    """
    Converts a multi-line FASTA file into a one-line-per-sequence format
    and writes the result to a new file at the specified path.

    Args:
        input_fasta: Path to the input multi-line FASTA file.
        output_fasta: Path to the output one-line FASTA file.

    Returns:
        FASTA file with single-line sequences.
    """
    os.makedirs(os.path.dirname(output_fasta), exist_ok = True)
    with open(input_fasta) as input_fasta, open(output_fasta, "w") as output_fasta:
        header = None
        sequence = []
        for line in input_fasta:
            line = line.strip()
            if len(line) == 0:
                continue
            if line.startswith(">"):
                if header is not None:
                    output_fasta.write(header + "\n" + "".join(sequence) + "\n")
                header = line
                sequence = []
            else:
                sequence.append(line)
        if header is not None:
                output_fasta.write(header + "\n" + "".join(sequence) + "\n")
                output_fasta.write(line.strip())


def parse_blast_output(
        input_file: str = 'data/example_blast_results.txt', 
        output_file: str = 'data/parsed_blast_results.txt'
) -> str:
    """
    Parses a BLAST output file, extracts protein descriptions from the
    "Description" column, sorts them alphabetically, and writes the results
    to a new text file.

    Args:
        input_file: Path to the input BLAST results file.
        output_file: Path to the output file with parsed protein descriptions.

    Returns:
        Text file containing the sorted protein descriptions.
    """
    os.makedirs(os.path.dirname(output_file), exist_ok = True)
    with open(input_file) as input_file:
        proteins = []
        for line in input_file:
            if line.startswith("Description"):
                protein = input_file.readline().split("  ")[0]
                proteins.append(protein)
    proteins.sort()
    with open(output_file, "w") as output_file:
        for protein in proteins:
            output_file.write(protein + "\n")

