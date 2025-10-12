import os


def convert_multiline_fasta_to_oneline(
        input_fasta: str = 'data/example_multiline_fasta.fasta',
        output_fasta: str = 'data/oneline_fasta.fasta'
) -> str:
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

