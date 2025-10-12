import os


def convert_multiline_fasta_to_oneline(
        input_fasta: str = 'data/example_multiline_fasta.fasta',
        output_fasta: str = 'data/oneline_fasta.fasta'
) -> str:
    os.makedirs(os.path.dirname(output_fasta), exist_ok = True)
    with open(input_fasta) as input_fasta, open(output_fasta, "a") as output_fasta:
        while True:
            header = input_fasta.readline()
            if header == "":
                break
            output_fasta.write(header)
            line = input_fasta.readline()
            while not line.startswith(">"):
                output_fasta.write(line.strip())

