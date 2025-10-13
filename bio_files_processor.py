import os


def convert_multiline_fasta_to_oneline(
        input_fasta: str = 'data/example_multiline_fasta.fasta',
        output_fasta: str = 'data/oneline_fasta.fasta'
) -> None:
    """
    Converts a multi-line FASTA file into a one-line-per-sequence format
    and writes the result to a new file at the specified path.

    Args:
        input_fasta: Path to the input multi-line FASTA file.
        output_fasta: Path to the output one-line FASTA file.

    Returns:
        None
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


def parse_blast_output(
        input_file: str = 'data/example_blast_results.txt', 
        output_file: str = 'data/parsed_blast_results.txt'
) -> None:
    """
    Parses a BLAST output file, extracts protein descriptions from the
    "Description" column, sorts them alphabetically, and writes the results
    to a new text file.

    Args:
        input_file: Path to the input BLAST results file.
        output_file: Path to the output file with parsed protein descriptions.

    Returns:
        None
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


def select_genes_from_gbk_to_fasta(
        input_gbk: str = "data/example_gbk.gbk",
        genes: str | list = ["ligB_1", "guaA"],
        n_before: int = 1,
        n_after: int = 1,
        output_fasta: str = 'data/selected_from_gbk.fasta'
) -> None:
    """
    Extracts neighboring gene sequences from a GenBank (.gbk) file and
    writes them to a FASTA file.

    Args:
        input_gbk: Path to the input GenBank file.
        genes: Target gene name or a list of gene names of interest.
        n_before: Number of neighboring genes to include before each target.
        n_after: Number of neighboring genes to include after each target.
        output_fasta: Path to the output FASTA file.

    Returns:
        None.
    """
    if type(genes) == "str":
        genes = [genes]

    with open(input_gbk) as input_gbk:
        locus = None
        gene = None
        features = {}
        indices = {}
        for line in input_gbk:
            if line.startswith("LOCUS"):
                locus = line.split()[1]
                features[locus] = [[], []]
                gene = None
            elif "/gene=" in line:
                gene = line.split("=")[1].replace('"', '').strip()
                features[locus][0].append(gene)
                if gene in genes:
                    index = features[locus][0].index(gene)
                    indices[locus] = [[], [], []]
                    indices[locus][0].append(index)
                    indices[locus][1].append(index - n_before)
                    indices[locus][2].append(index + n_after)
            elif "/translation" in line and gene is not None:
                translation = line.split("=")[1].replace('"', '').strip()
                features[locus][1].append(translation)
                gene = None
            elif line.startswith("CDS"):
                gene = None
    with open(output_fasta, "w") as output_fasta:
        for locus in indices:
            gene_indices = indices[locus][0]
            left_indices = indices[locus][1]
            right_indices = indices[locus][2]
            for i in range(len(gene_indices)):
                left_index = max(0, left_indices[i])
                right_index = min(len(features[locus][0]) - 1, right_indices[i])
                for j in range(left_index, right_index + 1):
                    gene_name = features[locus][0][j]
                    sequence = features[locus][1][j]
                    if gene_name in genes:
                        continue
                    output_fasta.write(f">{locus}|{gene_name}\n")
                    output_fasta.write(f"{sequence}\n")
                    