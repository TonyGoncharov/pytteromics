# pytteromics

`pytteromics` is a small (but proud) Python package for basic bioinformatics tasks.

The package includes tools for:

🐍 Filtering FASTQ reads by length, GC content, and quality thresholds.

🐍 Processing DNA/RNA sequences — transcription, reverse complement, and validation.

🐍 Converting multiline FASTA files into one-line-per-sequence format.

🐍 Parsing BLAST output files to extract and organize protein descriptions.

🐍 Extracting neighboring genes from a GenBank (.gbk) file for one or more target genes of interest.

## Installation

Clone the repository and use the package locally:

```bash
git clone https://github.com/TonyGoncharov/pytteromics.git
cd pytteromics
```
## Examples

```Python
from pytteromics import filter_fastq

# Filter reads by GC content, length, and quality
filtered_path = filter_fastq(
    input_fastq = "data/example.fastq",
    output_fastq = "data/filtered.fastq",
    gc_bounds = (40, 60),
    length_bounds = (100, 1000),
    quality_threshold = 30,
    output_mode = "rewrite"
)
```

```Python
from pytteromics import run_dna_rna_tools

# Example: transcription (DNA → RNA)
result = run_dna_rna_tools("transcribe", "ATGCTTAA")
print(result)  # Output: AUGCUUAA

# Example: reverse complement for multiple sequences
results = run_dna_rna_tools("reverse_complement", "ATGC", "GGAATT")
print(results)  # Output: ['GCAT', 'AATTCC']
```


```Python
from pytteromics import convert_multiline_fasta_to_oneline

output_path = convert_multiline_fasta_to_oneline(
    input_fasta = "data/example_multiline.fasta",
    output_fasta = "data/oneline_fasta.fasta"
)
```


```Python
from pytteromics import parse_blast_output

parsed_path = parse_blast_output(
    input_file = "data/example_blast_results.txt",
    output_file = "data/parsed_blast_results.txt"
)
```


```Python
from pytteromics import select_genes_from_gbk_to_fasta

select_genes_from_gbk_to_fasta(
    input_gbk = "data/example_gbk.gbk",
    genes = ["ligB_1", "guaA"], # target gene names
    n_before = 1, # number of neighboring genes before each target
    n_after = 1, # number of neighboring genes after each target
    output_fasta = "data/selected_from_gbk.fasta"
)
```