# pytteromics

`pytteromics` is a small (but proud) Python package for basic bioinformatics tasks on FASTQ data.
It provides simple utilities for analyzing nucleotide sequences, computing GC content, quality scores, and filtering reads based on sequence and quality thresholds.

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