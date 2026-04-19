# pytteromics

🎉 Meet the new implementation!

`pytteromics` is a small (but proud) Python package for basic bioinformatics tasks.

The package includes tools for:

🐍 FASTQ filtering by read length, GC content, and average quality threshold.

🐍 DNA sequence utilities: validation, reverse/complement, reverse-complement, and transcription to RNA.

🐍 RNA sequence utilities: validation, reverse/complement, and reverse-complement.

🐍 Protein sequence utilities: trypsin digestion into peptide fragments.

## Installation

Clone the repository and use the package locally:

```bash
git clone https://github.com/TonyGoncharov/pytteromics.git
cd pytteromics
```
## Examples

### FASTQ filtering
```python
from pytteromics import filter_fastq

# Filter reads by GC content, length, and quality
filter_fastq(
    input_fastq="data/example_fastq.fastq",
    output_fastq="data/filtered.fastq",
    gc_bounds=(40, 60),
    length_bounds=(50, 300),
    quality_threshold=30,
    output_mode="rewrite",
)
```

### DNA sequence operations
```python
from pytteromics import DNASequence

my_dna = DNASequence('ATGGCTGGTATTTGT')

print('complement: ', my_dna.complement())
print('reverse: ', my_dna.reverse())
print('reverse_complement: ', my_dna.reverse_complement())
print('transcribe: ', my_dna.transcribe())

# complement: TACCGACCATAAACA
# reverse: TGTTTATGGTCGGTA
# reverse_complement: ACAAATACCAGCCAT
# transcribe: AUGGCUGGUAUUUGU
```

### RNA sequence operations
```python
from pytteromics import RNASequence

my_rna = RNASequence('AUGGCUGGUAUUUGU')

print('complement: ', my_rna.complement())
print('reverse: ', my_rna.reverse())
print('reverse_complement: ', my_rna.reverse_complement())

# complement: UACCGACCAUAAACA
# reverse: UGUUUAUGGUCGGUA
# reverse_complement: ACAAAUACCAGCCAU
```

### Trypsin digestion
```python
from pytteromics import AminoAcidSequence

my_protein = AminoAcidSequence('MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPF')

print('tryptic peptides: ', my_protein.digest_with_trypsin())
# tryptic peptides: ['MK', 'WVTFISLLFLFSSAYSR', 'GVFR', 'R', 'DAHK', 'SEVAHR', 'FK', 'DLGEENFK', 'ALVLIAFAQYLQQCPF']
```
