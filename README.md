# pytteromics

`pytteromics` is a small (but proud) Python package for basic bioinformatics tasks on FASTQ data.
It provides simple utilities for analyzing nucleotide sequences, computing GC content, quality scores, and filtering reads based on sequence and quality thresholds.

## Installation

Clone the repository and use the package locally:

```bash
git clone https://github.com/TonyGoncharov/pytteromics.git
cd pytteromics
```

Then, in your Python script or Jupyter Notebook:

```Python
from modules.dna_rna_tools import is_nucleic_acid
from modules.filter_fastq import filter_fastq, gc_content, qscore
from tests.example_data import EXAMPLE_FASTQ
```

## Examples

```Python
filtered = filter_fastq(EXAMPLE_FASTQ, quality_threshold=30)
print(f"Reads before: {len(EXAMPLE_FASTQ)}, after: {len(filtered)}")

# Reads before: 12, after: 9
```

```Python
filtered = filter_fastq(EXAMPLE_FASTQ, gc_bounds=60)
for name in filtered.keys():
    print(name)

# @SRX079801
# @SRX079802
# @SRX079804
```


```Python
filtered = filter_fastq(
    EXAMPLE_FASTQ,
    gc_bounds=(40, 60),
    length_bounds=100,
    quality_threshold=30
)
```