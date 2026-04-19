import os
from abc import ABC, abstractmethod
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction


class BiologicalSequence(ABC):
    
    def __init__(self, sequence: str):
        self.sequence = sequence
        self._check_alphabet()

    def __len__(self) -> int:
        return len(self.sequence)   

    def __getitem__(self, key: int | slice) -> str:
        return self.sequence[key]
    
    def __str__(self) -> str:
        return self.sequence
    
    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self.sequence!r})'
    
    @property
    @abstractmethod
    def alphabet(self) -> set[str]:
        '''specific alphabet is defined in the descendant classes'''
        raise NotImplementedError
    
    def _check_alphabet(self) -> None:
        invalid = set()
        for char in self.sequence:
            if char not in self.alphabet:
                invalid.add(char)
        if invalid:
            raise ValueError(f'Invalid symbols for {self.__class__.__name__}: {sorted(invalid)}')
        

class NucleicAcidSequence(BiologicalSequence):

    @property
    @abstractmethod
    def _compl_rules(self) -> dict[str, str]:
        '''specific rules are defined in the descendant classes'''
        raise NotImplementedError

    def complement(self) -> 'NucleicAcidSequence':
        compl_seq = []
        for char in self.sequence:
            compl_seq.append(self._compl_rules[char])
        return self.__class__(''.join(compl_seq))
    
    def reverse(self) -> 'NucleicAcidSequence':
        return self.__class__(self.sequence[::-1])

    def reverse_complement(self) -> 'NucleicAcidSequence':
        return self.__class__(str(self.complement())[::-1])


class DNASequence(NucleicAcidSequence):

    @property
    def alphabet(self) -> set[str]:
        return set('ATCGatcg')

    @property
    def _compl_rules(self) -> dict[str, str]:
        return {
            'A': 'T',
            'T': 'A',
            'G': 'C',
            'C': 'G',
            'a': 't',
            't': 'a',
            'g': 'c',
            'c': 'g',
        }
    
    @property
    def _transc_rules(self) -> dict[str, str]:
        return {
            'A': 'A',
            'T': 'U',
            'G': 'G',
            'C': 'C',
            'a': 'a',
            't': 'u',
            'g': 'g',
            'c': 'c',
        }
    
    def transcribe(self) -> str:
        transc_seq = []
        for char in self.sequence:
            transc_seq.append(self._transc_rules[char])
        return RNASequence(''.join(transc_seq))
    

class RNASequence(NucleicAcidSequence):

    @property
    def alphabet(self) -> set[str]:
        return set('AUCGaucg')

    @property
    def _compl_rules(self) -> dict[str, str]:
        return {
            'A': 'U',
            'U': 'A',
            'G': 'C',
            'C': 'G',
            'a': 'u',
            'u': 'a',
            'g': 'c',
            'c': 'g',
        }
    

class AminoAcidSequence(BiologicalSequence):
    
    @property
    def alphabet(self) -> set[str]:
        return set('ACDEFGHIKLMNPQRSTVWYOU')
    
    def digest_with_trypsin(self) -> list[str]:
        peptides = []
        pep_start = 0

        for i in range(len(self.sequence)):
            aa = self.sequence[i]
            if aa in ('K', 'R'):
                peptides.append(self.sequence[pep_start:i+1])
                pep_start = i + 1

        if pep_start < len(self.sequence):
            peptides.append(self.sequence[pep_start:])

        return peptides


def _mean_quality(record: SeqRecord) -> float:
    qualities = record.letter_annotations.get('phred_quality')
    if qualities is None:
        raise ValueError('FASTQ record does not contain phred_quality annotations')
    return sum(qualities) / len(qualities)


def _validate_input(file_path: str) -> str:
    if not os.path.isfile(file_path):
        raise ValueError(f"File '{file_path}' does not exist")
    return file_path


def _prepare_output(file_path: str, output_mode: str) -> tuple[str, str]:
    output_dir = os.path.dirname(file_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    output_path = file_path

    if output_mode == 'append':
        mode = 'a'
    elif output_mode == 'rewrite':
        mode = 'w'
    else:
        raise ValueError("output_mode must be 'append' or 'rewrite'")
    return output_path, mode


def filter_fastq(
    input_fastq: str,
    output_fastq: str = 'output_fastq.fastq',
    gc_bounds: float | tuple[float, float] = (0, 100),
    length_bounds: int | tuple[int, int] = (0, 2**32),
    quality_threshold: int | float = 0,
    output_mode: str = 'append'
) -> None:
    """
    Filters reads in a FASTQ file by length, GC content, and quality score.

    Args:
        input_fastq: Path to the input FASTQ file.
        output_fastq: Path to the output FASTQ file.
        gc_bounds: Maximum or (min, max) GC content percentage.
        length_bounds: Minimum or (min, max) read length.
        quality_threshold: Minimum average quality score.
        output_mode: File writing mode ("append" or "rewrite").

    Returns:
        None
    """
    input_fastq = _validate_input(input_fastq)
    output_fastq, mode = _prepare_output(output_fastq, output_mode)

    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)
    gc_min, gc_max = gc_bounds

    if isinstance(length_bounds, (int, float)):
        length_bounds = (0, length_bounds)
    len_min, len_max = length_bounds

    passed_records: list[SeqRecord] = []

    for record in SeqIO.parse(input_fastq, 'fastq'):
        seq_len = len(record.seq)

        if not (len_min <= seq_len <= len_max):
            continue

        gc_percent = gc_fraction(record.seq) * 100
        if not (gc_min <= gc_percent <= gc_max):
            continue

        if _mean_quality(record) < quality_threshold:
            continue

        passed_records.append(record)

    with open(output_fastq, mode) as out_handle:
        SeqIO.write(passed_records, out_handle, 'fastq')
