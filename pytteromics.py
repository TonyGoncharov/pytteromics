import os
from abc import ABC, abstractmethod


class BiologicalSequence(ABC):
    
    def __init__(self, sequence: str):
        self.sequence = sequence
        self._check_alphabet()

    def __len__(self) -> int:
        return len(self.sequence)   

    def __getitem__(self, key):
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

    def complement(self) -> str:
        compl_seq = []
        for i in range(0, len(self.sequence)):
            compl_seq.append(self._compl_rules[self.sequence[i]])
        return ''.join(compl_seq)
    
    def reverse(self) -> str:
        return self.sequence[::-1]

    def reverse_complement(self) -> str:
        return self.complement()[::-1]


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
        for i in range(0, len(self.sequence)):
            transc_seq.append(self._transc_rules[self.sequence[i]])
        return ''.join(transc_seq)
    

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
    