import os

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pytteromics import filter_fastq

def make_record(seq: str, quality: int, name: str = 'read') -> SeqRecord:
    """Create a minimal FASTQ record with uniform phred quality."""
    record = SeqRecord(Seq(seq), id=name, description='')
    record.letter_annotations['phred_quality'] = [quality] * len(seq)
    return record


def write_fastq(path: str, records: list[SeqRecord]) -> None:
    with open(path, 'w') as fh:
        SeqIO.write(records, fh, 'fastq')


def read_fastq(path: str) -> list[SeqRecord]:
    return list(SeqIO.parse(path, 'fastq'))


REC_PASS = make_record('GCGCATGCAT' * 5, quality=30, name='pass')
REC_LOW_GC = make_record('ATATATATAT' * 5, quality=30, name='low_gc')
REC_SHORT = make_record('GCGCA', quality=30, name='short')
REC_LOW_Q = make_record('GCGCATGCAT' * 5, quality=5, name='low_quality')


@pytest.fixture
def output_path(tmp_path) -> str:
    """Provide a path for the output FASTQ file (file is not created yet)."""
    return str(tmp_path / 'output.fastq')


@pytest.fixture
def single_record_fastq(tmp_path) -> str:
    """Write a single passing record to a temporary FASTQ file."""
    path = str(tmp_path / 'single.fastq')
    write_fastq(path, [REC_PASS])
    return path


@pytest.fixture
def multi_record_fastq(tmp_path) -> str:
    """Write all four test records to a temporary FASTQ file."""
    path = str(tmp_path / 'multi.fastq')
    write_fastq(path, [REC_PASS, REC_LOW_GC, REC_SHORT, REC_LOW_Q])
    return path


# Testing errors

class TestFilterFastqErrors:

    def test_missing_input_file_raises(self, output_path):
        """filter_fastq raises ValueError when the input file does not exist."""
        with pytest.raises(ValueError):
            filter_fastq(
                input_fastq='nonexistent_file.fastq',
                output_fastq=output_path
            )

    def test_invalid_output_mode_raises(self, single_record_fastq, output_path):
        """filter_fastq raises ValueError for an unsupported output_mode value."""
        with pytest.raises(ValueError):
            filter_fastq(
                input_fastq=single_record_fastq,
                output_fastq=output_path,
                output_mode='overwrite'
            )


# Testing input and output

class TestFilterFastqIO:

    def test_output_file_is_created_with_passed_records(
        self, single_record_fastq, output_path
    ):
        """Output file is created and contains exactly the one record that passed."""
        filter_fastq(single_record_fastq, output_path, output_mode='rewrite')

        assert os.path.isfile(output_path)
        records = read_fastq(output_path)
        assert len(records) == 1
        assert records[0].id == 'pass'

    def test_append_mode_accumulates_records(self, single_record_fastq, output_path):
        """Running filter_fastq twice in append mode doubles the record count."""
        filter_fastq(single_record_fastq, output_path, output_mode='append')
        filter_fastq(single_record_fastq, output_path, output_mode='append')

        assert len(read_fastq(output_path)) == 2

    def test_rewrite_mode_replaces_content(self, single_record_fastq, output_path):
        """Running filter_fastq twice in rewrite mode keeps the record count at 1."""
        filter_fastq(single_record_fastq, output_path, output_mode='rewrite')
        filter_fastq(single_record_fastq, output_path, output_mode='rewrite')

        assert len(read_fastq(output_path)) == 1


# Testing fitlering

class TestFilterFastqFiltering:

    @pytest.mark.parametrize('filter_kwargs,should_pass,should_fail', [
        ({'gc_bounds': (40, 60)}, 'pass', 'low_gc'),
        ({'length_bounds': (20, 80)}, 'pass', 'short'),
        ({'quality_threshold': 20}, 'pass', 'low_quality')
    ])
    def test_filter_passes_correct_and_rejects_failing_read(
        self, multi_record_fastq, output_path,
        filter_kwargs, should_pass, should_fail
    ):
        """Each filter criterion passes the expected read and rejects the failing one."""
        filter_fastq(
            input_fastq=multi_record_fastq,
            output_fastq=output_path,
            output_mode='rewrite',
            **filter_kwargs
        )

        passed_ids = []
        for read in read_fastq(output_path):
            passed_ids.append(read.id)
        assert should_pass in passed_ids
        assert should_fail not in passed_ids


if __name__ == '__main__':
    pytest.main(['-v'])