import pytest

from sequana.fastq_splitter import FastqSplitter

from . import test_dir

data = f"{test_dir}/data/fastq/test.fastq"
datagz = f"{test_dir}/data/fastq/test.fastq.gz"


def test_fastq_splitter_by_size(tmp_path):
    pattern = str(tmp_path / "chunk_#.fastq")
    splitter = FastqSplitter(datagz)
    splitter.by_size(reads_per_chunk=100, pattern=pattern)

    parts = list(tmp_path.glob("chunk_*.fastq"))
    assert len(parts) >= 1


def test_fastq_splitter_by_part(tmp_path):
    pattern = str(tmp_path / "part_#.fastq")
    splitter = FastqSplitter(datagz)
    splitter.by_part(n_parts=2, pattern=pattern)

    parts = sorted(tmp_path.glob("part_*.fastq"))
    assert len(parts) == 2


def test_fastq_splitter_gz_input(tmp_path):
    pattern = str(tmp_path / "chunk_#.fastq")
    splitter = FastqSplitter(datagz)
    splitter.by_size(reads_per_chunk=100, pattern=pattern)

    parts = list(tmp_path.glob("chunk_*.fastq"))
    assert len(parts) >= 1


def test_fastq_splitter_gzip_output(tmp_path):
    pattern = str(tmp_path / "chunk_#.fastq.gz")
    splitter = FastqSplitter(datagz)
    splitter.by_size(reads_per_chunk=100, pattern=pattern, gzip_output=True)

    parts = list(tmp_path.glob("chunk_*.fastq.gz"))
    assert len(parts) >= 1
