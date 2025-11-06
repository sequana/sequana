#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2022 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
# sequana/fastq_tools/splitter.py
import gzip
import os
from pathlib import Path
from typing import Iterator, Union

import colorlog
from tqdm import tqdm

logger = colorlog.getLogger(__name__)


class FastqSplitter:
    def __init__(self, filename: Union[str, Path]):
        self.filename = Path(filename)
        self.is_gz = self.filename.suffix == ".gz"

    def _open(self):
        if self.is_gz:
            return gzip.open(self.filename, "rt")
        return open(self.filename, "r")

    def _read_fastq(self) -> Iterator[list[str]]:
        """Yield one FASTQ record (4 lines)."""
        with self._open() as fh:
            while True:
                record = [fh.readline() for _ in range(4)]
                if not record[0]:
                    break
                yield record

    def _open_outfile(self, pattern: str, part_num: int, gzip_output: bool):
        out_name = pattern.replace("#", str(part_num))
        out_path = Path(out_name)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        mode = "wt"
        if gzip_output and not out_path.suffix == ".gz":
            out_path = out_path.with_suffix(out_path.suffix + ".gz")
        if gzip_output:
            return gzip.open(out_path, mode, compresslevel=5)
        return open(out_path, mode)

    def get_nreads(self):
        logger.info("Computing read counts")
        from sequana.tools import GZLineCounter

        count = GZLineCounter(self.filename)
        count = count._use_gzip() / 4
        print(count)
        return count

    def by_size(self, reads_per_chunk: int, pattern: str, gzip_output: bool = False, buffer_size: int = 10000):
        """Split FASTQ into chunks of N reads with file-based progress."""
        # rough estimate of total reads
        nreads = self.get_nreads()
        est_total_files = max(1, nreads // reads_per_chunk)

        chunk_id = 1
        out_fh = self._open_outfile(pattern, chunk_id, gzip_output)
        buffer = []
        reads_written = 0

        progress = tqdm(total=est_total_files, desc="Splitting FASTQ", unit="file")

        for record in self._read_fastq():
            buffer.extend(record)
            reads_written += 1

            if len(buffer) >= 4 * buffer_size:
                out_fh.writelines(buffer)
                buffer.clear()

            if reads_written >= reads_per_chunk:
                if buffer:
                    out_fh.writelines(buffer)
                    buffer.clear()
                out_fh.close()
                chunk_id += 1
                out_fh = self._open_outfile(pattern, chunk_id, gzip_output)
                reads_written = 0
                progress.update(1)

        if buffer:
            out_fh.writelines(buffer)
        out_fh.close()
        progress.update(1)  # last file
        progress.close()

    def by_part(self, n_parts: int, pattern: str, gzip_output: bool = False, buffer_size: int = 10000):
        """Split FASTQ into N parts with approximate read counting."""
        nreads = self.get_nreads()
        reads_per_part = max(1, nreads // n_parts)

        chunk_id = 1
        out_fh = self._open_outfile(pattern, chunk_id, gzip_output)
        buffer = []
        reads_written = 0

        progress = tqdm(total=n_parts, desc="Splitting FASTQ", unit="file")

        for record in self._read_fastq():
            buffer.extend(record)
            reads_written += 1

            # write buffer for efficiency
            if len(buffer) >= 4 * buffer_size:
                out_fh.writelines(buffer)
                buffer.clear()

            # switch file when reads_per_part is reached, except for last chunk
            if reads_written >= reads_per_part and chunk_id < n_parts:
                if buffer:
                    out_fh.writelines(buffer)
                    buffer.clear()
                out_fh.close()
                chunk_id += 1
                out_fh = self._open_outfile(pattern, chunk_id, gzip_output)
                reads_written = 0
                progress.update(1)

        # write remaining reads to the last chunk
        if buffer:
            out_fh.writelines(buffer)
        out_fh.close()
        progress.update(1)
        progress.close()
