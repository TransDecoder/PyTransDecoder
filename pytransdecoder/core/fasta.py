"""
FASTA reading shared by all PyTransDecoder steps.

Text before the first '>' header line (e.g. ';' comment lines, as in FAIR Header Reference
files) is skipped with a warning, as Perl TransDecoder's Fasta_reader does. Biopython's "fasta"
parser rejects such files, and older Perl versions read the text as an extra sequence.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import IO, Iterator, Union

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)


def read_fasta(source: Union[str, Path, IO[str]]) -> Iterator[SeqRecord]:
    """Yield SeqRecords (id = first word of the header, description = whole header) from a
    path or an open text handle. Whitespace inside sequences is removed."""
    if isinstance(source, (str, Path)):
        with open(source) as handle:
            yield from _parse(handle, str(source))
    else:
        yield from _parse(source, getattr(source, "name", "FASTA input"))


def _parse(handle: IO[str], label: str) -> Iterator[SeqRecord]:
    header = None
    chunks = []
    skipped_text = False
    for line in handle:
        if line.startswith(">"):
            if header is not None:
                yield _record(header, chunks)
            elif skipped_text:
                logger.warning("%s: skipped text before the first '>' header line", label)
            header, chunks = line[1:].rstrip("\r\n"), []
        elif header is None:
            skipped_text = skipped_text or bool(line.strip())
        else:
            chunks.append("".join(line.split()))
    if header is None and skipped_text:
        logger.warning("%s: no '>' header line found; no sequences read", label)
    if header is not None:
        yield _record(header, chunks)


def _record(header: str, chunks) -> SeqRecord:
    words = header.split()
    seq_id = words[0] if words else ""
    return SeqRecord(Seq("".join(chunks)), id=seq_id, name=seq_id, description=header)
