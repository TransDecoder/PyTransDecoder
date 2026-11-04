"""External programs that some options run, and a check that they are on PATH.

pyTransdecoder runs these only when asked: gffread (--genome/--gtf), DIAMOND or NCBI BLAST+
(--blast-search-pep), HMMER (--pfam-search-db) and samtools (BAM structures). None of them can be
installed by pip. They are checked before any work starts, so a missing program stops the run at once
instead of after LongOrfs.
"""

import shutil
import sys
from typing import Iterable, List, Tuple

# program -> bioconda package that provides it
PACKAGE = {
    'gffread': 'gffread',
    'diamond': 'diamond',
    'makeblastdb': 'blast',
    'blastp': 'blast',
    'hmmsearch': 'hmmer',
    'hmmpress': 'hmmer',
    'samtools': 'samtools',
}


def missing_programs(needed: Iterable[Tuple[str, str]]) -> List[Tuple[str, str]]:
    """(program, reason) pairs whose program is not found on PATH, first occurrence of each program."""
    seen, missing = set(), []
    for program, reason in needed:
        if program in seen:
            continue
        seen.add(program)
        if shutil.which(program) is None:
            missing.append((program, reason))
    return missing


def require_programs(needed: Iterable[Tuple[str, str]]) -> None:
    """Exit with status 1 and one message listing every missing program, if any is missing."""
    missing = missing_programs(needed)
    if not missing:
        return
    lines = ["Error: required program(s) not found on PATH:"]
    for program, reason in missing:
        lines.append(f"  {program}  (needed for {reason})")
    packages = sorted({PACKAGE.get(p, p) for p, _ in missing})
    lines.append(f"Install with, e.g.: conda install -c bioconda {' '.join(packages)}")
    if any(p == 'diamond' for p, _ in missing):
        lines.append("(or use --blast-tool blastp with NCBI BLAST+)")
    print("\n".join(lines), file=sys.stderr)
    sys.exit(1)


def structures_need_samtools(path) -> bool:
    return str(path).endswith('.bam')
