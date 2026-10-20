"""
Coding gene models (exons + a CDS range) and their GFF3 / BED representation.

Used for TransDecoder's transcript-level outputs (one exon spanning the transcript) and for
ORFs propagated to the genome through transcript structures. The GFF3 layout follows Perl
TransDecoder's Gene_obj::to_GFF3_format: gene, then per isoform mRNA, five_prime_UTR, exon (+
CDS) in transcription order, three_prime_UTR; a blank line after each gene.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

Segment = Tuple[int, int]


@dataclass
class CodingModel:
    """One transcript isoform with a CDS. Coordinates are 1-based inclusive on `seqid`;
    exons are (lend, rend) with lend <= rend; cds_lend/cds_rend bound the CDS (incl. stop)."""
    seqid: str
    strand: str
    exons: List[Segment]
    cds_lend: int
    cds_rend: int
    model_id: str
    name: str = ""
    extra_attrs: str = ""          # appended to the mRNA attributes (e.g. "score=...;...")

    def __post_init__(self):
        if self.strand not in ("+", "-"):
            raise ValueError(f"strand must be + or -, not {self.strand!r}")
        self.exons = sorted((min(a, b), max(a, b)) for a, b in self.exons)
        self.cds_lend, self.cds_rend = sorted((self.cds_lend, self.cds_rend))

    # -- ordering helpers -------------------------------------------------------------
    def _in_tx_order(self, segments: Sequence[Segment]) -> List[Segment]:
        segs = sorted(segments)
        return segs if self.strand == "+" else segs[::-1]

    def exons_tx_order(self) -> List[Segment]:
        return self._in_tx_order(self.exons)

    def cds_segments(self) -> List[Segment]:
        """CDS pieces per exon, in transcription order."""
        segs = [(max(a, self.cds_lend), min(b, self.cds_rend)) for a, b in self.exons
                if b >= self.cds_lend and a <= self.cds_rend]
        if not segs:
            raise ValueError(f"{self.model_id}: CDS {self.cds_lend}-{self.cds_rend} lies outside its exons")
        return self._in_tx_order(segs)

    def _utr(self, left: bool) -> List[Segment]:
        if left:
            segs = [(a, min(b, self.cds_lend - 1)) for a, b in self.exons if a < self.cds_lend]
        else:
            segs = [(max(a, self.cds_rend + 1), b) for a, b in self.exons if b > self.cds_rend]
        return self._in_tx_order(segs)

    def utr5(self) -> List[Segment]:
        return _min_utr(self._utr(left=(self.strand == "+")))

    def utr3(self) -> List[Segment]:
        return _min_utr(self._utr(left=(self.strand == "-")))

    def span(self) -> Segment:
        return self.exons[0][0], self.exons[-1][1]

    # -- output -----------------------------------------------------------------------
    def gff3_lines(self, gene_id: str, source: str) -> List[str]:
        lend, rend = self.span()
        name = _quoted(self.name)
        attrs = f"ID={self.model_id};Parent={gene_id}" + (f";Name={name}" if name else "") + \
            (f";{self.extra_attrs}" if self.extra_attrs else "")
        lines = [_row(self.seqid, source, "mRNA", lend, rend, self.strand, ".", attrs)]
        for i, (a, b) in enumerate(self.utr5(), 1):
            lines.append(_row(self.seqid, source, "five_prime_UTR", a, b, self.strand, ".",
                              f"ID={self.model_id}.utr5p{i};Parent={self.model_id}"))
        cds = dict()
        cds_len = 0
        for a, b in self.cds_segments():
            # GFF3 phase: bases to skip from this segment's 5' end to reach a codon start
            cds[(a, b)] = (3 - cds_len % 3) % 3
            cds_len += b - a + 1
        for i, (ea, eb) in enumerate(self.exons_tx_order(), 1):
            lines.append(_row(self.seqid, source, "exon", ea, eb, self.strand, ".",
                              f"ID={self.model_id}.exon{i};Parent={self.model_id}"))
            for (a, b), phase in cds.items():
                if ea <= a and b <= eb:
                    lines.append(_row(self.seqid, source, "CDS", a, b, self.strand, str(phase),
                                      f"ID=cds.{self.model_id};Parent={self.model_id}"))
        for i, (a, b) in enumerate(self.utr3(), 1):
            lines.append(_row(self.seqid, source, "three_prime_UTR", a, b, self.strand, ".",
                              f"ID={self.model_id}.utr3p{i};Parent={self.model_id}"))
        return lines

    def bed12(self, name: str) -> str:
        """BED12 line (as Perl gff3_file_to_bed.pl: thickStart/End = CDS span)."""
        lend, rend = self.span()
        sizes = ",".join(str(b - a + 1) for a, b in self.exons)
        starts = ",".join(str(a - lend) for a, _ in self.exons)
        return "\t".join(map(str, [self.seqid, lend - 1, rend, name, 0, self.strand,
                                   self.cds_lend - 1, self.cds_rend, 0, len(self.exons), sizes, starts]))


@dataclass
class CodingGene:
    gene_id: str
    models: List[CodingModel] = field(default_factory=list)
    name: str = ""

    def gff3_text(self, source: str = "transdecoder") -> str:
        m0 = self.models[0]
        lend = min(m.span()[0] for m in self.models)
        rend = max(m.span()[1] for m in self.models)
        name = _quoted(self.name if self.name else m0.name)
        gene_attrs = f"ID={_quote_id(self.gene_id)}" + (f";Name={name}" if name else "")
        lines = [_row(m0.seqid, source, "gene", lend, rend, m0.strand, ".", gene_attrs)]
        for m in self.models:
            lines += m.gff3_lines(_quote_id(self.gene_id), source)
        return "\n".join(lines) + "\n\n"


def _min_utr(segments: List[Segment]) -> List[Segment]:
    # as Perl Gene_obj::has_5prime_UTR / has_3prime_UTR: a UTR needs more than 2 nt in total,
    # so the 1-2 nt left by a partial ORF's frame offset is not reported
    return segments if sum(b - a + 1 for a, b in segments) > 2 else []


def _row(seqid, source, ftype, a, b, strand, phase, attrs) -> str:
    return f"{seqid}\t{source}\t{ftype}\t{a}\t{b}\t.\t{strand}\t{phase}\t{attrs}"


def _quoted(name: Optional[str]) -> str:
    if not name:
        return ""
    if name[0] in "'\"" and name[-1] in "'\"":
        return name
    return f'"{name}"'


def _quote_id(identifier: str) -> str:
    return f'"{identifier}"' if ";" in identifier else identifier
