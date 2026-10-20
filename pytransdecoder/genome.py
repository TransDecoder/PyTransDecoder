"""
Propagate ORFs predicted on transcript sequences to genome coordinates through the
transcripts' genome structures, producing gene / mRNA / exon / CDS / UTR models.

Transcript structures can come from
  - GTF (exon features with gene_id / transcript_id: StringTie, LRAA, Cufflinks, GENCODE ...),
  - alignment GFF3 (cDNA_match / match lines with Target=<transcript> <start> <end>: GMAP,
    PASA, util/gtf_to_alignment_gff3.pl),
  - SAM / BAM spliced alignments (minimap2 -ax splice; primary alignments; BAM read via samtools).

For GTF and alignment GFF3 the result is the same as Perl TransDecoder's
util/cdna_alignment_orf_to_genome_orf.pl, except that each mRNA is named by its own ORF.
Transcripts of unknown strand ('.', e.g. StringTie single-exon transcripts) are taken to be
on the strand their sequence was extracted from: '+' for gffread (default), '-' for
TransDecoder's util/gtf_genome_to_cdna_fasta.pl (unstranded_as='-'). Alignment blocks carry transcript coordinates, so
transcripts with indels or clipped ends are mapped base-accurately; an ORF whose CDS is not
entirely within the aligned part of its transcript is not propagated.
"""

from __future__ import annotations

import logging
import re
import subprocess
from collections import OrderedDict, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from pytransdecoder.core.gene_model import CodingGene, CodingModel

logger = logging.getLogger(__name__)

# (genome lend, genome rend, transcript start, transcript end); transcript coordinates run in
# the transcribed direction: on '-' structures they decrease with genome position
Block = Tuple[int, int, int, int]

# alignment blocks separated by at most this many genome bases (an alignment deletion) belong
# to the same exon; larger gaps are introns
MAX_INDEL = 20


@dataclass
class TranscriptStructure:
    transcript_id: str
    seqid: str
    strand: str                      # orientation of the transcript sequence on the genome
    blocks: List[Block] = field(default_factory=list)
    gene_id: str = ""
    exons: List[Tuple[int, int]] = field(default_factory=list)   # merged across indels

    def finalize(self, transcript_length: Optional[int] = None) -> None:
        """Order blocks; assign transcript coordinates when the source has none (GTF)."""
        self.blocks.sort()
        if self.blocks and self.blocks[0][2] == 0:
            # contiguous transcript: cumulative exon lengths in the transcribed direction
            ordered = self.blocks if self.strand != "-" else self.blocks[::-1]
            pos, blocks = 1, []
            for g0, g1, _, _ in ordered:
                n = g1 - g0 + 1
                blocks.append((g0, g1, pos, pos + n - 1) if self.strand != "-" else (g0, g1, pos + n - 1, pos))
                pos += n
            self.blocks = sorted(blocks)
        if not self.exons:
            self.exons = _merge(((g0, g1) for g0, g1, _, _ in self.blocks), max_gap=MAX_INDEL)
        self.exons = _merge(self.exons)

    def to_genome(self, tpos: int) -> Optional[int]:
        """Genome position of transcript position tpos (1-based), None if not aligned."""
        for g0, g1, t0, t1 in self.blocks:
            lo, hi = min(t0, t1), max(t0, t1)
            if lo <= tpos <= hi:
                return g0 + (tpos - t0) if t1 >= t0 else g0 + (t0 - tpos)
        return None


def _merge(segments: Iterable[Tuple[int, int]], max_gap: int = 0) -> List[Tuple[int, int]]:
    merged: List[List[int]] = []
    for a, b in sorted(segments):
        if merged and a <= merged[-1][1] + 1 + max_gap:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    return [(a, b) for a, b in merged]


# --------------------------------------------------------------------------------------
# readers
# --------------------------------------------------------------------------------------

def _gtf_attrs(text: str) -> Dict[str, str]:
    attrs = {}
    for part in text.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        key, _, val = part.partition(" ")
        attrs.setdefault(key, val.strip().strip('"'))
    return attrs


def read_gtf(path: Path, unstranded_as: str = "+") -> Dict[str, TranscriptStructure]:
    structs: Dict[str, TranscriptStructure] = OrderedDict()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            x = line.rstrip("\n").split("\t")
            if len(x) < 9 or x[2] != "exon":
                continue
            a = _gtf_attrs(x[8])
            tid = a.get("transcript_id")
            if not tid:
                raise ValueError(f"GTF exon without transcript_id: {line.strip()}")
            s = structs.get(tid)
            if s is None:
                strand = x[6] if x[6] in ("+", "-") else unstranded_as
                s = structs[tid] = TranscriptStructure(tid, x[0], strand, gene_id=a.get("gene_id", ""))
            s.blocks.append((int(x[3]), int(x[4]), 0, 0))
    return structs


def read_alignment_gff3(path: Path) -> Dict[str, TranscriptStructure]:
    structs: Dict[str, TranscriptStructure] = OrderedDict()
    target_re = re.compile(r"Target=([^;\s]+)(?:\s+(\d+)\s+(\d+))?")
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            x = line.rstrip("\n").split("\t")
            if len(x) < 9:
                continue
            m = target_re.search(x[8])
            if not m:
                raise ValueError(f"alignment GFF3 line without Target=: {line.strip()}")
            tid, gene_id = m.group(1), ""
            gt = re.match(r"GENE\^(\S+),TRANS\^(\S+)", tid)   # util/gtf_to_alignment_gff3.pl
            if gt:
                gene_id, tid = gt.group(1), gt.group(2)
            s = structs.get(tid)
            if s is None:
                s = structs[tid] = TranscriptStructure(tid, x[0], x[6], gene_id=gene_id)
            g0, g1 = int(x[3]), int(x[4])
            if m.group(2):
                t0, t1 = int(m.group(2)), int(m.group(3))
                t0, t1 = (min(t0, t1), max(t0, t1)) if x[6] != "-" else (max(t0, t1), min(t0, t1))
                s.blocks.append((g0, g1, t0, t1))
            else:
                s.blocks.append((g0, g1, 0, 0))
    return structs


_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def read_sam(path: Path) -> Dict[str, TranscriptStructure]:
    """Primary spliced alignments; BAM input is read through `samtools view`."""
    structs: Dict[str, TranscriptStructure] = OrderedDict()
    if str(path).endswith(".bam"):
        proc = subprocess.Popen(["samtools", "view", "-F", "0x904", str(path)],
                                stdout=subprocess.PIPE, text=True)
        lines = proc.stdout
    else:
        proc, lines = None, open(path)
    try:
        for line in lines:
            if line.startswith("@"):
                continue
            x = line.split("\t", 11)
            flag = int(x[1])
            if flag & 0x904 or x[2] == "*":
                continue
            tid, seqid, gpos, cigar = x[0], x[2], int(x[3]), x[5]
            ops = [(int(n), op) for n, op in _CIGAR_RE.findall(cigar)]
            qlen = sum(n for n, op in ops if op in "MIS=XH")
            minus = bool(flag & 16)
            strand = "-" if minus else "+"
            q, g = 0, gpos
            blocks, exons, exon_start = [], [], gpos
            for n, op in ops:
                if op in "M=X":
                    t0, t1 = (q + 1, q + n) if not minus else (qlen - q, qlen - q - n + 1)
                    blocks.append((g, g + n - 1, t0, t1))
                    q += n
                    g += n
                elif op in "ISH":
                    q += n
                elif op == "D":
                    g += n
                elif op == "N":
                    exons.append((exon_start, g - 1))
                    g += n
                    exon_start = g
            exons.append((exon_start, g - 1))
            s = TranscriptStructure(tid, seqid, strand, blocks=blocks, exons=exons)
            structs[tid] = s
    finally:
        if proc is not None:
            proc.stdout.close()
            if proc.wait() != 0:
                raise RuntimeError(f"samtools view {path} failed")
        else:
            lines.close()
    return structs


def read_structures(path: Path, fmt: Optional[str] = None,
                    unstranded_as: str = "+") -> Dict[str, TranscriptStructure]:
    path = Path(path)
    if fmt is None:
        name = path.name.lower()
        if name.endswith((".bam", ".sam")):
            fmt = "sam"
        elif name.endswith(".gtf"):
            fmt = "gtf"
        else:
            fmt = "gff3"
    structs = read_gtf(path, unstranded_as) if fmt == "gtf" else \
        {"gff3": read_alignment_gff3, "sam": read_sam}[fmt](path)
    for s in structs.values():
        s.finalize()
    return structs


@dataclass
class TranscriptOrf:
    transcript_id: str
    model_id: str
    gene_id: str
    name: str
    strand: str
    cds_lend: int
    cds_rend: int


def read_orfs(td_gff3: Path) -> "OrderedDict[str, List[TranscriptOrf]]":
    """ORFs of a TransDecoder transcript-level GFF3 (Perl or PyTransDecoder), by transcript."""
    gene_names: Dict[str, str] = {}
    mrnas: "OrderedDict[str, dict]" = OrderedDict()
    cds: Dict[str, List[int]] = defaultdict(list)
    with open(td_gff3) as fh:
        for line in fh:
            x = line.rstrip("\n").split("\t")
            if line.startswith("#") or len(x) != 9:
                continue
            attrs = dict(kv.split("=", 1) for kv in x[8].split(";") if "=" in kv)
            if x[2] == "gene":
                gene_names[attrs["ID"]] = attrs.get("Name", "").strip('"')
            elif x[2] == "mRNA":
                mrnas[attrs["ID"]] = {"seqid": x[0], "gene": attrs.get("Parent", ""), "strand": x[6]}
            elif x[2] == "CDS":
                cds[attrs["Parent"]] += [int(x[3]), int(x[4])]
    orfs: "OrderedDict[str, List[TranscriptOrf]]" = OrderedDict()
    for model_id, m in mrnas.items():
        if not cds.get(model_id):
            continue
        orf = TranscriptOrf(m["seqid"], model_id, m["gene"], gene_names.get(m["gene"], ""),
                            m["strand"], min(cds[model_id]), max(cds[model_id]))
        orfs.setdefault(m["seqid"], []).append(orf)
    return orfs


def read_fasta_lengths(path: Path) -> Dict[str, int]:
    lengths, name, n = {}, None, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = n
                name, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
    if name is not None:
        lengths[name] = n
    return lengths


# --------------------------------------------------------------------------------------
# propagation
# --------------------------------------------------------------------------------------

def place_orf(orf: TranscriptOrf, struct: TranscriptStructure, transcript_length: int) -> Optional[CodingModel]:
    """Genome model of one ORF, or None (logged) if it cannot be placed."""
    cds_lend, cds_rend = orf.cds_lend, min(orf.cds_rend, transcript_length)
    tstrand = struct.strand if struct.strand in "+-" else "+"
    if orf.strand == "+":
        strand = tstrand
    else:
        if len(struct.exons) > 1:
            logger.warning("%s: ORF on the opposite strand of spliced transcript %s; skipped",
                           orf.model_id, struct.transcript_id)
            return None
        strand = "-" if tstrand == "+" else "+"   # single exon: orient by the ORF
    g_a, g_b = struct.to_genome(cds_lend), struct.to_genome(cds_rend)
    if g_a is None or g_b is None:
        logger.warning("%s: CDS %d-%d of %s is not entirely in its aligned part; skipped",
                       orf.model_id, cds_lend, cds_rend, struct.transcript_id)
        return None
    return CodingModel(seqid=struct.seqid, strand=strand, exons=list(struct.exons),
                       cds_lend=min(g_a, g_b), cds_rend=max(g_a, g_b), model_id=orf.model_id, name=orf.name)


def propagate(orfs_by_transcript, structures: Dict[str, TranscriptStructure],
              transcript_lengths: Dict[str, int]) -> Tuple[List[CodingGene], int, int]:
    """Returns (genes sorted by id, number propagated, total ORFs)."""
    genes: Dict[str, CodingGene] = {}
    n_done = n_total = 0
    for tid in sorted(orfs_by_transcript):
        for orf in orfs_by_transcript[tid]:
            n_total += 1
            struct = structures.get(tid)
            if struct is None:
                logger.warning("%s has no genome structure; skipped", tid)
                continue
            if tid not in transcript_lengths:
                raise ValueError(f"no sequence length for {tid}: provide the transcript FASTA")
            model = place_orf(orf, struct, transcript_lengths[tid])
            if model is None:
                continue
            base_gene = struct.gene_id or re.split(r"~~|::", orf.gene_id)[0]
            # gene ids can repeat at different loci or strands
            gene_key = f"{base_gene}^{model.seqid}^{model.strand}"
            genes.setdefault(gene_key, CodingGene(gene_key)).models.append(model)
            n_done += 1
    return [genes[k] for k in sorted(genes)], n_done, n_total


def run(td_gff3: Path, structures_file: Path, transcripts_fasta: Path, out_gff3: Path,
        structures_format: Optional[str] = None, source: str = "transdecoder",
        unstranded_as: str = "+") -> Tuple[int, int]:
    structures = read_structures(structures_file, structures_format, unstranded_as)
    genes, n_done, n_total = propagate(read_orfs(td_gff3), structures, read_fasta_lengths(transcripts_fasta))
    with open(out_gff3, "w") as out:
        for gene in genes:
            out.write(gene.gff3_text(source))
    logger.info("%d / %d transcript ORFs propagated to the genome", n_done, n_total)
    return n_done, n_total
