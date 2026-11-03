"""
TransDecoder.LongOrfs - Extract long ORFs from transcript sequences

Phase 1 of TransDecoder pipeline.
"""

import sys
from pathlib import Path
from typing import Optional, Dict
from .core.fasta import read_fasta
from tqdm import tqdm

from .core.orf_finder import ORFFinder
from .core.writers import fasta_entry
from .core.base_freqs import write_base_frequencies

# Minimum ORF length (aa). Predict's decoy q-values (stratified by ORF length) decide how much coding
# evidence short ORFs need.
DEFAULT_MIN_PROTEIN_LENGTH = 50


def run_longorfs(
    transcripts_file: Path,
    min_protein_length: int = DEFAULT_MIN_PROTEIN_LENGTH,
    genetic_code: str = "universal",
    strand_specific: bool = False,
    output_dir: Optional[Path] = None,
    gene_trans_map_file: Optional[Path] = None,
    complete_orfs_only: bool = False,
    verbose: bool = False,
    cpu: int = 1
):
    """
    Extract long ORFs from transcripts (Phase 1).
    
    Args:
        transcripts_file: Input FASTA file
        min_protein_length: Minimum protein length in amino acids
        genetic_code: Genetic code name
        strand_specific: Only analyze top strand
        output_dir: Output directory
        gene_trans_map_file: Gene to transcript mapping file
        complete_orfs_only: Only output complete ORFs
        verbose: Verbose output
        cpu: Worker processes for ORF extraction (output does not depend on it)
    """
    
    # Set up output directory
    if output_dir is None:
        output_dir = Path.cwd()
    output_dir = Path(output_dir)
    
    # Create transdecoder_dir
    workdir = output_dir / f"{transcripts_file.name}.transdecoder_dir"
    workdir.mkdir(parents=True, exist_ok=True)
    
    if verbose:
        print(f"Output directory: {workdir}", file=sys.stderr)
    
    # Load gene-transcript mapping if provided
    gene_trans_map = {}
    if gene_trans_map_file:
        with open(gene_trans_map_file) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    gene_id, trans_id = parts[0], parts[1]
                    gene_trans_map[trans_id] = gene_id
    
    # LongOrfs settings, for Predict's decoy ORFs (same finder settings on shuffled transcripts)
    import json
    with open(workdir / "longorfs.params.json", "w") as fh:
        json.dump({"min_protein_length": min_protein_length, "genetic_code": genetic_code,
                   "strand_specific": bool(strand_specific),
                   "complete_orfs_only": bool(complete_orfs_only)}, fh, indent=1)

    # Output files
    prefix = workdir / "longest_orfs"
    cds_file = prefix.with_suffix('.cds')
    gff3_file = prefix.with_suffix('.gff3')
    pep_file = prefix.with_suffix('.pep')

    finder_args = dict(
        min_protein_length=min_protein_length,
        allow_5prime_partial=not complete_orfs_only,
        allow_3prime_partial=not complete_orfs_only,
        forward_strand=True,
        reverse_strand=not strand_specific,
        genetic_code=genetic_code,
        complete_orfs_only=complete_orfs_only
    )
    settings = (finder_args, gene_trans_map, genetic_code, strand_specific)
    cpu = max(1, int(cpu))

    if verbose:
        print("\n\n- extracting ORFs from transcripts.\n", file=sys.stderr)

    if cpu == 1:
        # one pass over the input, streamed (a counting pass for the progress bar only when verbose)
        total_transcripts = sum(1 for _ in read_fasta(transcripts_file)) if verbose else None
        chunks = _batches(((r.id, str(r.seq)) for r in read_fasta(transcripts_file)), SERIAL_BATCH)
    else:
        records = [(r.id, str(r.seq)) for r in read_fasta(transcripts_file)]
        total_transcripts = len(records)
        # Model IDs are <acc>.p1, .p2, ... numbered within a transcript; IDs of different accessions
        # never coincide (the last '.p' of an ID starts its numeric suffix), so chunks can be numbered
        # independently. A repeated accession continues the numbering of its earlier copy, so inputs
        # with repeated accessions are processed in one chunk.
        n_repeated = total_transcripts - len({acc for acc, _ in records})
        if n_repeated:
            print(f"Warning: {n_repeated} repeated transcript accessions; extracting ORFs in a single "
                  f"process", file=sys.stderr)
            chunks, cpu = [records], 1
        else:
            chunks = _chunks_by_size(records, cpu * 16)
        del records

    if verbose:
        print(f"-total transcripts to examine: {total_transcripts}\n", file=sys.stderr)

    base_counts = {base: 0 for base in 'ACGT'}
    progress = tqdm(total=total_transcripts, desc="Processing transcripts", disable=not verbose,
                    file=sys.stderr)
    with open(gff3_file, 'w') as gff_fh, open(cds_file, 'w') as cds_fh, open(pep_file, 'w') as pep_fh:
        gff_fh.write("##gff-version 3\n")
        for n, counts, gff, cds, pep in _map_chunks(chunks, settings, cpu):
            for base in 'ACGT':
                base_counts[base] += counts[base]
            gff_fh.write(gff)
            cds_fh.write(cds)
            pep_fh.write(pep)
            progress.update(n)
    progress.close()

    # Base frequencies (the per-chunk counts are integers, so their sum is exact)
    write_base_frequencies(base_counts, str(workdir / "base_freqs.dat"))

    # Print summary
    if verbose:
        print(f"\n\n#################################", file=sys.stderr)
        print(f"### Done preparing long ORFs  ###", file=sys.stderr)
        print(f"##################################", file=sys.stderr)
        print(f"\n\tUse file: {pep_file}  for Pfam and/or BlastP searches to enable homology-based coding region identification.\n\n", file=sys.stderr)
        print(f"\tThen, run TransDecoder.Predict\n\n", file=sys.stderr)


# -- ORF extraction of chunks of transcripts (in this process or in worker processes) -------------

_SETTINGS = None
_SEEN_PROT_IDS: set = set()
SERIAL_BATCH = 2000  # transcripts per batch when streaming with --cpu 1


def _init_worker(settings):
    """Settings, and one set of issued model IDs per process for the whole run."""
    global _SETTINGS, _SEEN_PROT_IDS
    _SETTINGS = settings
    _SEEN_PROT_IDS = set()


def _extract_chunk(records):
    """(n transcripts, base counts, gff3 text, cds text, pep text) of an iterable of (acc, sequence),
    numbering model IDs as the original single loop does."""
    finder_args, gene_trans_map, genetic_code, strand_specific = _SETTINGS
    orf_finder = ORFFinder(**finder_args)
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    counts = {base: 0 for base in 'ACGT'}
    gff, cds, pep = [], [], []
    seen_prot_ids = _SEEN_PROT_IDS
    n = 0
    for acc, sequence in records:
        n += 1
        # base counts, as core.base_freqs.compute_base_counts
        useq = sequence.upper()
        for base in 'ACGT':
            k = useq.count(base)
            counts[base] += k
            if not strand_specific:
                counts[complement[base]] += k

        for orf in orf_finder.find_all_orfs(sequence, acc):
            # Generate unique model ID
            pcounter = 1
            model_id = f"{acc}.p{pcounter}"
            while model_id in seen_prot_ids:
                pcounter += 1
                model_id = f"{acc}.p{pcounter}"
            seen_prot_ids.add(model_id)

            gene_id = gene_trans_map.get(acc, f"GENE.{acc}")
            gene_id = f"{gene_id}~~{model_id}"

            orf.model_id = model_id
            orf.gene_id = gene_id
            orf.description = f"ORF type:{orf.orf_type} ({orf.strand})"

            gff.append(orf.to_gff3() + "\n")
            cds.append(fasta_entry(
                f">{model_id} type:{orf.orf_type} {acc}:{orf.start}-{orf.end}({orf.strand})", orf.sequence))
            pep.append(fasta_entry(
                f">{model_id} type:{orf.orf_type} gc:{genetic_code} {acc}:{orf.start}-{orf.end}({orf.strand})",
                orf.protein))
    return n, counts, "".join(gff), "".join(cds), "".join(pep)


def _batches(it, size):
    """Consecutive lists of up to `size` items of an iterator."""
    from itertools import islice
    while True:
        b = list(islice(it, size))
        if not b:
            return
        yield b


def _chunks_by_size(records, n_chunks):
    """Consecutive chunks of roughly equal total sequence length (order preserved)."""
    total = sum(len(s) for _, s in records) or 1
    target = total / n_chunks
    chunks, cur, acc = [], [], 0
    for rec in records:
        cur.append(rec)
        acc += len(rec[1])
        if acc >= target:
            chunks.append(cur)
            cur, acc = [], 0
    if cur:
        chunks.append(cur)
    return chunks


def _map_chunks(chunks, settings, cpu):
    """_extract_chunk results in chunk order."""
    if cpu <= 1:
        _init_worker(settings)
        try:
            for c in chunks:
                yield _extract_chunk(c)
        finally:
            _init_worker(None)  # release the ID set
        return
    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing as mp
    ctx = mp.get_context("fork") if "fork" in mp.get_all_start_methods() else None
    with ProcessPoolExecutor(max_workers=cpu, mp_context=ctx, initializer=_init_worker,
                             initargs=(settings,)) as ex:
        yield from ex.map(_extract_chunk, chunks)
