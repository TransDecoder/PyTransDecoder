"""
TransDecoder.LongOrfs - Extract long ORFs from transcript sequences

Phase 1 of TransDecoder pipeline.
"""

import sys
from pathlib import Path
from typing import Optional, Dict
from collections import defaultdict
from .core.fasta import read_fasta
from tqdm import tqdm

from .core.orf_finder import ORFFinder
from .core.models import ORF
from .core.writers import GFF3Writer, FastaWriter
from .core.base_freqs import compute_base_counts, write_base_frequencies

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
    verbose: bool = False
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

    # Step 1: Compute base frequencies
    if verbose:
        print("\n\n-first extracting base frequencies, we'll need them later.\n", file=sys.stderr)
    
    base_freqs_file = workdir / "base_freqs.dat"
    base_counts = compute_base_counts(str(transcripts_file), strand_specific)
    write_base_frequencies(base_counts, str(base_freqs_file))
    
    # Step 2: Extract ORFs
    if verbose:
        print("\n\n- extracting ORFs from transcripts.\n", file=sys.stderr)
    
    # Initialize ORF finder
    orf_finder = ORFFinder(
        min_protein_length=min_protein_length,
        allow_5prime_partial=not complete_orfs_only,
        allow_3prime_partial=not complete_orfs_only,
        forward_strand=True,
        reverse_strand=not strand_specific,
        genetic_code=genetic_code,
        complete_orfs_only=complete_orfs_only
    )
    
    # Output files
    prefix = workdir / "longest_orfs"
    cds_file = prefix.with_suffix('.cds')
    gff3_file = prefix.with_suffix('.gff3')
    pep_file = prefix.with_suffix('.pep')
    
    # Count total transcripts for progress bar
    total_transcripts = sum(1 for _ in read_fasta(transcripts_file))
    
    if verbose:
        print(f"-total transcripts to examine: {total_transcripts}\n", file=sys.stderr)
    
    # Process transcripts
    model_counter = 0
    trans_counter = 0
    seen_prot_ids = set()
    
    with GFF3Writer(gff3_file) as gff_writer, \
         FastaWriter(cds_file) as cds_writer, \
         FastaWriter(pep_file) as pep_writer:
        
        progress = tqdm(
            read_fasta(transcripts_file),
            total=total_transcripts,
            desc="Processing transcripts",
            disable=not verbose,
            file=sys.stderr
        )
        
        for record in progress:
            trans_counter += 1
            
            acc = record.id
            sequence = str(record.seq)
            
            # Find all ORFs
            orfs = orf_finder.find_all_orfs(sequence, acc)
            
            # Process each ORF
            for orf in orfs:
                model_counter += 1
                
                # Generate unique model ID
                pcounter = 1
                model_id = f"{acc}.p{pcounter}"
                while model_id in seen_prot_ids:
                    pcounter += 1
                    model_id = f"{acc}.p{pcounter}"
                seen_prot_ids.add(model_id)
                
                # Generate gene ID
                if acc in gene_trans_map:
                    gene_id = gene_trans_map[acc]
                else:
                    gene_id = f"GENE.{acc}"
                
                gene_id = f"{gene_id}~~{model_id}"
                
                # Update ORF with IDs
                orf.model_id = model_id
                orf.gene_id = gene_id
                orf.description = f"ORF type:{orf.orf_type} ({orf.strand})"
                
                # Write outputs
                gff_writer.write_orf(orf)
                
                # Write CDS
                cds_header = f">{model_id} type:{orf.orf_type} {acc}:{orf.start}-{orf.end}({orf.strand})"
                cds_writer.write_sequence(cds_header, orf.sequence)
                
                # Write protein
                pep_header = f">{model_id} type:{orf.orf_type} gc:{genetic_code} {acc}:{orf.start}-{orf.end}({orf.strand})"
                pep_writer.write_sequence(pep_header, orf.protein)
    
    # Print summary
    if verbose:
        print(f"\n\n#################################", file=sys.stderr)
        print(f"### Done preparing long ORFs  ###", file=sys.stderr)
        print(f"##################################", file=sys.stderr)
        print(f"\n\tUse file: {pep_file}  for Pfam and/or BlastP searches to enable homology-based coding region identification.\n\n", file=sys.stderr)
        print(f"\tThen, run TransDecoder.Predict\n\n", file=sys.stderr)
