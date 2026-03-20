"""
Command-line interface for PyTransDecoder
"""

import sys
import argparse
from pathlib import Path
from . import __version__
from .longorfs import run_longorfs


def create_longorfs_parser(subparsers):
    """Create parser for longorfs command"""
    parser = subparsers.add_parser(
        'longorfs',
        help='Extract long ORFs from transcripts (Phase 1)',
        description='''Extract long ORFs from transcripts (Phase 1)
        
This command identifies all potential ORFs in the input transcripts and
outputs them in GFF3 format along with CDS and protein sequences.

Example:
    pytransdecoder longorfs -t transcripts.fasta
        ''',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        '-t', '--transcripts',
        required=True,
        type=Path,
        help='Transcripts FASTA file'
    )
    parser.add_argument(
        '-m', '--min-protein-length',
        type=int,
        default=100,
        help='Minimum protein length in amino acids (default: 100)'
    )
    parser.add_argument(
        '-G', '--genetic-code',
        default='universal',
        help='Genetic code (universal, Euplotes, Tetrahymena, Candida, etc.) (default: universal)'
    )
    parser.add_argument(
        '-S', '--strand-specific',
        action='store_true',
        help='Only analyze top strand'
    )
    parser.add_argument(
        '-O', '--output-dir', '--output_dir',
        dest='output_dir',
        type=Path,
        default=None,
        help='Output directory (default: current directory)'
    )
    parser.add_argument(
        '--gene-trans-map', '--gene_trans_map',
        dest='gene_trans_map',
        type=Path,
        help='Gene-to-transcript mapping file (tab-delimited: gene_id<tab>trans_id)'
    )
    parser.add_argument(
        '--complete-orfs-only', '--complete_orfs_only',
        dest='complete_orfs_only',
        action='store_true',
        help='Only output complete ORFs (with start and stop codons)'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )
    parser.add_argument(
        '--version',
        action='store_true',
        help='Show version and exit'
    )
    
    parser.set_defaults(func=longorfs_cmd)
    return parser


def create_predict_parser(subparsers):
    """Create parser for predict command"""
    parser = subparsers.add_parser(
        'predict',
        help='Predict likely coding regions (Phase 2)',
        description='''Predict likely coding regions (Phase 2)
        
This command takes the output from TransDecoder.LongOrfs and identifies
the most likely coding regions using:
- Hexamer composition scoring
- Homology support (BLAST and Pfam)
- ORF length
- Start codon refinement

Example:
    pytransdecoder predict -t transcripts.fasta --retain_pfam_hits pfam.domtblout
        ''',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        '-t', '--transcripts',
        required=True,
        type=Path,
        help='Transcripts FASTA file'
    )
    parser.add_argument(
        '-O', '--output-dir', '--output_dir',
        dest='output_dir',
        type=Path,
        help='Output directory (default: current directory)'
    )
    parser.add_argument(
        '-T', '--top-orfs-train', '--top_orfs_train',
        dest='top_orfs_train',
        type=int,
        default=500,
        help='Number of top ORFs to use for training Markov model (default: 500)'
    )
    parser.add_argument(
        '--retain-long-orfs-mode', '--retain_long_orfs_mode',
        type=str,
        choices=['dynamic', 'strict'],
        default='dynamic',
        dest='retain_long_orfs_mode',
        help='Mode for retaining long ORFs (default: dynamic)'
    )
    parser.add_argument(
        '--retain-long-orfs-length', '--retain_long_orfs_length',
        type=int,
        default=1000000,
        dest='retain_long_orfs_length',
        help='Under strict mode, minimum length to auto-retain (default: 1000000)'
    )
    parser.add_argument(
        '--retain-pfam-hits', '--retain_pfam_hits',
        type=Path,
        dest='retain_pfam_hits',
        help='Pfam domain hits file from hmmscan'
    )
    parser.add_argument(
        '--retain-blastp-hits', '--retain_blastp_hits',
        type=Path,
        dest='retain_blastp_hits',
        help='BLASTP hits file in outfmt 6 format'
    )
    parser.add_argument(
        '--single-best-only', '--single_best_only',
        action='store_true',
        dest='single_best_only',
        help='Retain only single best ORF per transcript'
    )
    parser.add_argument(
        '--no-refine-starts', '--no_refine_starts',
        action='store_true',
        dest='no_refine_starts',
        help='Skip start codon refinement'
    )
    parser.add_argument(
        '-G', '--genetic-code',
        default='Standard',
        help='Genetic code (default: Standard)'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )
    
    parser.set_defaults(func=predict_cmd)
    return parser


def longorfs_cmd(args):
    """Execute longorfs command"""
    if args.version:
        print(f"TransDecoder.LongOrfs {__version__}")
        sys.exit(0)
    
    # Check that transcripts file exists
    if not args.transcripts.exists():
        print(f"Error: Transcripts file not found: {args.transcripts}", file=sys.stderr)
        sys.exit(1)
    
    # Check gene_trans_map if provided
    if args.gene_trans_map and not args.gene_trans_map.exists():
        print(f"Error: Gene-transcript map file not found: {args.gene_trans_map}", file=sys.stderr)
        sys.exit(1)
    
    try:
        run_longorfs(
            transcripts_file=args.transcripts,
            min_protein_length=args.min_protein_length,
            genetic_code=args.genetic_code,
            strand_specific=args.strand_specific,
            output_dir=args.output_dir,
            gene_trans_map_file=args.gene_trans_map,
            complete_orfs_only=args.complete_orfs_only,
            verbose=args.verbose
        )
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


def predict_cmd(args):
    """Execute predict command"""
    from .predict import TransDecoderPredict
    import logging
    
    # Setup logging
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # Check that transcripts file exists
    if not args.transcripts.exists():
        print(f"Error: Transcripts file not found: {args.transcripts}", file=sys.stderr)
        sys.exit(1)
    
    # Check optional input files
    if args.retain_pfam_hits and not args.retain_pfam_hits.exists():
        print(f"Error: Pfam hits file not found: {args.retain_pfam_hits}", file=sys.stderr)
        sys.exit(1)
    
    if args.retain_blastp_hits and not args.retain_blastp_hits.exists():
        print(f"Error: BLASTP hits file not found: {args.retain_blastp_hits}", file=sys.stderr)
        sys.exit(1)
    
    try:
        predictor = TransDecoderPredict(
            transcripts_file=str(args.transcripts),
            output_dir=str(args.output_dir) if args.output_dir else None,
            top_orfs_train=args.top_orfs_train,
            retain_long_orfs_mode=args.retain_long_orfs_mode,
            retain_long_orfs_length=args.retain_long_orfs_length,
            retain_pfam_hits=str(args.retain_pfam_hits) if args.retain_pfam_hits else None,
            retain_blastp_hits=str(args.retain_blastp_hits) if args.retain_blastp_hits else None,
            single_best_only=args.single_best_only,
            no_refine_starts=args.no_refine_starts,
            genetic_code=args.genetic_code
        )
        predictor.run()
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


def create_pipeline_parser():
    """Create standalone parser for pyTransdecoder pipeline command"""
    parser = argparse.ArgumentParser(
        prog='pyTransdecoder',
        description='pyTransdecoder: Full pipeline - run LongOrfs, optional BLAST/Pfam searches, then Predict',
        epilog='Python port of TransDecoder (https://github.com/TransDecoder/TransDecoder)',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('-t', '--transcripts', required=False, default=None, type=Path,
                        help='Transcripts FASTA file (required unless --genome and --gtf are provided)')
    parser.add_argument('--genome', type=Path, default=None, dest='genome',
                        help='Genome FASTA file (use with --gtf to extract cDNA sequences via gffread)')
    parser.add_argument('--gtf', type=Path, default=None, dest='gtf',
                        help='Genome annotation GTF file (use with --genome to extract cDNA sequences via gffread)')
    parser.add_argument('-m', '--min-protein-length', type=int, default=100, dest='min_protein_length',
                        help='Minimum protein length in amino acids (default: 100)')
    parser.add_argument('-G', '--genetic-code', default='universal', dest='genetic_code',
                        help='Genetic code (default: universal)')
    parser.add_argument('-S', '--strand-specific', action='store_true', dest='strand_specific',
                        help='Only analyze top strand')
    parser.add_argument('-O', '--output-dir', type=Path, default=None, dest='output_dir',
                        help='Output directory (default: current directory)')
    parser.add_argument('--gene-trans-map', type=Path, default=None, dest='gene_trans_map',
                        help='Gene-to-transcript mapping file (tab-delimited: gene_id<tab>trans_id)')
    parser.add_argument('--complete-orfs-only', action='store_true', dest='complete_orfs_only',
                        help='Only output complete ORFs (with start and stop codons)')

    # BLAST args
    parser.add_argument('--blast-search-pep', type=Path, default=None, dest='blast_search_pep',
                        help='Protein FASTA to search against; triggers homology search step '
                             '(database will be built automatically for the selected --blast-tool)')
    parser.add_argument('--blast-tool', choices=['diamond', 'blastp'], default='diamond',
                        dest='blast_tool',
                        help='Homology search tool: diamond or blastp (default: diamond)')
    parser.add_argument('--blast-evalue', type=float, default=1e-5, dest='blast_evalue',
                        help='E-value cutoff for homology search (default: 1e-5)')
    parser.add_argument('--blast-threads', type=int, default=1, dest='blast_threads',
                        help='Number of threads for homology search (default: 1)')
    parser.add_argument('--pfam-search-db', type=Path, default=None, dest='pfam_search_db',
                        help='Pfam HMM database to search with hmmsearch; '
                             'hmmpress will be run automatically if needed')

    # Predict args
    parser.add_argument('-T', '--top-orfs-train', type=int, default=500, dest='top_orfs_train',
                        help='Number of top ORFs to use for training Markov model (default: 500)')
    parser.add_argument('--retain-long-orfs-mode', type=str, choices=['dynamic', 'strict'],
                        default='dynamic', dest='retain_long_orfs_mode',
                        help='Mode for retaining long ORFs (default: dynamic)')
    parser.add_argument('--retain-long-orfs-length', type=int, default=1000000,
                        dest='retain_long_orfs_length',
                        help='Under strict mode, minimum length to auto-retain (default: 1000000)')
    parser.add_argument('--single-best-only', action='store_true', dest='single_best_only',
                        help='Retain only single best ORF per transcript')
    parser.add_argument('--no-refine-starts', action='store_true', dest='no_refine_starts',
                        help='Skip start codon refinement')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Verbose output')
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')

    return parser


def pipeline_cmd(args):
    """Orchestrate the full pyTransdecoder pipeline"""
    import subprocess
    import logging
    from .predict import TransDecoderPredict

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    # Validate inputs and resolve transcripts file
    genome_mode = args.genome is not None or args.gtf is not None
    if genome_mode:
        if args.genome is None or args.gtf is None:
            print("Error: --genome and --gtf must be provided together", file=sys.stderr)
            sys.exit(1)
        if not args.genome.exists():
            print(f"Error: Genome FASTA not found: {args.genome}", file=sys.stderr)
            sys.exit(1)
        if not args.gtf.exists():
            print(f"Error: GTF file not found: {args.gtf}", file=sys.stderr)
            sys.exit(1)
        # Extract cDNA sequences with gffread
        out_dir = args.output_dir if args.output_dir else Path('.')
        out_dir.mkdir(parents=True, exist_ok=True)
        cdna_fasta = out_dir / (args.gtf.stem + ".cDNA.fasta")
        gffread_cmd = ['gffread', str(args.gtf), '-g', str(args.genome), '-w', str(cdna_fasta)]
        print(f"Extracting cDNA sequences: {' '.join(gffread_cmd)}", file=sys.stderr)
        result = subprocess.run(gffread_cmd)
        if result.returncode != 0:
            print(f"Error: gffread failed with return code {result.returncode}", file=sys.stderr)
            sys.exit(result.returncode)
        args.transcripts = cdna_fasta
    elif args.transcripts is None:
        print("Error: provide -t/--transcripts or both --genome and --gtf", file=sys.stderr)
        sys.exit(1)

    if not args.transcripts.exists():
        print(f"Error: Transcripts file not found: {args.transcripts}", file=sys.stderr)
        sys.exit(1)

    if args.blast_search_pep and not args.blast_search_pep.exists():
        print(f"Error: Protein FASTA not found: {args.blast_search_pep}", file=sys.stderr)
        sys.exit(1)

    if args.pfam_search_db and not args.pfam_search_db.exists():
        print(f"Error: Pfam search database not found: {args.pfam_search_db}", file=sys.stderr)
        sys.exit(1)

    if args.gene_trans_map and not args.gene_trans_map.exists():
        print(f"Error: Gene-transcript map file not found: {args.gene_trans_map}", file=sys.stderr)
        sys.exit(1)

    # Phase 1: LongOrfs
    try:
        run_longorfs(
            transcripts_file=args.transcripts,
            min_protein_length=args.min_protein_length,
            genetic_code=args.genetic_code,
            strand_specific=args.strand_specific,
            output_dir=args.output_dir,
            gene_trans_map_file=args.gene_trans_map,
            complete_orfs_only=args.complete_orfs_only,
            verbose=args.verbose
        )
    except Exception as e:
        print(f"Error during LongOrfs phase: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)

    # Derive workdir path (mirrors longorfs.py logic)
    output_dir = args.output_dir if args.output_dir else Path('.')
    workdir = output_dir / f"{args.transcripts.name}.transdecoder_dir"

    # Phase 1.5: Optional homology search (diamond or blastp)
    retain_blastp_hits = None
    if args.blast_search_pep:
        pep_file = workdir / "longest_orfs.pep"
        blast_out = workdir / "blastp.outfmt6"

        # Build the database
        if args.blast_tool == 'diamond':
            db_path = workdir / "blast_db"
            makedb_cmd = [
                'diamond', 'makedb',
                '--in', str(args.blast_search_pep),
                '-d', str(db_path),
                '-p', str(args.blast_threads)
            ]
        else:
            db_path = workdir / "blast_db"
            makedb_cmd = [
                'makeblastdb',
                '-in', str(args.blast_search_pep),
                '-dbtype', 'prot',
                '-out', str(db_path)
            ]
        print(f"Building {args.blast_tool} database: {' '.join(makedb_cmd)}", file=sys.stderr)
        result = subprocess.run(makedb_cmd)
        if result.returncode != 0:
            print(f"Error: database build failed with return code {result.returncode}", file=sys.stderr)
            sys.exit(result.returncode)

        # Run the search
        if args.blast_tool == 'diamond':
            blast_cmd = [
                'diamond', 'blastp',
                '-q', str(pep_file),
                '-d', str(db_path),
                '-k', '1',
                '-f', '6',
                '-e', str(args.blast_evalue),
                '-p', str(args.blast_threads),
                '-o', str(blast_out)
            ]
        else:
            blast_cmd = [
                'blastp',
                '-query', str(pep_file),
                '-db', str(db_path),
                '-max_target_seqs', '1',
                '-outfmt', '6',
                '-evalue', str(args.blast_evalue),
                '-num_threads', str(args.blast_threads),
                '-out', str(blast_out)
            ]
        print(f"Running homology search ({args.blast_tool}): {' '.join(blast_cmd)}", file=sys.stderr)
        result = subprocess.run(blast_cmd)
        if result.returncode != 0:
            print(f"Error: {args.blast_tool} failed with return code {result.returncode}", file=sys.stderr)
            sys.exit(result.returncode)
        retain_blastp_hits = str(blast_out)

    retain_pfam_hits = None
    if args.pfam_search_db:
        pep_file = workdir / "longest_orfs.pep"
        pfam_out = workdir / "pfam.domtblout"
        pfam_db = args.pfam_search_db
        hmmpress_outputs = [Path(str(pfam_db) + ext) for ext in ('.h3f', '.h3i', '.h3m', '.h3p')]

        if not all(path.exists() for path in hmmpress_outputs):
            hmmpress_cmd = ['hmmpress', '-f', str(pfam_db)]
            print(f"Preparing Pfam database: {' '.join(hmmpress_cmd)}", file=sys.stderr)
            result = subprocess.run(hmmpress_cmd)
            if result.returncode != 0:
                print(f"Error: hmmpress failed with return code {result.returncode}", file=sys.stderr)
                sys.exit(result.returncode)

        hmmsearch_cmd = ['hmmsearch', '--domtblout', str(pfam_out), str(pfam_db), str(pep_file)]
        print(f"Running Pfam search: {' '.join(hmmsearch_cmd)}", file=sys.stderr)
        result = subprocess.run(hmmsearch_cmd)
        if result.returncode != 0:
            print(f"Error: hmmsearch failed with return code {result.returncode}", file=sys.stderr)
            sys.exit(result.returncode)
        retain_pfam_hits = str(pfam_out)

    # Phase 2: Predict
    try:
        predictor = TransDecoderPredict(
            transcripts_file=str(args.transcripts),
            output_dir=str(args.output_dir) if args.output_dir else None,
            top_orfs_train=args.top_orfs_train,
            retain_long_orfs_mode=args.retain_long_orfs_mode,
            retain_long_orfs_length=args.retain_long_orfs_length,
            retain_pfam_hits=retain_pfam_hits,
            retain_blastp_hits=retain_blastp_hits,
            single_best_only=args.single_best_only,
            no_refine_starts=args.no_refine_starts,
            genetic_code=args.genetic_code
        )
        predictor.run()
    except Exception as e:
        print(f"Error during Predict phase: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)

    # Phase 3: Propagate ORFs to genome coordinates (only when genome+gtf provided)
    if genome_mode:
        util_dir = Path(__file__).parent.parent / 'util'
        gtf_to_gff3_pl = util_dir / 'gtf_to_alignment_gff3.pl'
        cdna_to_genome_pl = util_dir / 'cdna_alignment_orf_to_genome_orf.pl'

        # Convert GTF to alignment GFF3
        alignment_gff3 = output_dir / (args.gtf.stem + ".gff3")
        cmd = ['perl', str(gtf_to_gff3_pl), str(args.gtf)]
        print(f"Converting GTF to alignment GFF3: {' '.join(cmd)}", file=sys.stderr)
        with open(alignment_gff3, 'w') as fh:
            result = subprocess.run(cmd, stdout=fh)
        if result.returncode != 0:
            print(f"Error: gtf_to_alignment_gff3.pl failed with return code {result.returncode}", file=sys.stderr)
            sys.exit(result.returncode)

        # Map transcript ORF coords to genome coords
        td_gff3 = output_dir / f"{args.transcripts.name}.transdecoder.gff3"
        genome_gff3 = output_dir / f"{args.transcripts.name}.transdecoder.genome.gff3"
        cmd = ['perl', str(cdna_to_genome_pl), str(td_gff3), str(alignment_gff3), str(args.transcripts)]
        print(f"Propagating ORFs to genome coordinates: {' '.join(cmd)}", file=sys.stderr)
        with open(genome_gff3, 'w') as fh:
            result = subprocess.run(cmd, stdout=fh)
        if result.returncode != 0:
            print(f"Error: cdna_alignment_orf_to_genome_orf.pl failed with return code {result.returncode}", file=sys.stderr)
            sys.exit(result.returncode)
        print(f"Genome-coordinate ORF annotations written to: {genome_gff3}", file=sys.stderr)


def pipeline_main():
    """Entry point for pyTransdecoder pipeline command"""
    parser = create_pipeline_parser()
    args = parser.parse_args()
    pipeline_cmd(args)


def main():
    """Main entry point for CLI"""
    parser = argparse.ArgumentParser(
        prog='pytransdecoder',
        description='PyTransDecoder: Identify candidate coding regions within transcript sequences',
        epilog='Python port of TransDecoder (https://github.com/TransDecoder/TransDecoder)'
    )
    
    parser.add_argument(
        '--version',
        action='version',
        version=f'%(prog)s {__version__}'
    )
    
    subparsers = parser.add_subparsers(
        title='commands',
        description='Valid commands',
        help='Use "pytransdecoder <command> -h" for command-specific help',
        dest='command'
    )
    
    # Create subcommands
    create_longorfs_parser(subparsers)
    create_predict_parser(subparsers)
    
    # Parse arguments
    args = parser.parse_args()
    
    # Execute command
    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()
        sys.exit(1)


def cli():
    """Alias for main() to maintain compatibility"""
    main()


if __name__ == '__main__':
    main()
