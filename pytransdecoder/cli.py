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
