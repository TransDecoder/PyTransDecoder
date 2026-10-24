"""
Command-line interface for PyTransDecoder
"""

import sys
import argparse
from pathlib import Path
from . import __version__
from .longorfs import run_longorfs, DEFAULT_MIN_PROTEIN_LENGTH
from .predict import (DEFAULT_START_CONTRAST_ALPHA, DEFAULT_START_KOZAK_WEIGHT,
                      DEFAULT_SHORT_ORF_LENGTH, DEFAULT_SHORT_ORF_MIN_CONTRAST,
                      DEFAULT_SHORT_ORF_RULE)


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
        default=DEFAULT_MIN_PROTEIN_LENGTH,
        help="Minimum protein length in amino acids; ORFs under 100 codons must also pass Predict's "
             f"short-ORF frame-contrast test (default: {DEFAULT_MIN_PROTEIN_LENGTH})"
    )
    parser.add_argument(
        '-G', '--genetic-code',
        default='universal',
        help='Genetic code (universal, Euplotid, Tetrahymena, Candida, etc.) (default: universal)'
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
        '--retain-multiple-orfs', '--retain_multiple_orfs',
        action='store_false',
        dest='single_best_only',
        help='Report every non-overlapping ORF that passes selection, not only the best per '
             'transcript (the default reports one ORF per transcript)'
    )
    parser.add_argument(
        '--single-best-only', '--single_best_only',
        action='store_true',
        dest='single_best_only',
        help=argparse.SUPPRESS  # now the default; accepted for compatibility
    )
    parser.add_argument(
        '--no-refine-starts', '--no_refine_starts',
        action='store_true',
        dest='no_refine_starts',
        help='Skip start codon refinement'
    )
    parser.add_argument(
        '--start-contrast-alpha', type=float, default=DEFAULT_START_CONTRAST_ALPHA,
        dest='start_contrast_alpha',
        help="Start refinement: UTR score per base of the segment before an open ORF's first in-frame "
             "ATG, i.e. the fraction of typical coding signal the segment may show and still be "
             "called UTR. Higher values trim more open 5' ends (suits full-length/long-read "
             "transcripts); lower values keep more 5' partial (suits fragmentary assemblies) "
             f"(default: {DEFAULT_START_CONTRAST_ALPHA})"
    )
    parser.add_argument(
        '--start-kozak-weight', type=float, default=DEFAULT_START_KOZAK_WEIGHT,
        dest='start_kozak_weight',
        help="Start refinement: bases of coding signal the first ATG may overrule per bit of its "
             "Kozak score (log2-odds vs internal ATGs); 0 ignores start context "
             f"(default: {DEFAULT_START_KOZAK_WEIGHT:g})"
    )
    parser.add_argument(
        '--short-orf-length', type=int, default=DEFAULT_SHORT_ORF_LENGTH, dest='short_orf_length',
        help="ORFs with fewer codons than this (from LongOrfs -m below it) need --short-orf-min-contrast "
             f"unless supported by homology (default: {DEFAULT_SHORT_ORF_LENGTH})"
    )
    parser.add_argument(
        '--short-orf-min-contrast', type=float, default=DEFAULT_SHORT_ORF_MIN_CONTRAST,
        dest='short_orf_min_contrast',
        help="Frame contrast (in-frame score minus the better shifted frame) a short ORF needs, in bases "
             "of the training ORFs' median per-base contrast; 0 applies only the usual frame test "
             f"(default: {DEFAULT_SHORT_ORF_MIN_CONTRAST:g})"
    )
    parser.add_argument(
        '--short-orf-rule', choices=['cliff', 'fixed', 'inverse', 'taper', 'curve', 'decoy'], default=DEFAULT_SHORT_ORF_RULE,
        dest='short_orf_rule',
        help="Required frame contrast vs ORF length (experimental): cliff = only below --short-orf-length; "
             "fixed = all lengths; inverse = scaled by short-orf-length/length; taper = falling linearly "
             f"to 0 at the long-ORF auto-accept length (default: {DEFAULT_SHORT_ORF_RULE})"
    )
    parser.add_argument('--short-orf-curve', default=None, dest='short_orf_curve',
                        help="For --short-orf-rule curve: TSV of codons and log prior odds (experimental)")
    parser.add_argument('--decoy-q', type=float, default=0.1, dest='decoy_q',
                        help="For --short-orf-rule decoy: q-value cutoff (experimental; default: %(default)s)")
    parser.add_argument('--decoy-shuffles', type=int, default=0, dest='decoy_shuffles',
                        help="Shuffled copies per transcript for decoys; 0 = auto (~20 Mb of decoys, max 50)")
    parser.add_argument('--decoy-seed', type=int, default=42, dest='decoy_seed', help="Decoy shuffle seed")
    parser.add_argument('--decoy-strata', choices=['len', 'len_tx', 'len_tx_type'], default='len_tx_type',
                        dest='decoy_strata', help="Decoy FDR strata (default: %(default)s)")
    parser.add_argument('--decoy-long-orfs', choices=['none', 'gc', 'decoy'], default='none',
                        dest='decoy_long_orfs',
                        help="Decoy rule: keep ORFs on length alone by none, the GC table, or the 0.999 quantile "
                             "of decoy ORF lengths (default: %(default)s)")
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
            genetic_code=args.genetic_code,
            start_contrast_alpha=args.start_contrast_alpha,
            start_kozak_weight=args.start_kozak_weight,
            short_orf_length=args.short_orf_length,
            short_orf_min_contrast=args.short_orf_min_contrast,
            short_orf_rule=args.short_orf_rule,
            short_orf_curve=args.short_orf_curve,
            decoy_q=args.decoy_q,
            decoy_shuffles=args.decoy_shuffles,
            decoy_seed=args.decoy_seed,
            decoy_strata=args.decoy_strata,
            decoy_long_orfs=args.decoy_long_orfs,
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
    parser.add_argument('--structures', type=Path, default=None, dest='structures',
                        help='Genome structures of the -t transcripts for annotating the genome: GTF, alignment '
                             'GFF3 (GMAP/PASA, Target=), or SAM/BAM spliced alignments (minimap2). Writes '
                             '<transcripts>.transdecoder.genome.gff3')
    parser.add_argument('-m', '--min-protein-length', type=int, default=DEFAULT_MIN_PROTEIN_LENGTH,
                        dest='min_protein_length',
                        help='Minimum protein length in amino acids (default: %(default)s)')
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
    parser.add_argument('--retain-multiple-orfs', action='store_false', dest='single_best_only',
                        help='Report every non-overlapping ORF that passes selection (default: '
                             'one ORF per transcript)')
    parser.add_argument('--single-best-only', action='store_true', dest='single_best_only',
                        help=argparse.SUPPRESS)  # now the default; accepted for compatibility
    parser.add_argument('--no-refine-starts', action='store_true', dest='no_refine_starts',
                        help='Skip start codon refinement')
    parser.add_argument('--start-contrast-alpha', type=float, default=DEFAULT_START_CONTRAST_ALPHA,
                        dest='start_contrast_alpha',
                        help="Start refinement UTR score per base of segment (default: %(default)s)")
    parser.add_argument('--start-kozak-weight', type=float, default=DEFAULT_START_KOZAK_WEIGHT,
                        dest='start_kozak_weight',
                        help="Start refinement bases per bit of Kozak score (default: %(default)s)")
    parser.add_argument('--short-orf-length', type=int, default=DEFAULT_SHORT_ORF_LENGTH,
                        dest='short_orf_length',
                        help="ORFs with fewer codons need --short-orf-min-contrast (default: %(default)s)")
    parser.add_argument('--short-orf-min-contrast', type=float, default=DEFAULT_SHORT_ORF_MIN_CONTRAST,
                        dest='short_orf_min_contrast',
                        help="Frame contrast required of short ORFs, in bases of median coding contrast "
                             "(default: %(default)s)")
    parser.add_argument('--short-orf-rule', choices=['cliff', 'fixed', 'inverse', 'taper', 'curve', 'decoy'],
                        default=DEFAULT_SHORT_ORF_RULE, dest='short_orf_rule',
                        help="Required frame contrast vs ORF length (experimental; default: %(default)s)")
    parser.add_argument('--short-orf-curve', default=None, dest='short_orf_curve',
                        help="For --short-orf-rule curve: TSV of codons and log prior odds (experimental)")
    parser.add_argument('--decoy-q', type=float, default=0.1, dest='decoy_q',
                        help="For --short-orf-rule decoy: q-value cutoff (experimental; default: %(default)s)")
    parser.add_argument('--decoy-shuffles', type=int, default=0, dest='decoy_shuffles',
                        help="Shuffled copies per transcript for decoys; 0 = auto (~20 Mb of decoys, max 50)")
    parser.add_argument('--decoy-seed', type=int, default=42, dest='decoy_seed', help="Decoy shuffle seed")
    parser.add_argument('--decoy-strata', choices=['len', 'len_tx', 'len_tx_type'], default='len_tx_type',
                        dest='decoy_strata', help="Decoy FDR strata (default: %(default)s)")
    parser.add_argument('--decoy-long-orfs', choices=['none', 'gc', 'decoy'], default='none',
                        dest='decoy_long_orfs',
                        help="Decoy rule: keep ORFs on length alone by none, the GC table, or the 0.999 quantile "
                             "of decoy ORF lengths (default: %(default)s)")
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
    if genome_mode and args.structures is not None:
        print("Error: use --structures with -t, not with --genome/--gtf", file=sys.stderr)
        sys.exit(1)
    if args.structures is not None and not args.structures.exists():
        print(f"Error: structures file not found: {args.structures}", file=sys.stderr)
        sys.exit(1)
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
            genetic_code=args.genetic_code,
            start_contrast_alpha=args.start_contrast_alpha,
            start_kozak_weight=args.start_kozak_weight,
            short_orf_length=args.short_orf_length,
            short_orf_min_contrast=args.short_orf_min_contrast,
            short_orf_rule=args.short_orf_rule,
            short_orf_curve=args.short_orf_curve,
            decoy_q=args.decoy_q,
            decoy_shuffles=args.decoy_shuffles,
            decoy_seed=args.decoy_seed,
            decoy_strata=args.decoy_strata,
            decoy_long_orfs=args.decoy_long_orfs,
        )
        predictor.run()
    except Exception as e:
        print(f"Error during Predict phase: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)

    # Phase 3: Propagate ORFs to genome coordinates (with --genome/--gtf or --structures)
    structures = args.gtf if genome_mode else args.structures
    if structures is not None:
        from .genome import run as genome_run
        td_gff3 = output_dir / f"{args.transcripts.name}.transdecoder.gff3"
        genome_gff3 = output_dir / f"{args.transcripts.name}.transdecoder.genome.gff3"
        n_done, n_total = genome_run(td_gff3, structures, args.transcripts, genome_gff3)
        print(f"{n_done} / {n_total} transcript ORFs propagated to the genome; "
              f"genome-coordinate annotations written to: {genome_gff3} "
              f"and {genome_gff3.with_suffix('.bed')}", file=sys.stderr)


def pipeline_main():
    """Entry point for pyTransdecoder pipeline command"""
    parser = create_pipeline_parser()
    args = parser.parse_args()
    pipeline_cmd(args)


def create_genome_map_parser(subparsers):
    """Create parser for genome-map command"""
    parser = subparsers.add_parser(
        'genome-map',
        help='Propagate predicted ORFs to the genome through transcript structures',
        description="""Propagate ORFs predicted on transcripts to genome coordinates, writing gene / mRNA /
exon / CDS / five_prime_UTR / three_prime_UTR models (GFF3).

Transcript structures: GTF (StringTie, LRAA, Cufflinks, GENCODE ...), alignment GFF3 (GMAP, PASA;
Target=<transcript> <start> <end>) or SAM/BAM spliced alignments (minimap2; primary alignments).

Example:
    pytransdecoder genome-map --orfs transcripts.fa.transdecoder.gff3 --structures transcripts.gtf \\
        -t transcripts.fa -o transcripts.fa.transdecoder.genome.gff3
""",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--orfs', required=True, type=Path,
                        help='TransDecoder transcript-level GFF3 (<transcripts>.transdecoder.gff3)')
    parser.add_argument('--structures', required=True, type=Path,
                        help='transcript genome structures: .gtf, alignment .gff3, .sam or .bam')
    parser.add_argument('--format', choices=['gtf', 'gff3', 'sam'], default=None,
                        help='structures format (default: from the file extension)')
    parser.add_argument('-t', '--transcripts', required=True, type=Path, help='transcripts FASTA')
    parser.add_argument('--unstranded-as', choices=['plus', 'minus'], default='plus', dest='unstranded_as',
                        help="strand of GTF transcripts with strand '.' (the strand their sequence was "
                             "extracted from): plus for gffread (default), minus for TransDecoder's "
                             "util/gtf_genome_to_cdna_fasta.pl")
    parser.add_argument('-o', '--output', type=Path, default=None,
                        help='output GFF3 (default: <orfs without .gff3>.genome.gff3); a BED12 file is written '
                             'next to it (.gff3 replaced by .bed)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    parser.set_defaults(func=genome_map_cmd)
    return parser


def genome_map_cmd(args):
    """Execute genome-map command"""
    import logging
    from .genome import run as genome_run
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    for f in (args.orfs, args.structures, args.transcripts):
        if not f.exists():
            print(f"Error: file not found: {f}", file=sys.stderr)
            sys.exit(1)
    out = args.output or Path(str(args.orfs).removesuffix('.gff3') + '.genome.gff3')
    n_done, n_total = genome_run(args.orfs, args.structures, args.transcripts, out, args.format,
                                 unstranded_as='+' if args.unstranded_as == 'plus' else '-')
    print(f"{n_done} / {n_total} transcript ORFs propagated to the genome: {out} "
          f"(BED: {str(out).removesuffix('.gff3') + '.bed'})", file=sys.stderr)


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
    create_genome_map_parser(subparsers)
    
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
