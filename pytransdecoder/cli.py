"""
Command-line interface for PyTransDecoder
"""

import sys
import click
from pathlib import Path
from . import __version__
from .longorfs import run_longorfs
from .predict import main as predict_main


@click.group()
@click.version_option(version=__version__)
def cli():
    """
    PyTransDecoder: Identify candidate coding regions within transcript sequences
    
    Python port of TransDecoder (https://github.com/TransDecoder/TransDecoder)
    """
    pass


@cli.command('longorfs')
@click.option(
    '-t', '--transcripts',
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help='Transcripts FASTA file'
)
@click.option(
    '-m', '--min-protein-length',
    default=100,
    type=int,
    show_default=True,
    help='Minimum protein length in amino acids'
)
@click.option(
    '-G', '--genetic-code',
    default='universal',
    show_default=True,
    help='Genetic code (universal, Euplotes, Tetrahymena, Candida, etc.)'
)
@click.option(
    '-S', '--strand-specific',
    is_flag=True,
    help='Only analyze top strand'
)
@click.option(
    '-O', '--output-dir',
    type=click.Path(path_type=Path),
    default=None,
    help='Output directory (default: current directory)'
)
@click.option(
    '--gene-trans-map',
    type=click.Path(exists=True, path_type=Path),
    help='Gene-to-transcript mapping file (tab-delimited: gene_id<tab>trans_id)'
)
@click.option(
    '--complete-orfs-only',
    is_flag=True,
    help='Only output complete ORFs (with start and stop codons)'
)
@click.option(
    '-v', '--verbose',
    is_flag=True,
    help='Verbose output'
)
@click.option(
    '--version',
    is_flag=True,
    help='Show version and exit'
)
def longorfs_cmd(
    transcripts,
    min_protein_length,
    genetic_code,
    strand_specific,
    output_dir,
    gene_trans_map,
    complete_orfs_only,
    verbose,
    version
):
    """
    Extract long ORFs from transcripts (Phase 1)
    
    This command identifies all potential ORFs in the input transcripts and
    outputs them in GFF3 format along with CDS and protein sequences.
    
    Example:
    
        pytransdecoder longorfs -t transcripts.fasta
    """
    if version:
        click.echo(f"TransDecoder.LongOrfs {__version__}")
        sys.exit(0)
    
    try:
        run_longorfs(
            transcripts_file=transcripts,
            min_protein_length=min_protein_length,
            genetic_code=genetic_code,
            strand_specific=strand_specific,
            output_dir=output_dir,
            gene_trans_map_file=gene_trans_map,
            complete_orfs_only=complete_orfs_only,
            verbose=verbose
        )
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        if verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


@cli.command('predict')
@click.option('-t', '--transcripts', 'transcripts_file', required=True, type=click.Path(exists=True),
              help='Transcripts FASTA file')
@click.option('-O', '--output-dir', type=click.Path(),
              help='Output directory (default: current directory)')
@click.option('-T', '--top-orfs-train', type=int, default=500,
              help='Number of top ORFs to use for training Markov model (default: 500)')
@click.option('--retain-long-orfs-mode', type=click.Choice(['dynamic', 'strict']), default='dynamic',
              help='Mode for retaining long ORFs (default: dynamic)')
@click.option('--retain-long-orfs-length', type=int, default=1000000,
              help='Under strict mode, minimum length to auto-retain (default: 1000000)')
@click.option('--retain-pfam-hits', type=click.Path(exists=True),
              help='Pfam domain hits file from hmmscan')
@click.option('--retain-blastp-hits', type=click.Path(exists=True),
              help='BLASTP hits file in outfmt 6 format')
@click.option('--single-best-only', is_flag=True,
              help='Retain only single best ORF per transcript')
@click.option('--no-refine-starts', is_flag=True,
              help='Skip start codon refinement')
@click.option('-G', '--genetic-code', default='Standard',
              help='Genetic code (default: Standard)')
@click.option('-v', '--verbose', is_flag=True,
              help='Verbose output')
def predict_cmd(
    transcripts_file,
    output_dir,
    top_orfs_train,
    retain_long_orfs_mode,
    retain_long_orfs_length,
    retain_pfam_hits,
    retain_blastp_hits,
    single_best_only,
    no_refine_starts,
    genetic_code,
    verbose
):
    """
    Predict likely coding regions (Phase 2)
    
    This command takes the output from TransDecoder.LongOrfs and identifies
    the most likely coding regions using:
    - Hexamer composition scoring
    - Homology support (BLAST and Pfam)
    - ORF length
    - Start codon refinement
    """
    from .predict import TransDecoderPredict
    import logging
    
    # Setup logging
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    try:
        predictor = TransDecoderPredict(
            transcripts_file=transcripts_file,
            output_dir=output_dir,
            top_orfs_train=top_orfs_train,
            retain_long_orfs_mode=retain_long_orfs_mode,
            retain_long_orfs_length=retain_long_orfs_length,
            retain_pfam_hits=retain_pfam_hits,
            retain_blastp_hits=retain_blastp_hits,
            single_best_only=single_best_only,
            no_refine_starts=no_refine_starts,
            genetic_code=genetic_code
        )
        predictor.run()
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        if verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    cli()
