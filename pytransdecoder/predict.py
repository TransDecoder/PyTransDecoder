"""
TransDecoder.Predict - Final coding region predictions

This module implements the prediction phase of TransDecoder, which:
1. Trains a Markov model on longest ORFs
2. Scores all candidate ORFs
3. Incorporates homology information (BLAST/Pfam)
4. Selects best ORFs per transcript
5. Refines start codons
6. Generates final output files
"""

import logging
import math
import re
from pathlib import Path
from typing import Optional, Dict, List, Tuple
from Bio import SeqIO

from .core.fasta import read_fasta

from pytransdecoder.core.base_freqs import read_base_frequencies

logger = logging.getLogger(__name__)

# Start refinement: an ORF open at its 5' end stays 5'-partial only when the segment before its
# first in-frame ATG has at least this fraction of the training ORFs' per-base frame contrast
# (TransDecoderPaper/benchmark_v2; chosen on the 2026 benchmark, alpha 0.4).
DEFAULT_START_CONTRAST_ALPHA = 0.4

# Kozak log-odds (nats) of the candidate ATG add to the frame-contrast cutoff, independent of the
# segment's length: an open 5' end of L bases is trimmed when its contrast <
# (alpha * L + weight * kozak) * the training ORFs' per-base contrast. A weight of 25 equates
# one nat of Kozak evidence with 25 bases of the median coding contrast (benchmark_v2 and GM12878
# MANE starts, 2026).
DEFAULT_START_KOZAK_WEIGHT = 25.0

# ORFs shorter than this many codons (LongOrfs -m below 100) are kept only with homology support or
# when their frame contrast (in-frame score minus the better shifted frame) reaches
# DEFAULT_SHORT_ORF_MIN_CONTRAST bases' worth of the training ORFs' median per-base contrast. Short
# ORFs in noncoding RNA often pass the plain frame test; a fixed total contrast, rather than a
# per-base one, asks for more coding evidence the shorter the ORF (benchmark_v2 lncRNA and random
# negatives, 2026).
DEFAULT_SHORT_ORF_LENGTH = 100
DEFAULT_SHORT_ORF_MIN_CONTRAST = 200.0

# Stop codons the GC-dependent long-ORF length table was simulated with
STANDARD_STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})


def _perl_reverse_sort(items, key):
    """Longest-first ordering exactly as Perl's `reverse sort {$a <=> $b}` produces it.

    Perl's sort is stable, so reversing an ascending sort places tied items in
    reverse input order (unlike sorted(..., reverse=True), which keeps input order).
    The training-set selection applies this ordering at each of its three steps,
    and tie order determines which equal-length ORFs make the training set.
    """
    return list(reversed(sorted(items, key=key)))


class TransDecoderPredict:
    """Main class for TransDecoder prediction phase"""
    
    def __init__(
        self,
        transcripts_file: Path,
        output_dir: Optional[Path] = None,
        top_orfs_train: int = 500,
        retain_long_orfs_mode: str = 'dynamic',
        retain_long_orfs_length: int = 1000000,
        retain_pfam_hits: Optional[Path] = None,
        retain_blastp_hits: Optional[Path] = None,
        single_best_only: bool = True,
        no_refine_starts: bool = False,
        genetic_code: str = "Standard",
        start_contrast_alpha: float = DEFAULT_START_CONTRAST_ALPHA,
        start_kozak_weight: float = DEFAULT_START_KOZAK_WEIGHT,
        short_orf_length: int = DEFAULT_SHORT_ORF_LENGTH,
        short_orf_min_contrast: float = DEFAULT_SHORT_ORF_MIN_CONTRAST,
    ):
        """
        Initialize TransDecoder.Predict
        
        Args:
            transcripts_file: Input transcripts FASTA file
            output_dir: Output directory (defaults to current directory)
            top_orfs_train: Number of top ORFs to use for training (default: 500)
            retain_long_orfs_mode: 'dynamic' or 'strict' (default: dynamic)
            retain_long_orfs_length: Min length to auto-retain under strict mode
            retain_pfam_hits: Pfam domain hits file
            retain_blastp_hits: BLASTP hits file  
            single_best_only: Retain only the best ORF per transcript (default); False keeps
                every non-overlapping ORF that passes selection
            no_refine_starts: Skip start codon refinement
            genetic_code: Genetic code to use
            start_contrast_alpha: Start refinement keeps an open 5' end only when the segment
                before the first in-frame ATG has at least this fraction of the training ORFs'
                per-base frame contrast
            start_kozak_weight: bases of median per-base frame contrast that one unit of the
                first ATG's Kozak log-odds is worth in that decision (0 uses frame contrast alone)
            short_orf_length: ORFs with fewer codons than this need short_orf_min_contrast
            short_orf_min_contrast: frame contrast, in bases of the training ORFs' median per-base
                contrast, that an ORF shorter than short_orf_length needs without homology support
                (0 applies only the usual frame test)
        """
        self.transcripts_file = Path(transcripts_file)
        self.output_dir = Path(output_dir) if output_dir else Path.cwd()
        self.top_orfs_train = top_orfs_train
        self.retain_long_orfs_mode = retain_long_orfs_mode
        self.retain_long_orfs_length = retain_long_orfs_length
        self.retain_pfam_hits = Path(retain_pfam_hits) if retain_pfam_hits else None
        self.retain_blastp_hits = Path(retain_blastp_hits) if retain_blastp_hits else None
        self.single_best_only = single_best_only
        self.no_refine_starts = no_refine_starts
        self.genetic_code = genetic_code
        self.start_contrast_alpha = start_contrast_alpha
        self.start_kozak_weight = start_kozak_weight
        self.short_orf_length = short_orf_length
        self.short_orf_min_contrast = short_orf_min_contrast
        self._coding_contrast_cache = None
        # Terminal stops excluded from Markov training and scoring
        from pytransdecoder.core.translator import Translator
        self.stop_codons = frozenset(Translator(genetic_code).get_stop_codons())
        
        # Derived paths
        self.workdir = self.output_dir / f"{self.transcripts_file.name}.transdecoder_dir"
        self.checkpoints_dir = self.workdir / "__checkpoints_TDpredict"
        
        # Input files from LongOrfs phase
        self.longest_orfs_prefix = self.workdir / "longest_orfs"
        self.cds_file = self.longest_orfs_prefix.with_suffix(self.longest_orfs_prefix.suffix + ".cds")
        self.gff3_file = self.longest_orfs_prefix.with_suffix(self.longest_orfs_prefix.suffix + ".gff3")
        self.pep_file = self.longest_orfs_prefix.with_suffix(self.longest_orfs_prefix.suffix + ".pep")
        self.base_freqs_file = self.workdir / "base_freqs.dat"
        
        # GC% -> minimum ORF length (nt) retained regardless of coding score: the 0.999
        # quantile of LongOrfs ORF lengths in random sequence of that GC content
        # (TransDecoderPaper/analyses/long_orf_threshold; same table as TransDecoder.Predict)
        self.gc_to_min_length = [
            (25, 519), (30, 552), (35, 573), (40, 615),
            (45, 663), (50, 726), (55, 804), (60, 927),
            (65, 1083), (70, 1347), (75, 1746), (80, 2436)
        ]
        
    def validate_inputs(self):
        """Validate that required input files exist"""
        if not self.transcripts_file.exists():
            raise FileNotFoundError(f"Transcripts file not found: {self.transcripts_file}")
        
        if not self.workdir.exists():
            raise FileNotFoundError(
                f"Working directory not found: {self.workdir}\n"
                "Please run TransDecoder.LongOrfs first!"
            )
        
        for required_file in [self.cds_file, self.gff3_file, self.pep_file]:
            if not required_file.exists():
                raise FileNotFoundError(
                    f"Required file not found: {required_file}\n"
                    "Please run TransDecoder.LongOrfs first!"
                )
    
    def run(self):
        """Run the complete prediction pipeline"""
        logger.info("Starting TransDecoder.Predict")
        
        # Validate inputs
        self.validate_inputs()
        
        # Create checkpoints directory
        self.checkpoints_dir.mkdir(exist_ok=True)
        
        # Step 1: Get top longest ORFs for training
        logger.info(f"Selecting top {self.top_orfs_train} longest ORFs for training...")
        top_cds_file = self._select_training_orfs()
        
        # Step 2: Calculate hexamer scores (Markov model)
        logger.info("Training hexamer scoring model...")
        hexamer_scores_file = self._train_hexamer_model(top_cds_file)
        
        # Step 3: Score all ORFs
        logger.info("Scoring all candidate ORFs...")
        cds_scores_file = self._score_all_orfs(hexamer_scores_file)
        
        # Step 4: Select best ORFs
        logger.info("Selecting best ORFs per transcript...")
        best_orfs_gff3 = self._select_best_orfs(cds_scores_file, top_cds_file, hexamer_scores_file)
        
        # Step 5: Refine start codons (if enabled)
        if not self.no_refine_starts:
            logger.info("Refining start codon predictions...")
            best_orfs_gff3 = self._refine_start_codons(best_orfs_gff3, top_cds_file, hexamer_scores_file)
        
        # Step 6: Generate final outputs
        logger.info("Generating final output files...")
        self._generate_final_outputs(best_orfs_gff3)
        
        logger.info("TransDecoder.Predict completed successfully!")
        
    def _select_training_orfs(self) -> Path:
        """
        Select top longest ORFs for training, removing redundancy
        
        Process:
        1. Get top N*10 longest ORFs (up to max protein length)
        2. Remove redundant/similar sequences
        3. Select top N from the non-redundant set
        
        Returns:
            Path to file containing training ORF sequences
        """
        checkpoint = self.checkpoints_dir / "training_orfs.ok"
        output_file = self.checkpoints_dir / "top_training_orfs.cds"
        
        if checkpoint.exists():
            logger.info("Training ORFs already selected (checkpoint exists)")
            return output_file
        
        # Parameters from Perl code
        red_num = self.top_orfs_train * 10
        max_cds_length = 5000  # From Perl: my $max_prot_length_for_training = 5000 (CDS nucleotide length, not protein aa)
        
        # Step 1: Get top N*10 longest ORFs (up to max CDS length)
        logger.info(f"Getting top {red_num} longest ORFs (max CDS length: {max_cds_length} nt)...")
        longest_orfs = []
        
        with open(self.cds_file) as f:
            for record in read_fasta(f):
                seq_len = len(record.seq)
                
                # Filter by CDS nucleotide length (not protein length)
                if seq_len <= max_cds_length:
                    longest_orfs.append((record, seq_len))
        
        # Sort by length (longest first)
        longest_orfs = _perl_reverse_sort(longest_orfs, key=lambda x: x[1])
        
        # Take top red_num
        longest_orfs = longest_orfs[:red_num]
        logger.info(f"Selected {len(longest_orfs)} ORFs for redundancy filtering")
        
        # Step 2: Remove redundant sequences
        logger.info("Removing redundant sequences...")
        nr_orfs = self._exclude_similar_proteins(longest_orfs)
        logger.info(f"After redundancy filtering: {len(nr_orfs)} ORFs retained "
                   f"({len(nr_orfs)/len(longest_orfs)*100:.1f}%)")
        
        # Step 3: Select top N from non-redundant set
        nr_orfs_sorted = _perl_reverse_sort(nr_orfs, key=lambda x: len(x.seq))
        training_orfs = nr_orfs_sorted[:self.top_orfs_train]
        logger.info(f"Final training set: {len(training_orfs)} ORFs")
        
        # Write output
        with open(output_file, 'w') as f:
            SeqIO.write(training_orfs, f, "fasta")
        
        checkpoint.write_text("OK")
        return output_file
    
    def _exclude_similar_proteins(self, orf_records: List[Tuple]) -> List:
        """
        Exclude similar proteins using k-mer profiling
        
        Implementation of exclude_similar_proteins.pl logic:
        - Uses 5-mer profiles of translated sequences
        - Filters low complexity sequences (< 30% unique kmers)
        - Excludes sequences with > 80% similarity to already selected sequences
        
        Args:
            orf_records: List of (SeqRecord, length) tuples
            
        Returns:
            List of non-redundant SeqRecord objects
        """
        from pytransdecoder.core.translator import Translator
        
        WMER_SIZE = 5
        MAX_PER_ID = 80
        captured_wmers = {}  # wmer -> list of protein accessions
        retained_orfs = []
        
        # Create translator
        translator = Translator(self.genetic_code)
        
        # exclude_similar_proteins.pl re-sorts its input longest-first, which reverses ties again
        for record, _ in _perl_reverse_sort(orf_records, key=lambda x: x[1]):
            acc = record.id
            sequence = str(record.seq)
            
            # Translate the CDS sequence
            protein = translator.translate(sequence, frame=0)
            
            # Get spaced 5-mers for this protein (spaced by WMER_SIZE)
            wmers_spaced = {}
            for i in range(0, len(protein) - WMER_SIZE + 1, WMER_SIZE):
                wmer = protein[i:i + WMER_SIZE]
                wmers_spaced[wmer] = 1
            
            num_unique_wmers = len(wmers_spaced)
            if num_unique_wmers == 0:
                continue
                
            pct_unique_wmers = num_unique_wmers / (len(protein) / WMER_SIZE) * 100
            
            # Skip low complexity sequences
            if pct_unique_wmers < 30:
                logger.debug(f"Skipping {acc} as likely low complexity sequence")
                continue
            
            # Check if this protein is too similar to already captured ones
            prots_seen = {}
            for wmer in wmers_spaced:
                if wmer in captured_wmers:
                    for prot_acc in captured_wmers[wmer]:
                        prots_seen[prot_acc] = prots_seen.get(prot_acc, 0) + 1
            
            proxy_per_id = 0
            if prots_seen:
                max_seen_acc = max(prots_seen, key=prots_seen.get)
                max_seen_count = prots_seen[max_seen_acc]
                proxy_per_id = max_seen_count / num_unique_wmers * 100
                logger.debug(f"Analyzing {acc}: max_seen={max_seen_count} in {max_seen_acc}, "
                           f"proxy_per_id={proxy_per_id:.2f}")
            
            if proxy_per_id <= MAX_PER_ID:
                # Keep this ORF - it's unique enough
                
                # Recompute wmers as overlapping (not spaced) for storage
                wmers_overlapping = {}
                for i in range(len(protein) - WMER_SIZE + 1):
                    wmer = protein[i:i + WMER_SIZE]
                    wmers_overlapping[wmer] = 1
                
                # Store wmers for this protein
                for wmer in wmers_overlapping:
                    if wmer not in captured_wmers:
                        captured_wmers[wmer] = []
                    captured_wmers[wmer].append(acc)
                
                retained_orfs.append(record)
            else:
                logger.debug(f"Skipping training candidate: {acc}, not unique enough "
                           f"(proxy_per_id={proxy_per_id:.2f})")
        
        return retained_orfs
    
    def _train_hexamer_model(self, training_orfs_file: Path) -> Path:
        """
        Train hexamer scoring model on training ORFs
        
        This implements a Markov chain model that computes log-likelihood scores
        for hexamers (6-mers) in each reading frame.
        
        The algorithm:
        1. Count k-mers (k=1 to 6) in each frame of training ORFs
        2. Calculate Markov probabilities: P(base | k-1mer, frame)
        3. Compare to background base frequencies
        4. Output log-likelihood ratios
        
        Args:
            training_orfs_file: File containing training ORF sequences
            
        Returns:
            Path to hexamer scores file
        """
        checkpoint = self.checkpoints_dir / "hexamer_scores.ok"
        output_file = self.workdir / "hexamer.scores"
        
        if checkpoint.exists():
            logger.info("Hexamer model already trained (checkpoint exists)")
            return output_file
        
        # Parse training ORFs and count framed k-mers
        logger.info("Counting framed k-mers in training ORFs...")
        framed_kmers = {}  # key: "kmer-frame", value: count
        
        with open(training_orfs_file) as f:
            for record in read_fasta(f):
                sequence = str(record.seq).upper()
                seq_len = len(sequence)
                
                # Count k-mers for Markov orders 0-5 (producing 1-mers to 6-mers)
                for markov_order in range(6):
                    for i in range(markov_order, seq_len):
                        frame = i % 3
                        
                        # Skip stop codons at the end
                        if (i == seq_len - 3 and frame == 0):
                            codon = sequence[i:i+3]
                            if codon in self.stop_codons:
                                break
                        
                        # Count frame positions (needed for Markov chain)
                        if markov_order == 0:
                            key = f"FRAME-{frame}"
                            framed_kmers[key] = framed_kmers.get(key, 0) + 1
                        
                        # Count k-mer in this frame
                        if i >= markov_order:
                            kmer = sequence[i - markov_order:i + 1]
                            key = f"{kmer}-{frame}"
                            framed_kmers[key] = framed_kmers.get(key, 0) + 1
        
        logger.info(f"Counted {len(framed_kmers)} unique framed k-mers")
        
        # Load background base probabilities
        logger.info("Loading background base frequencies...")
        if not self.base_freqs_file.exists():
            raise FileNotFoundError(f"{self.base_freqs_file} not found; run LongOrfs first")
        background_probs = read_base_frequencies(str(self.base_freqs_file))
        logger.info(f"Background base frequencies: {background_probs}")

        # Calculate log-likelihood ratios
        logger.info("Computing log-likelihood scores...")
        scores = []
        
        for framed_kmer, kmer_count in sorted(framed_kmers.items()):
            if '-' not in framed_kmer or framed_kmer.startswith('FRAME-'):
                continue
                
            kmer, frame_str = framed_kmer.rsplit('-', 1)
            frame = int(frame_str)
            
            # Skip k-mers with non-GATC bases
            if any(c not in 'GATC' for c in kmer):
                continue
            
            kmer_length = len(kmer)
            
            # Get k-1mer count in previous frame
            kminus1_frame = (frame - 1) % 3
            
            if kmer_length > 1:
                kminus1mer = kmer[:-1]
                kminus1_key = f"{kminus1mer}-{kminus1_frame}"
                kminus1_count = framed_kmers.get(kminus1_key, 0)
            else:
                kminus1_key = f"FRAME-{kminus1_frame}"
                kminus1_count = framed_kmers.get(kminus1_key, 0)
            
            # Markov probability with pseudocounts
            markov_prob = (kmer_count + 1) / (kminus1_count + 4)
            
            # Background probability
            last_base = kmer[-1]
            background_prob = background_probs[last_base]
            
            # Log-likelihood ratio
            import math
            loglikelihood = math.log(markov_prob / background_prob)
            
            scores.append((framed_kmer, kmer_count, kminus1_count, loglikelihood))
        
        # Write scores file
        logger.info(f"Writing {len(scores)} hexamer scores to {output_file}")
        with open(output_file, 'w') as f:
            f.write("#framed_kmer\tkmer_count\tkminus1_prefix_count\tloglikelihood\n")
            for framed_kmer, kmer_count, kminus1_count, loglikelihood in scores:
                f.write(f"{framed_kmer}\t{kmer_count}\t{kminus1_count}\t{loglikelihood}\n")
        
        checkpoint.write_text("OK")
        return output_file
    
    def _score_all_orfs(self, hexamer_scores_file: Path) -> Path:
        """
        Score all ORFs using hexamer model
        
        Scores each ORF in all 6 reading frames using the trained hexamer model.
        The score for each frame is the sum of log-likelihood scores for each
        hexamer in that frame.
        
        Args:
            hexamer_scores_file: File containing hexamer log-likelihood scores
            
        Returns:
            Path to CDS scores file
        """
        checkpoint = self.checkpoints_dir / "cds_scores.ok"
        output_file = self.workdir / "longest_orfs.cds.scores"
        
        if checkpoint.exists():
            logger.info("ORF scores already computed (checkpoint exists)")
            return output_file
        
        # Load hexamer scores
        logger.info("Loading hexamer scores...")
        scores = {}
        with open(hexamer_scores_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    framed_kmer = parts[0]
                    loglikelihood = float(parts[3])
                    scores[framed_kmer] = loglikelihood
        
        logger.info(f"Loaded {len(scores)} hexamer scores")
        
        # Score all ORFs in all 6 frames
        logger.info("Scoring all candidate ORFs...")
        from pytransdecoder.core.sequence import reverse_complement
        
        results = []
        with open(self.cds_file) as f:
            for record in read_fasta(f):
                accession = record.id
                sequence = str(record.seq).upper()
                
                # Use Markov order 5 (hexamers = 6-mers)
                markov_order = 5
                
                # Score forward strand frames
                score1 = self._score_cds_via_markov(sequence, markov_order, scores)
                score2 = self._score_cds_via_markov(sequence[1:], markov_order, scores)
                score3 = self._score_cds_via_markov(sequence[2:], markov_order, scores)
                
                # Score reverse strand frames
                rev_seq = reverse_complement(sequence)
                score4 = self._score_cds_via_markov(rev_seq, markov_order, scores)
                score5 = self._score_cds_via_markov(rev_seq[1:], markov_order, scores)
                score6 = self._score_cds_via_markov(rev_seq[2:], markov_order, scores)
                
                results.append((accession, markov_order, len(sequence),
                              score1, score2, score3, score4, score5, score6))
        
        # Write scores file
        logger.info(f"Writing scores for {len(results)} ORFs to {output_file}")
        with open(output_file, 'w') as f:
            f.write("#acc\tMarkov_order\tseq_length\tscore_1\tscore_2\tscore_3\tscore_4\tscore_5\tscore_6\n")
            for row in results:
                acc, order, length, s1, s2, s3, s4, s5, s6 = row
                f.write(f"{acc}\t{order}\t{length}\t{s1:.2f}\t{s2:.2f}\t{s3:.2f}\t{s4:.2f}\t{s5:.2f}\t{s6:.2f}\n")
        
        checkpoint.write_text("OK")
        return output_file
    
    def _score_cds_via_markov(self, sequence: str, markov_order: int, scores: Dict[str, float]) -> float:
        """
        Score a sequence using Markov model
        
        Args:
            sequence: DNA sequence to score
            markov_order: Markov order (5 for hexamers)
            scores: Dictionary of framed k-mer scores
            
        Returns:
            Total log-likelihood score
        """
        seq_length = len(sequence)
        
        if seq_length < markov_order + 1:
            return 0.0
        
        total_score = 0.0
        
        for i in range(seq_length):
            frame = i % 3
            
            # Use appropriate k-mer size (up to markov_order)
            markov_use = min(i, markov_order)
            kmer = sequence[i - markov_use:i + 1]
            
            # Skip stop codons at the end
            if i == seq_length - 3 and frame == 0:
                codon = sequence[i:i+3]
                if codon in self.stop_codons:
                    break
            
            # Look up score
            framed_kmer = f"{kmer}-{frame}"
            loglikelihood = scores.get(framed_kmer, 0.0)
            
            total_score += loglikelihood
        
        return total_score
    
    def _select_best_orfs(self, cds_scores_file: Path, training_orfs_file: Optional[Path] = None,
                          hexamer_scores_file: Optional[Path] = None) -> Path:
        """
        Select best ORFs per transcript based on scores and homology
        
        Selection criteria (any of):
        1. Has BLAST hit
        2. Has Pfam hit
        3. Score[frame0] > 0 AND score[frame0] > max(score[frame1], score[frame2]); for ORFs
           shorter than short_orf_length codons, score[frame0] - max(score[frame1], score[frame2])
           must also reach short_orf_min_contrast * the training ORFs' median per-base contrast
        4. Length >= min_length_auto_accept (for dynamic/strict mode)
        
        Prioritization:
        - homology_count (BLAST + Pfam hits)
        - frame score[0]
        - ORF length
        
        Args:
            cds_scores_file: File containing CDS scores
            
        Returns:
            Path to GFF3 file with selected ORFs
        """
        checkpoint = self.checkpoints_dir / "best_orfs.ok"
        output_file = self.workdir / "best_candidates.gff3"
        
        if checkpoint.exists():
            logger.info("Best ORFs already selected (checkpoint exists)")
            return output_file
        
        # Parse homology data
        blast_hits = {}
        if self.retain_blastp_hits:
            logger.info(f"Parsing BLAST hits from {self.retain_blastp_hits}")
            blast_hits = self._parse_blastp_hits()
        
        pfam_hits = {}
        if self.retain_pfam_hits:
            logger.info(f"Parsing Pfam hits from {self.retain_pfam_hits}")
            pfam_hits = self._parse_pfam_hits()
        
        # Parse CDS scores
        logger.info(f"Parsing CDS scores from {cds_scores_file}")
        cds_scores = self._parse_cds_scores(cds_scores_file)
        
        # Calculate min length for auto-accept based on GC content
        min_length_auto_accept = self._min_length_auto_accept()

        # Frame contrast required of short ORFs
        short_orf_nt = 3 * self.short_orf_length
        short_min_contrast = 0.0
        if self.short_orf_min_contrast > 0 and training_orfs_file and hexamer_scores_file:
            mu = self._coding_contrast(training_orfs_file, hexamer_scores_file)
            if mu:
                short_min_contrast = self.short_orf_min_contrast * mu
                logger.info("ORFs shorter than %d codons need frame contrast >= %.2f (%g bases x %.4f) "
                            "without homology support", self.short_orf_length, short_min_contrast,
                            self.short_orf_min_contrast, mu)
        
        # Parse GFF3 and select best ORFs
        logger.info(f"Selecting best ORFs from {self.gff3_file}")
        from pytransdecoder.core.gff3_parser import parse_gff3
        
        orfs_by_transcript = {}
        for orf in parse_gff3(self.gff3_file):
            # one candidate per ORF: its mRNA feature carries the ORF ID used for
            # scores and homology hits (gene/exon/CDS lines would otherwise enter as
            # unscored duplicates and could displace the ORF in overlap filtering)
            if orf['type'] != 'mRNA':
                continue
            transcript_id = orf['seqid']
            if transcript_id not in orfs_by_transcript:
                orfs_by_transcript[transcript_id] = []
            orfs_by_transcript[transcript_id].append(orf)
        
        # Select best ORFs per transcript
        selected_orfs = []
        for transcript_id in sorted(orfs_by_transcript.keys()):
            orfs = orfs_by_transcript[transcript_id]
            
            # Filter and score ORFs
            candidates = []
            for orf in orfs:
                orf_id = orf['attributes'].get('ID', [''])[0]
                orf_length = orf['end'] - orf['start'] + 1
                
                # Get homology count
                homology_count = 0
                blast_info = ""
                pfam_info = ""
                if orf_id in blast_hits:
                    homology_count += 1
                    blast_info = blast_hits[orf_id]
                if orf_id in pfam_hits:
                    homology_count += 1
                    pfam_info = pfam_hits[orf_id]
                
                # Get CDS scores (6 frames)
                scores = cds_scores.get(orf_id, [0.0] * 6)
                
                # Apply selection criteria
                frame_contrast = scores[0] - max(scores[1], scores[2])
                meets_criteria = (
                    homology_count > 0 or
                    orf_length >= min_length_auto_accept or
                    (scores[0] > 0 and frame_contrast > 0 and
                     (orf_length >= short_orf_nt or frame_contrast >= short_min_contrast))
                )
                
                if meets_criteria:
                    candidates.append({
                        'orf': orf,
                        'orf_id': orf_id,
                        'length': orf_length,
                        'homology_count': homology_count,
                        'scores': scores,
                        'blast_info': blast_info,
                        'pfam_info': pfam_info
                    })
            
            if not candidates:
                continue
            
            # Prioritize ORFs
            candidates.sort(key=lambda x: (
                -x['homology_count'],
                -x['scores'][0],
                -x['length'],
                # exact ties (e.g. identical repeated ORFs): leftmost, then plus strand
                min(x['orf']['start'], x['orf']['end']),
                max(x['orf']['start'], x['orf']['end']),
                x['orf']['strand']
            ))
            
            # Select ORFs
            if self.single_best_only:
                # the top candidate by homology, then coding score, then length (the priority
                # used for overlap removal); ranking by length alone favours long antisense ORFs
                selected_orfs.append(candidates[0])
            else:
                # Remove overlapping ORFs
                selected_orfs.extend(self._remove_overlapping_orfs(candidates))
        
        # Write output GFF3
        logger.info(f"Writing {len(selected_orfs)} selected ORFs to {output_file}")
        with open(output_file, 'w') as f:
            f.write("##gff-version 3\n")
            for candidate in selected_orfs:
                orf = candidate['orf']
                
                # Update attributes with score and homology info
                attrs = orf['attributes'].copy()
                attrs['score'] = [str(candidate['scores'][0])]
                if candidate['blast_info']:
                    attrs['blast'] = [candidate['blast_info']]
                if candidate['pfam_info']:
                    attrs['pfam'] = [candidate['pfam_info']]
                
                # Write GFF3 line
                attr_str = ';'.join([f"{k}={','.join(v)}" for k, v in attrs.items()])
                f.write(f"{orf['seqid']}\ttransdecoder\t{orf['type']}\t{orf['start']}\t{orf['end']}\t"
                       f"{orf['score']}\t{orf['strand']}\t{orf['phase']}\t{attr_str}\n")
        
        checkpoint.write_text("OK")
        return output_file
    
    def _parse_blastp_hits(self) -> Dict[str, str]:
        """Parse BLASTP hits file (outfmt 6 format)"""
        hits = {}
        with open(self.retain_blastp_hits) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 12:
                    query_id = parts[0]
                    subject_id = parts[1]
                    evalue = parts[10]
                    bitscore = parts[11]
                    if query_id not in hits:
                        hits[query_id] = f",blast:{subject_id}|{evalue}|{bitscore}"
        return hits
    
    def _parse_pfam_hits(self) -> Dict[str, str]:
        """Parse Pfam hits file (domtblout format from hmmscan/hmmsearch)"""
        hits = {}
        with open(self.retain_pfam_hits) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) >= 13:
                    # Try hmmscan format first (pfam_acc in column 1)
                    if parts[1].startswith('PF'):
                        orf_id = parts[3]
                        pfam_name = parts[0]
                        pfam_acc = parts[1]
                        domain_evalue = parts[12]
                    # Try hmmsearch format (pfam_acc in column 4)
                    elif len(parts) >= 13 and parts[4].startswith('PF'):
                        orf_id = parts[0]
                        pfam_name = parts[3]
                        pfam_acc = parts[4]
                        domain_evalue = parts[12]
                    else:
                        continue
                    
                    if orf_id not in hits:
                        hits[orf_id] = ""
                    hits[orf_id] += f",pfam:{pfam_name}|{pfam_acc}|{domain_evalue}"
        return hits
    
    def _parse_cds_scores(self, cds_scores_file: Path) -> Dict[str, List[float]]:
        """Parse CDS scores file"""
        scores = {}
        with open(cds_scores_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 9:
                    orf_id = parts[0]
                    # scores: markov_order, seq_length, score_1..6
                    frame_scores = [float(parts[i]) for i in range(3, 9)]
                    scores[orf_id] = frame_scores
        return scores
    
    def _min_length_auto_accept(self) -> float:
        """ORF length (nt) at or above which an ORF is retained regardless of its coding score."""
        if self.retain_long_orfs_mode != 'dynamic':
            logger.info(f"Strict mode: min auto-accept length = {self.retain_long_orfs_length} nt")
            return self.retain_long_orfs_length
        if self.stop_codons != STANDARD_STOP_CODONS:
            # the GC table was simulated with TAA/TAG/TGA; with fewer stops random ORFs are much
            # longer (e.g. TGA-only ciliate code), so no ORF is retained on length alone
            logger.info("Dynamic mode: genetic code stops %s differ from TAA/TAG/TGA; "
                        "no length-based auto-accept", ",".join(sorted(self.stop_codons)))
            return float('inf')
        min_length = self._get_min_length_by_gc()
        logger.info(f"Dynamic mode: min auto-accept length = {min_length} nt")
        return min_length

    def _get_min_length_by_gc(self) -> float:
        """
        Calculate minimum ORF length based on GC content (dynamic mode)
        
        Returns:
            Minimum length threshold in nucleotides
        """
        # Calculate GC content from transcripts
        gc_count = 0
        total_count = 0
        
        with open(self.transcripts_file) as f:
            for record in read_fasta(f):
                sequence = str(record.seq).upper()
                gc_count += sequence.count('G') + sequence.count('C')
                # as Perl's base_freqs.dat: A/C/G/T only
                total_count += sum(sequence.count(b) for b in 'ACGT')
        
        gc_pct = (gc_count / total_count * 100) if total_count > 0 else 50
        logger.info(f"Transcriptome GC content: {gc_pct:.1f}%")
        
        # Find appropriate threshold from lookup table
        min_length = 100  # default
        for gc_threshold, length_threshold in self.gc_to_min_length:
            if gc_pct <= gc_threshold:
                min_length = length_threshold
                break
        else:
            # GC > 80%: as Perl, rely on the coding score alone (no length-based auto-accept)
            min_length = float('inf')
        
        return min_length
    
    def _remove_overlapping_orfs(self, candidates: List[Dict]) -> List[Dict]:
        """
        Remove overlapping ORFs, keeping higher priority ones
        
        ORFs with homology support are always retained, even if they overlap.
        This matches Perl's behavior: $gene_entry->{homology_count} || ! &has_sufficient_overlap()
        
        Args:
            candidates: List of candidate ORFs (already sorted by priority)
            
        Returns:
            List of non-overlapping ORFs (plus all ORFs with homology)
        """
        MAX_PCT_OVERLAP = 10
        
        selected = []
        for candidate in candidates:
            # Always include ORFs with homology support (Perl behavior)
            if candidate['homology_count'] > 0:
                selected.append(candidate)
                continue
            
            orf = candidate['orf']
            start = orf['start']
            end = orf['end']
            strand = orf['strand']
            length = end - start + 1
            
            # Check for overlap with already selected ORFs
            # NOTE: Perl checks overlap regardless of strand to avoid selecting
            # multiple ORFs at the same genomic location (even on opposite strands)
            overlaps = False
            for selected_candidate in selected:
                sel_orf = selected_candidate['orf']
                # Perl does NOT check strand - it filters overlaps regardless of strand
                # This prevents selecting ORFs on opposite strands at same location
                
                sel_start = sel_orf['start']
                sel_end = sel_orf['end']
                sel_length = sel_end - sel_start + 1
                
                # Calculate overlap
                overlap_start = max(start, sel_start)
                overlap_end = min(end, sel_end)
                
                if overlap_start <= overlap_end:
                    overlap_length = overlap_end - overlap_start + 1
                    pct_overlap_curr = (overlap_length / length) * 100
                    pct_overlap_sel = (overlap_length / sel_length) * 100
                    
                    if pct_overlap_curr > MAX_PCT_OVERLAP or pct_overlap_sel > MAX_PCT_OVERLAP:
                        overlaps = True
                        break
            
            if not overlaps:
                selected.append(candidate)
        
        return selected
    
    def _coding_contrast(self, training_orfs_file: Path, hexamer_scores_file: Path) -> Optional[float]:
        """Median per-base frame contrast of the training ORFs (None without training ORFs)."""
        if self._coding_contrast_cache is None:
            scores = self._load_hexamer_scores(hexamer_scores_file)
            per_base = []
            with open(training_orfs_file) as f:
                for record in read_fasta(f):
                    seq = str(record.seq).upper()
                    if len(seq) >= 30:
                        per_base.append(self._frame_contrast(seq, scores) / len(seq))
            if not per_base:
                return None
            per_base.sort()
            mid = len(per_base) // 2
            self._coding_contrast_cache = (per_base[mid] if len(per_base) % 2
                                           else (per_base[mid - 1] + per_base[mid]) / 2)
        return self._coding_contrast_cache

    def _refine_start_codons(self, gff3_file: Path, training_orfs_file: Path,
                             hexamer_scores_file: Path) -> Path:
        """
        Decide, for each ORF open at its 5' end (no start codon), whether it is a true
        5'-partial ORF or extends into a 5' UTR and should start at its first in-frame ATG.

        Two signals are combined:
        - frame contrast of the segment between the open end and the ATG: the coding Markov
          score in frame minus the better of the two shifted frames, as in ORF selection.
          Coding sequence scores well above its shifted frames; UTR does not.
        - Kozak context of the ATG: a log-odds model of -6..-1 and +4..+5 trained on the
          training ORFs' starts against in-frame ATGs inside them.
        The ORF is trimmed to the ATG when the segment's contrast is below
        (start_contrast_alpha * length + start_kozak_weight * kozak) times the training ORFs'
        median per-base contrast: the Kozak evidence counts the same whatever the segment's length.

        Returns:
            Path to the revised GFF3 file
        """
        checkpoint = self.checkpoints_dir / "refine_starts.ok"
        revised_gff3 = Path(str(gff3_file) + ".revised_starts.gff3")

        if checkpoint.exists():
            logger.info("Start codon refinement already done (checkpoint exists)")
            return revised_gff3 if revised_gff3.exists() else gff3_file

        scores = self._load_hexamer_scores(hexamer_scores_file)
        coding_contrast = self._coding_contrast(training_orfs_file, hexamer_scores_file)
        if coding_contrast is None:
            logger.warning("No training ORFs for start refinement; skipping")
            checkpoint.write_text("OK")
            return gff3_file
        logger.info("Training ORFs: median frame contrast %.4f per base; open 5' ends kept when "
                    "the upstream segment reaches %.2f of it", coding_contrast, self.start_contrast_alpha)

        transcripts = self._load_transcripts()
        kozak = self._train_kozak_model(training_orfs_file, transcripts)

        num_revised = self._revise_start_sites(
            gff3_file, revised_gff3, transcripts, scores, coding_contrast, kozak
        )
        logger.info("Revised %s start codon positions", num_revised)

        checkpoint.write_text("OK")
        return revised_gff3

    def _load_hexamer_scores(self, hexamer_scores_file: Path) -> Dict[str, float]:
        scores = {}
        with open(hexamer_scores_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    scores[parts[0]] = float(parts[3])
        return scores

    def _frame_contrast(self, sequence: str, scores: Dict[str, float]) -> float:
        """In-frame Markov score minus the better of the two shifted-frame scores."""
        return self._score_cds_via_markov(sequence, 5, scores) - max(
            self._score_cds_via_markov(sequence[1:], 5, scores),
            self._score_cds_via_markov(sequence[2:], 5, scores),
        )

    def _train_kozak_model(self, training_orfs_file: Path, transcripts: Dict[str, str]):
        """
        Start-codon (Kozak) context model (core/kozak.py): log-odds of -6..-1 and +4..+5 for the
        training ORFs' starts vs in-frame ATGs inside them. Writes start_refinement.kozak.tsv,
        .kozak.counts.tsv and .kozak.logo.svg (or .kozak.SKIPPED) in the work directory.

        Returns the model, or None if fewer than kozak.MIN_TRAINING_STARTS training ORFs start
        with ATG.
        """
        from pytransdecoder.core import kozak
        from pytransdecoder.core.sequence import reverse_complement

        prefix = self.workdir / "start_refinement"
        pos_counts, neg_counts = kozak.new_counts(), kozak.new_counts()
        num_starts = 0
        pattern = re.compile(r"(\S+):(\d+)-(\d+)\(([+-])\)")
        with open(training_orfs_file) as handle:
            for record in read_fasta(handle):
                match = pattern.search(record.description)
                if not match or match.group(1) not in transcripts:
                    continue
                seq = transcripts[match.group(1)].upper()
                lo, hi = sorted((int(match.group(2)), int(match.group(3))))
                if match.group(4) == "-":
                    seq = reverse_complement(seq)
                    lo, hi = len(seq) - hi + 1, len(seq) - lo + 1
                start, end = lo - 1, hi
                if seq[start:start + 3] != "ATG":
                    continue
                num_starts += 1
                kozak.add_site(pos_counts, seq, start)
                for i in range(start + 3, end - 3, 3):
                    if seq[i:i + 3] == "ATG":
                        kozak.add_site(neg_counts, seq, i)

        skipped = Path(f"{prefix}.kozak.SKIPPED")
        skipped.unlink(missing_ok=True)
        if num_starts < kozak.MIN_TRAINING_STARTS:
            msg = (f"only {num_starts} training ORFs start with ATG (minimum {kozak.MIN_TRAINING_STARTS}); "
                   "start refinement uses frame contrast only")
            logger.warning(msg)
            skipped.write_text(msg + "\n")
            return None

        model = kozak.log_odds_model(pos_counts, neg_counts)
        kozak.write_model(model, num_starts, Path(f"{prefix}.kozak.tsv"))
        kozak.write_counts(pos_counts, neg_counts, Path(f"{prefix}.kozak.counts.tsv"))
        kozak.write_logo_svg(pos_counts, model, num_starts, Path(f"{prefix}.kozak.logo.svg"))
        logger.info("Kozak model from %d training starts: %s.kozak.tsv, %s.kozak.logo.svg",
                    num_starts, prefix, prefix)
        # score with the model as written (3 decimals), as Perl TransDecoder does
        return kozak.read_model(Path(f"{prefix}.kozak.tsv"))

    @staticmethod
    def _kozak_score(model, sequence: str, atg_index: int) -> float:
        from pytransdecoder.core import kozak
        return kozak.score(model, sequence, atg_index)

    def _revise_start_sites(
        self,
        gff3_file: Path,
        revised_gff3_file: Path,
        transcripts: Dict[str, str],
        scores: Dict[str, float],
        coding_contrast: float,
        kozak=None,
    ) -> int:
        """Trim open 5' ends that extend into a UTR to their first in-frame ATG
        (see _refine_start_codons)."""
        from pytransdecoder.core.sequence import reverse_complement

        num_revised = 0
        with open(gff3_file) as src, open(revised_gff3_file, "w") as dst:
            for line in src:
                parts = line.rstrip("\n").split("\t")
                if line.startswith("#") or len(parts) != 9 or parts[2] != "mRNA":
                    dst.write(line)
                    continue

                seqid, start, end, strand = parts[0], int(parts[3]), int(parts[4]), parts[6]
                transcript_seq = transcripts.get(seqid)
                if transcript_seq is None:
                    dst.write(line)
                    continue

                if strand == "-":
                    oriented = reverse_complement(transcript_seq)
                    orf_start = len(transcript_seq) - end
                else:
                    oriented = transcript_seq
                    orf_start = start - 1
                orf = oriented[orf_start:orf_start + end - start + 1]
                if orf[:3] == "ATG":
                    dst.write(line)
                    continue

                atg = next((i for i in range(3, len(orf) - 3, 3) if orf[i:i + 3] == "ATG"), None)
                if atg is None:
                    dst.write(line)
                    continue

                contrast = self._frame_contrast(orf[:atg], scores)
                kozak_score = self._kozak_score(kozak, oriented, orf_start + atg) if kozak else 0.0
                cutoff = self.start_contrast_alpha * atg + self.start_kozak_weight * kozak_score
                if contrast >= cutoff * coding_contrast:
                    dst.write(line)  # upstream segment looks coding: 5' partial
                    continue

                if strand == "-":
                    parts[4] = str(end - atg)
                else:
                    parts[3] = str(start + atg)
                attrs = self._parse_attrs(parts[8])
                attrs["start_revised"] = ["true"]
                attrs["upstream_frame_contrast"] = [f"{contrast / (coding_contrast * atg):.3f}"]
                attrs["kozak_score"] = [f"{kozak_score:.2f}"]
                parts[8] = self._format_attrs(attrs)
                dst.write("\t".join(parts) + "\n")
                num_revised += 1

        return num_revised

    def _load_transcripts(self) -> Dict[str, str]:
        transcripts = {}
        with open(self.transcripts_file) as handle:
            for record in read_fasta(handle):
                transcripts[record.id] = str(record.seq).upper()
        return transcripts

    def _parse_attrs(self, attrs_str: str) -> Dict[str, List[str]]:
        attrs = {}
        for entry in attrs_str.split(";"):
            if "=" not in entry:
                continue
            key, value = entry.split("=", 1)
            attrs[key] = value.split(",")
        return attrs

    def _format_attrs(self, attrs: Dict[str, List[str]]) -> str:
        return ";".join(f"{key}={','.join(values)}" for key, values in attrs.items())

    def _translate_brief(self, sequence: str) -> str:
        from pytransdecoder.core.translator import Translator

        translator = Translator(self.genetic_code)
        return translator.translate(sequence, frame=0)
    
    def _generate_final_outputs(self, gff3_file: Path):
        """
        Generate final output files (GFF3, BED, PEP, CDS)
        
        Creates four output files in the same directory as the input transcripts:
        - transcripts.transdecoder.gff3 - GFF3 format predictions
        - transcripts.transdecoder.bed - BED format predictions  
        - transcripts.transdecoder.pep - Predicted peptide sequences
        - transcripts.transdecoder.cds - Predicted CDS sequences
        
        Args:
            gff3_file: Input GFF3 file with selected/refined ORFs
        """
        checkpoint = self.checkpoints_dir / "final_outputs.ok"
        
        if checkpoint.exists():
            logger.info("Final outputs already generated (checkpoint exists)")
            return
        
        # Output file paths
        # When output_dir is specified, write final outputs there
        # Otherwise write to the directory containing the transcripts file
        base_name = self.transcripts_file.name
        output_base = self.output_dir / f"{base_name}.transdecoder"
        
        final_gff3 = Path(str(output_base) + ".gff3")
        final_bed = Path(str(output_base) + ".bed")
        final_pep = Path(str(output_base) + ".pep")
        final_cds = Path(str(output_base) + ".cds")
        
        # Write final GFF3 with full gene/mRNA/CDS hierarchy
        logger.info(f"Creating {final_gff3}")
        self._write_final_gff3(gff3_file, final_gff3)
        
        # Generate BED file
        logger.info(f"Creating {final_bed}")
        self._gff3_to_bed(gff3_file, final_bed)
        
        # Generate peptide file
        logger.info(f"Creating {final_pep}")
        self._gff3_to_proteins(gff3_file, final_pep, seq_type='pep')
        
        # Generate CDS file
        logger.info(f"Creating {final_cds}")
        self._gff3_to_proteins(gff3_file, final_cds, seq_type='cds')
        
        logger.info(f"Final outputs written to: {base_name}.transdecoder.*")
        checkpoint.write_text("OK")
    
    def _final_models(self, gff3_file: Path):
        """
        Transcript-level coding models for the selected ORFs: one exon spanning the transcript,
        the CDS inside it, oriented by the ORF strand (as Perl TransDecoder reports them).

        Yields (gene_id, CodingModel).
        """
        from pytransdecoder.core.gene_model import CodingModel
        from pytransdecoder.core.gff3_parser import parse_gff3
        from pytransdecoder.core.sequence import reverse_complement
        from pytransdecoder.core.translator import Translator

        translator = Translator(self.genetic_code)
        transcripts = self._load_transcripts()
        for feature in parse_gff3(gff3_file):
            if feature['type'] != 'mRNA' or feature['seqid'] not in transcripts:
                continue
            attrs = feature['attributes']
            seq = transcripts[feature['seqid']]
            start, end, strand = feature['start'], min(feature['end'], len(seq)), feature['strand']
            orf_seq = seq[start - 1:end]
            if strand == '-':
                orf_seq = reverse_complement(orf_seq)
            orf_type = self._get_orf_type_from_protein(translator.translate(orf_seq, frame=0))
            extra = ";".join(f"{k}={','.join(v)}" for k, v in attrs.items()
                             if k not in ("ID", "Parent", "Name", "blast", "pfam"))
            model = CodingModel(seqid=feature['seqid'], strand=strand, exons=[(1, len(seq))],
                                cds_lend=start, cds_rend=end, model_id=attrs.get('ID', [''])[0],
                                name=self._orf_com_name(attrs, strand, orf_type), extra_attrs=extra)
            yield attrs.get('Parent', [''])[0], model

    def _write_final_gff3(self, src_gff3: Path, out_gff3: Path):
        """
        Final GFF3 as in Perl TransDecoder: per ORF a gene and mRNA spanning the transcript,
        the transcript as one exon, the CDS, and five_prime_UTR / three_prime_UTR features.
        """
        from pytransdecoder.core.gene_model import CodingGene

        with open(out_gff3, 'w') as out:
            for gene_id, model in self._final_models(src_gff3):
                out.write(CodingGene(gene_id=gene_id, models=[model]).gff3_text())

    def _gff3_to_bed(self, gff3_file: Path, bed_file: Path):
        """BED12 of the final ORFs as in Perl TransDecoder (thickStart/End = CDS)."""
        with open(bed_file, 'w') as out:
            for gene_id, model in self._final_models(gff3_file):
                name = f"ID={model.model_id};{gene_id};{model.name.replace(' ', '_')}"
                out.write(model.bed12(name) + "\n")

    def _orf_com_name(self, attrs: Dict[str, List[str]], strand: str, orf_type: str) -> str:
        """Perl-style ORF description: 'ORF type:<type> (<strand>),score=<s>[,blast...][,pfam...]'."""
        com_name = f"ORF type:{orf_type} ({strand})"
        score = attrs.get('score', [''])[0]
        if score:
            try:
                score = f"{float(score):.2f}"
            except ValueError:
                pass
            com_name += f",score={score}"
        for key in ('blast', 'pfam'):
            info = attrs.get(key, [''])[0]
            if info:
                com_name += info if info.startswith(',') else f",{info}"
        return com_name

    @staticmethod
    def _get_orf_type_from_protein(protein: str) -> str:
        """Infer ORF completeness from the translated protein sequence."""
        has_start = protein.startswith("M")
        has_stop = protein.endswith("*")

        if has_start and has_stop:
            return "complete"
        if has_start:
            return "3prime_partial"
        if has_stop:
            return "5prime_partial"
        return "internal"

    def _format_transdecoder_fasta_header(
        self,
        feature: Dict[str, object],
        gene_id: str,
        sequence: str,
        protein: str,
        seq_type: str,
    ) -> str:
        """Reconstruct the Perl-style FASTA header for final outputs."""
        attrs = feature['attributes']
        model_id = attrs.get('ID', [''])[0]
        transcript_id = feature['seqid']
        start = feature['start']
        end = feature['end']
        strand = feature['strand']

        orf_type = self._get_orf_type_from_protein(protein)
        com_name = self._orf_com_name(attrs, strand, orf_type)

        seq_len = len(protein.rstrip('*')) if seq_type == 'pep' else len(sequence)

        header_parts = [model_id]
        if gene_id:
            header_parts.append(gene_id)
        header_parts.extend([
            com_name,
            f"len:{seq_len}",
            f"{transcript_id}:{start}-{end}({strand})",
        ])

        return " ".join(part for part in header_parts if part)
    
    def _gff3_to_proteins(self, gff3_file: Path, output_file: Path, seq_type: str = 'pep'):
        """
        Extract protein or CDS sequences from GFF3
        
        Args:
            gff3_file: Input GFF3 file
            output_file: Output FASTA file
            seq_type: 'pep' for protein sequences, 'cds' for CDS sequences
        """
        from pytransdecoder.core.gff3_parser import parse_gff3
        from pytransdecoder.core.translator import Translator
        
        # Create translator
        translator = Translator(self.genetic_code)
        
        # Load transcript sequences
        transcripts = {}
        with open(self.transcripts_file) as f:
            for record in read_fasta(f):
                transcripts[record.id] = str(record.seq).upper()
        
        # Extract sequences
        with open(output_file, 'w') as out:
            for feature in parse_gff3(gff3_file):
                transcript_id = feature['seqid']
                start = feature['start'] - 1  # Convert to 0-based
                end = feature['end']
                strand = feature['strand']
                orf_id = feature['attributes'].get('ID', [''])[0]
                gene_id = feature['attributes'].get('Parent', [''])[0]
                
                if transcript_id not in transcripts:
                    logger.warning(f"Transcript {transcript_id} not found in input file")
                    continue
                
                transcript_seq = transcripts[transcript_id]
                transcript_len = len(transcript_seq)
                
                # Handle invalid coordinates (negative or beyond transcript boundaries)
                # This can happen with partial ORFs at transcript boundaries
                if start < 0:
                    logger.warning(f"ORF {orf_id} has negative start coordinate ({feature['start']}), clamping to 0")
                    start = 0
                if end > transcript_len:
                    logger.warning(f"ORF {orf_id} end coordinate ({end}) exceeds transcript length ({transcript_len}), clamping")
                    end = transcript_len
                if start >= end:
                    logger.warning(f"ORF {orf_id} has invalid coordinates (start={start}, end={end}), skipping")
                    continue
                
                # Extract ORF sequence
                orf_seq = transcript_seq[start:end]
                
                # Reverse complement if on minus strand
                if strand == '-':
                    from pytransdecoder.core.sequence import reverse_complement
                    orf_seq = reverse_complement(orf_seq)

                protein_seq = translator.translate(orf_seq, frame=0)
                header = self._format_transdecoder_fasta_header(
                    feature=feature,
                    gene_id=gene_id,
                    sequence=orf_seq,
                    protein=protein_seq,
                    seq_type=seq_type,
                )
                
                # Generate output sequence
                if seq_type == 'pep':
                    output_seq = protein_seq
                else:
                    # CDS sequence
                    output_seq = orf_seq
                
                # Write FASTA record
                out.write(f">{header}\n")
                
                # Write sequence (60 chars per line)
                for i in range(0, len(output_seq), 60):
                    out.write(output_seq[i:i+60] + '\n')
