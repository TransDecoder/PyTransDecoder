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
import random
import re
from pathlib import Path
from typing import Optional, Dict, List, Tuple
from Bio import SeqIO

logger = logging.getLogger(__name__)


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
        single_best_only: bool = False,
        no_refine_starts: bool = False,
        genetic_code: str = "Standard"
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
            single_best_only: Retain only single best ORF per transcript
            no_refine_starts: Skip start codon refinement
            genetic_code: Genetic code to use
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
        
        # Derived paths
        self.workdir = self.output_dir / f"{self.transcripts_file.name}.transdecoder_dir"
        self.checkpoints_dir = self.workdir / "__checkpoints_TDpredict"
        
        # Input files from LongOrfs phase
        self.longest_orfs_prefix = self.workdir / "longest_orfs"
        self.cds_file = self.longest_orfs_prefix.with_suffix(self.longest_orfs_prefix.suffix + ".cds")
        self.gff3_file = self.longest_orfs_prefix.with_suffix(self.longest_orfs_prefix.suffix + ".gff3")
        self.pep_file = self.longest_orfs_prefix.with_suffix(self.longest_orfs_prefix.suffix + ".pep")
        self.base_freqs_file = self.workdir / "base_freqs.dat"
        
        # GC content to minimum ORF length mapping (from Perl code)
        # Based on 0.999 quantile of random ORFs
        self.gc_to_min_length = [
            (25, 465), (30, 510), (35, 555), (40, 590),
            (45, 645), (50, 749), (55, 797), (60, 927),
            (65, 1086), (70, 1358), (75, 1743), (80, 2422)
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
        best_orfs_gff3 = self._select_best_orfs(cds_scores_file)
        
        # Step 5: Refine start codons (if enabled)
        if not self.no_refine_starts:
            logger.info("Refining start codon predictions...")
            best_orfs_gff3 = self._refine_start_codons(best_orfs_gff3, top_cds_file)
        
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
            for record in SeqIO.parse(f, "fasta"):
                seq_len = len(record.seq)
                
                # Filter by CDS nucleotide length (not protein length)
                if seq_len <= max_cds_length:
                    longest_orfs.append((record, seq_len))
        
        # Sort by length (longest first)
        longest_orfs.sort(key=lambda x: x[1], reverse=True)
        
        # Take top red_num
        longest_orfs = longest_orfs[:red_num]
        logger.info(f"Selected {len(longest_orfs)} ORFs for redundancy filtering")
        
        # Step 2: Remove redundant sequences
        logger.info("Removing redundant sequences...")
        nr_orfs = self._exclude_similar_proteins(longest_orfs)
        logger.info(f"After redundancy filtering: {len(nr_orfs)} ORFs retained "
                   f"({len(nr_orfs)/len(longest_orfs)*100:.1f}%)")
        
        # Step 3: Select top N from non-redundant set
        nr_orfs_sorted = sorted(nr_orfs, key=lambda x: len(x.seq), reverse=True)
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
        
        # Process in order of longest first (already sorted)
        for record, _ in orf_records:
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
            for record in SeqIO.parse(f, "fasta"):
                sequence = str(record.seq).upper()
                seq_len = len(sequence)
                
                # Count k-mers for Markov orders 0-5 (producing 1-mers to 6-mers)
                for markov_order in range(6):
                    for i in range(markov_order, seq_len):
                        frame = i % 3
                        
                        # Skip stop codons at the end
                        if (i == seq_len - 3 and frame == 0):
                            codon = sequence[i:i+3]
                            if codon in ('TAA', 'TAG', 'TGA'):
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
        background_probs = {}
        if not self.base_freqs_file.exists():
            # If base_freqs.dat doesn't exist, calculate it now
            logger.warning(f"Base frequencies file not found: {self.base_freqs_file}")
            logger.info("Calculating base frequencies from transcripts...")
            background_probs = self._compute_base_frequencies()
        else:
            with open(self.base_freqs_file) as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        base, count, ratio = parts[0], parts[1], parts[2]
                        background_probs[base] = float(ratio)
        
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
            background_prob = background_probs.get(last_base, 0.25)
            
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
    
    def _compute_base_frequencies(self) -> Dict[str, float]:
        """
        Compute background base frequencies from transcripts
        
        Returns:
            Dictionary mapping base (A, C, G, T) to frequency
        """
        base_counter = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        
        with open(self.transcripts_file) as f:
            for record in SeqIO.parse(f, "fasta"):
                sequence = str(record.seq).upper()
                for base in sequence:
                    if base in base_counter:
                        base_counter[base] += 1
        
        total = sum(base_counter.values())
        background_probs = {base: count / total for base, count in base_counter.items()}
        
        # Write to file for future use
        with open(self.base_freqs_file, 'w') as f:
            for base in sorted(background_probs.keys()):
                count = base_counter[base]
                ratio = background_probs[base]
                f.write(f"{base}\t{count}\t{ratio:.3f}\n")
        
        logger.info(f"Base frequencies: {background_probs}")
        return background_probs
    
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
            for record in SeqIO.parse(f, "fasta"):
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
                if codon in ('TAA', 'TAG', 'TGA'):
                    break
            
            # Look up score
            framed_kmer = f"{kmer}-{frame}"
            loglikelihood = scores.get(framed_kmer, 0.0)
            
            total_score += loglikelihood
        
        return total_score
    
    def _select_best_orfs(self, cds_scores_file: Path) -> Path:
        """
        Select best ORFs per transcript based on scores and homology
        
        Selection criteria (any of):
        1. Has BLAST hit
        2. Has Pfam hit
        3. Score[frame0] > 0 AND score[frame0] > max(score[frame1], score[frame2])
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
        if self.retain_long_orfs_mode == 'dynamic':
            min_length_auto_accept = self._get_min_length_by_gc()
            logger.info(f"Dynamic mode: min auto-accept length = {min_length_auto_accept} nt")
        else:
            min_length_auto_accept = self.retain_long_orfs_length
            logger.info(f"Strict mode: min auto-accept length = {min_length_auto_accept} nt")
        
        # Parse GFF3 and select best ORFs
        logger.info(f"Selecting best ORFs from {self.gff3_file}")
        from pytransdecoder.core.gff3_parser import parse_gff3
        
        orfs_by_transcript = {}
        for orf in parse_gff3(self.gff3_file):
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
                meets_criteria = (
                    homology_count > 0 or
                    orf_length >= min_length_auto_accept or
                    (scores[0] > 0 and scores[0] > max(scores[1], scores[2]))
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
                -x['length']
            ))
            
            # Select ORFs
            if self.single_best_only:
                # Re-sort by homology then length
                candidates.sort(key=lambda x: (
                    -x['homology_count'],
                    -x['length']
                ))
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
    
    def _get_min_length_by_gc(self) -> int:
        """
        Calculate minimum ORF length based on GC content (dynamic mode)
        
        Returns:
            Minimum length threshold in nucleotides
        """
        # Calculate GC content from transcripts
        gc_count = 0
        total_count = 0
        
        with open(self.transcripts_file) as f:
            for record in SeqIO.parse(f, "fasta"):
                sequence = str(record.seq).upper()
                gc_count += sequence.count('G') + sequence.count('C')
                total_count += len(sequence)
        
        gc_pct = (gc_count / total_count * 100) if total_count > 0 else 50
        logger.info(f"Transcriptome GC content: {gc_pct:.1f}%")
        
        # Find appropriate threshold from lookup table
        min_length = 100  # default
        for gc_threshold, length_threshold in self.gc_to_min_length:
            if gc_pct <= gc_threshold:
                min_length = length_threshold
                break
        else:
            # GC > 80%, use last entry
            min_length = self.gc_to_min_length[-1][1]
        
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
    
    def _refine_start_codons(self, gff3_file: Path, training_orfs_file: Path) -> Path:
        """
        Refine start codon predictions using PWM
        
        This is a complex feature that trains a Position Weight Matrix (PWM)
        on the sequence context around start codons and uses it to refine
        predictions for 5' partial ORFs.
        
        For now, this is a placeholder that returns the input file unchanged.
        Full PWM-based refinement can be implemented later if needed.
        
        Args:
            gff3_file: GFF3 file with selected ORFs
            training_orfs_file: Training ORFs file
            
        Returns:
            Path to refined GFF3 file (currently just returns input)
        """
        checkpoint = self.checkpoints_dir / "refine_starts.ok"
        revised_gff3 = Path(str(gff3_file) + ".revised_starts.gff3")
        
        if checkpoint.exists():
            logger.info("Start codon refinement already done (checkpoint exists)")
            return revised_gff3 if revised_gff3.exists() else gff3_file

        training_outputs = self._train_start_pwm(training_orfs_file)
        if training_outputs is None:
            logger.info("Insufficient training data for start codon refinement - skipping")
            checkpoint.write_text("OK")
            return gff3_file

        transcripts = self._load_transcripts()
        pwm_plus, pwm_minus, pwm_range, min_threshold = training_outputs
        num_revised = self._revise_start_sites(
            gff3_file, revised_gff3, transcripts, pwm_plus, pwm_minus, pwm_range, min_threshold
        )
        logger.info("Revised %s start codon positions", num_revised)
        
        checkpoint.write_text("OK")
        return revised_gff3

    def _train_start_pwm(self, training_orfs_file: Path):
        """
        Train PWM assets used by start codon refinement.

        Returns:
            Tuple of (pwm_plus, pwm_minus, pwm_range, min_threshold) or None if the
            training set does not contain enough valid start-site features.
        """
        from pytransdecoder.core.pwm import PWM, build_pwm, trapezoid_auc

        pwm_left = 20
        pwm_right = 10
        pwm_length = pwm_left + 3 + pwm_right
        out_prefix = self.workdir / "start_refinement"

        transcripts = self._load_transcripts()
        positive_features, negative_features = self._extract_pwm_features(
            training_orfs_file, transcripts, pwm_left, pwm_right
        )

        if not positive_features or not negative_features:
            return None

        self._write_feature_file(out_prefix.with_name(out_prefix.name + ".+.features"), positive_features)
        self._write_feature_file(out_prefix.with_name(out_prefix.name + ".-.features"), negative_features)

        pwm_plus = build_pwm(positive_features)
        pwm_minus = build_pwm(negative_features)
        pwm_plus_file = out_prefix.with_name(out_prefix.name + ".+.pwm")
        pwm_minus_file = out_prefix.with_name(out_prefix.name + ".-.pwm")
        pwm_plus.write(pwm_plus_file)
        pwm_minus.write(pwm_minus_file)

        enhanced_features, enhanced_pwm = self._deplete_feature_noise(positive_features, pwm_minus)
        enhanced_features_file = out_prefix.with_name(out_prefix.name + ".enhanced.+.features")
        enhanced_pwm_file = out_prefix.with_name(out_prefix.name + ".enhanced.+.pwm")
        self._write_feature_file(enhanced_features_file, enhanced_features)
        enhanced_pwm.write(enhanced_pwm_file)

        features_minus_file = out_prefix.with_name(out_prefix.name + ".-.features")
        feature_scores_file = out_prefix.with_name(out_prefix.name + ".feature.scores")
        feature_roc_file = out_prefix.with_name(out_prefix.name + ".feature.scores.roc")
        feature_auc_file = out_prefix.with_name(out_prefix.name + ".feature.scores.roc.auc")
        enhanced_scores_file = out_prefix.with_name(out_prefix.name + ".enhanced.feature.scores")
        enhanced_roc_file = out_prefix.with_name(out_prefix.name + ".enhanced.feature.scores.roc")
        enhanced_auc_file = out_prefix.with_name(out_prefix.name + ".enhanced.feature.scores.roc.auc")

        self._score_feature_sets(
            enhanced_features, negative_features, pwm_length, pwm_left, feature_scores_file=enhanced_scores_file
        )
        self._feature_scores_to_roc(enhanced_scores_file, enhanced_roc_file)
        self._compute_auc(enhanced_roc_file, enhanced_auc_file)

        # Preserve the non-enhanced files too for compatibility/debugging.
        self._score_feature_sets(
            positive_features, negative_features, pwm_length, pwm_left, feature_scores_file=feature_scores_file
        )
        self._feature_scores_to_roc(feature_scores_file, feature_roc_file)
        self._compute_auc(feature_roc_file, feature_auc_file)

        best_range, min_threshold = self._select_pwm_range_and_threshold(
            enhanced_auc_file, enhanced_roc_file
        )
        return enhanced_pwm, pwm_minus, best_range, min_threshold

    def _extract_pwm_features(
        self,
        training_orfs_file: Path,
        transcripts: Dict[str, str],
        pwm_left: int,
        pwm_right: int,
    ) -> Tuple[List[str], List[str]]:
        from pytransdecoder.core.sequence import reverse_complement

        positive_features = []
        negative_features = []
        pwm_length = pwm_left + 3 + pwm_right
        pattern = re.compile(r"(\S+):(\d+)-(\d+)\(([+-])\)")

        with open(training_orfs_file) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                match = pattern.search(record.description)
                if not match:
                    continue

                transcript_id = match.group(1)
                start_coord = int(match.group(2))
                end_coord = int(match.group(3))
                strand = match.group(4)

                transcript_seq = transcripts.get(transcript_id)
                if transcript_seq is None:
                    continue

                if strand == "-":
                    transcript_seq = reverse_complement(transcript_seq)
                    start_coord = len(transcript_seq) - max(start_coord, end_coord) + 1
                else:
                    start_coord = min(start_coord, end_coord)

                start_index = start_coord - 1
                if transcript_seq[start_index:start_index + 3] != "ATG":
                    continue

                feature_seq = self._extract_feature_seq(
                    transcript_seq, start_index, pwm_left, pwm_right, pwm_length
                )
                if feature_seq:
                    positive_features.append(feature_seq)

                downstream_seq = transcript_seq[start_coord + 1 :]
                negative_features.extend(
                    self._extract_all_start_features(
                        downstream_seq, pwm_left, pwm_right, pwm_length
                    )
                )

        return positive_features, negative_features

    def _extract_feature_seq(
        self,
        sequence: str,
        start_index: int,
        pwm_left: int,
        pwm_right: int,
        pwm_length: int,
    ) -> Optional[str]:
        begin = start_index - pwm_left
        end = begin + pwm_length
        if begin < 0 or end > len(sequence):
            return None
        feature_seq = sequence[begin:end].upper()
        if any(base not in "GATC" for base in feature_seq):
            return None
        return feature_seq

    def _extract_all_start_features(
        self,
        sequence: str,
        pwm_left: int,
        pwm_right: int,
        pwm_length: int,
    ) -> List[str]:
        features = []
        start = 0
        while True:
            pos = sequence.find("ATG", start)
            if pos == -1:
                break
            feature_seq = self._extract_feature_seq(sequence, pos, pwm_left, pwm_right, pwm_length)
            if feature_seq:
                features.append(feature_seq)
            start = pos + 1
        return features

    def _write_feature_file(self, filename: Path, features: List[str]) -> None:
        with open(filename, "w") as handle:
            for feature in features:
                handle.write(f"{feature}\n")

    def _deplete_feature_noise(self, positive_features: List[str], pwm_minus):
        from pytransdecoder.core.pwm import build_pwm

        rng = random.Random(1)
        features = positive_features[:]
        rng.shuffle(features)

        num_incorporate = max(1, int(len(features) * 30 / 100))
        init_features = features[:num_incorporate]
        remaining_features = features[num_incorporate:]

        pwm_plus = build_pwm(init_features)
        scored_features = []
        for feature in init_features:
            score = pwm_plus.score_plus_minus(feature, pwm_minus)
            if score is not None:
                scored_features.append({"score": score, "seq": feature})
        scored_features.sort(key=lambda item: item["score"])

        for feature in remaining_features:
            if not scored_features:
                break
            score = pwm_plus.score_plus_minus(feature, pwm_minus)
            if score is None or score <= scored_features[0]["score"]:
                continue

            purge_feature = scored_features.pop(0)
            pwm_plus.remove_feature(purge_feature["seq"])
            pwm_plus.add_feature(feature)
            pwm_plus.build()
            scored_features.append({"score": score, "seq": feature})

            for scored in scored_features:
                rescored = pwm_plus.score_plus_minus(scored["seq"], pwm_minus)
                scored["score"] = float("-inf") if rescored is None else rescored
            scored_features.sort(key=lambda item: item["score"])

        retained_features = []
        for scored_feature in scored_features:
            if scored_feature["score"] <= 0:
                pwm_plus.remove_feature(scored_feature["seq"])
            else:
                retained_features.append(scored_feature["seq"])

        pwm_plus.build()
        return retained_features, pwm_plus

    def _score_feature_sets(
        self,
        positive_features: List[str],
        negative_features: List[str],
        pwm_length: int,
        atg_position: int,
        feature_scores_file: Path,
    ) -> None:
        from pytransdecoder.core.pwm import build_pwm

        rng = random.Random(1)
        pwm_upstream_max = atg_position
        pwm_downstream_max = pwm_length - (atg_position + 3)
        up_down_combos = [
            (up, down)
            for up in range(1, pwm_upstream_max + 1)
            for down in range(1, pwm_downstream_max + 1)
        ]

        num_rounds = 5
        fraction_train = 0.75
        max_feature_select = 1000

        with open(feature_scores_file, "w") as handle:
            for _ in range(num_rounds):
                plus_train, plus_test = self._sample_features(positive_features, fraction_train, rng)
                minus_train, minus_test = self._sample_features(negative_features, fraction_train, rng)
                if not plus_train or not minus_train or not plus_test or not minus_test:
                    continue

                pwm_plus = build_pwm(plus_train)
                pwm_minus = build_pwm(minus_train)

                self._score_features(
                    handle, plus_test, pwm_plus, pwm_minus, up_down_combos, atg_position, "pos", max_feature_select
                )
                self._score_features(
                    handle, minus_test, pwm_plus, pwm_minus, up_down_combos, atg_position, "neg", max_feature_select
                )

    def _sample_features(
        self, features: List[str], fraction_train: float, rng: random.Random
    ) -> Tuple[List[str], List[str]]:
        features = features[:]
        rng.shuffle(features)
        num_train = int(fraction_train * len(features))
        num_train = min(max(num_train, 1), len(features) - 1) if len(features) > 1 else len(features)
        return features[:num_train], features[num_train:]

    def _score_features(
        self,
        handle,
        features: List[str],
        pwm_plus,
        pwm_minus,
        up_down_combos: List[Tuple[int, int]],
        atg_position: int,
        feature_set_type: str,
        max_feature_select: int,
    ) -> None:
        for up, down in up_down_combos:
            range_left = atg_position - up
            range_right = atg_position + 2 + down
            local_pwm_len = up + down
            for feature in features[:max_feature_select]:
                score = pwm_plus.score_plus_minus(feature, pwm_minus, pwm_range=(range_left, range_right))
                if score is None:
                    score_str = "NA"
                else:
                    score_str = f"{score / local_pwm_len:.3f}"
                handle.write(f"{up},{down}\t{feature_set_type}\t{score_str}\n")

    def _feature_scores_to_roc(self, feature_scores_file: Path, roc_file: Path) -> None:
        scores_by_category = {}
        min_score = None
        max_score = None

        with open(feature_scores_file) as handle:
            for line in handle:
                category, pos_or_neg, score_str = line.rstrip().split("\t")
                if score_str == "NA":
                    continue
                score = float(score_str)
                scores_by_category.setdefault(category, []).append((pos_or_neg, score))
                min_score = score if min_score is None else min(min_score, score)
                max_score = score if max_score is None else max(max_score, score)

        if min_score is None or max_score is None:
            roc_file.write_text("cat\tthresh\tTP\tTN\tFP\tFN\tTPR\tFPR\tF1\n")
            return

        delta = (max_score - min_score) / 10 if max_score != min_score else 1.0
        with open(roc_file, "w") as handle:
            handle.write("cat\tthresh\tTP\tTN\tFP\tFN\tTPR\tFPR\tF1\n")
            for category, score_entries in scores_by_category.items():
                threshold = min_score
                while threshold < max_score:
                    tp = tn = fp = fn = 0
                    for pos_or_neg, score in score_entries:
                        if pos_or_neg == "pos":
                            if score >= threshold:
                                tp += 1
                            else:
                                fn += 1
                        else:
                            if score >= threshold:
                                fp += 1
                            else:
                                tn += 1
                    tpr = tp / (tp + fn) if (tp + fn) else 0.0
                    fpr = fp / (fp + tn) if (fp + tn) else 0.0
                    denom = (2 * tp + fp + fn)
                    f1 = (2 * tp) / denom if denom else 0.0
                    handle.write(
                        f"{category}\t{threshold}\t{tp}\t{tn}\t{fp}\t{fn}\t{tpr}\t{fpr}\t{f1}\n"
                    )
                    threshold += delta

    def _compute_auc(self, roc_file: Path, auc_file: Path) -> None:
        from pytransdecoder.core.pwm import trapezoid_auc

        points_by_category = {}
        with open(roc_file) as handle:
            next(handle, None)
            for line in handle:
                parts = line.rstrip().split("\t")
                if len(parts) < 8:
                    continue
                category = parts[0]
                tpr = float(parts[6])
                fpr = float(parts[7])
                points_by_category.setdefault(category, []).append((fpr, tpr))

        with open(auc_file, "w") as handle:
            for category in sorted(points_by_category):
                auc = trapezoid_auc(points_by_category[category])
                handle.write(f"{category}\t{auc}\n")

    def _select_pwm_range_and_threshold(self, auc_file: Path, roc_file: Path) -> Tuple[Tuple[int, int], float]:
        best_range = None
        best_auc = float("-inf")
        with open(auc_file) as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                range_str, auc_str = line.split("\t")
                auc = float(auc_str)
                if auc > best_auc:
                    best_auc = auc
                    best_range = range_str

        if best_range is None:
            raise ValueError("Unable to determine best PWM range")

        best_threshold = None
        best_f1 = float("-inf")
        with open(roc_file) as handle:
            next(handle, None)
            for line in handle:
                parts = line.rstrip().split("\t")
                if len(parts) < 9 or parts[0] != best_range:
                    continue
                threshold = float(parts[1])
                f1 = float(parts[8])
                if f1 > best_f1:
                    best_f1 = f1
                    best_threshold = threshold

        if best_threshold is None:
            raise ValueError("Unable to determine best PWM threshold")

        range_left, range_right = [int(value) for value in best_range.split(",")]
        return (range_left, range_right), best_threshold

    def _revise_start_sites(
        self,
        gff3_file: Path,
        revised_gff3_file: Path,
        transcripts: Dict[str, str],
        pwm_plus,
        pwm_minus,
        pwm_range: Tuple[int, int],
        min_threshold: float,
    ) -> int:
        from pytransdecoder.core.sequence import reverse_complement

        atg_pwm_pos = 20
        adj_dist = 30
        adj_pct = 15
        num_revised = 0
        alt_start_scores_file = self.workdir / "start_refinement.alt_start_scores"

        with open(gff3_file) as src, open(revised_gff3_file, "w") as dst, open(
            alt_start_scores_file, "w"
        ) as scores_handle:
            for line in src:
                stripped = line.rstrip("\n")
                if not stripped or stripped.startswith("#"):
                    dst.write(line)
                    continue

                parts = stripped.split("\t")
                if len(parts) != 9:
                    dst.write(line)
                    continue

                seqid, source, feature_type = parts[0], parts[1], parts[2]
                if feature_type != "mRNA":
                    dst.write(line)
                    continue

                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attrs = self._parse_attrs(parts[8])
                transcript_seq = transcripts.get(seqid)
                if transcript_seq is None:
                    dst.write(line)
                    continue

                oriented_seq = transcript_seq
                start_pos = start
                if strand == "-":
                    oriented_seq = reverse_complement(transcript_seq)
                    start_pos = len(oriented_seq) - end + 1

                start_index = start_pos - 1
                if oriented_seq[start_index:start_index + 3] == "ATG":
                    dst.write(line)
                    continue

                orf_len = end - start + 1
                max_search_pos = max(start_index + adj_dist, start_index + int(adj_pct * orf_len / 100))
                best_alt_start = None
                best_alt_score = None
                alt_score_entries = []

                search_pos = 0
                while True:
                    pos = oriented_seq.find("ATG", search_pos)
                    if pos == -1 or pos > max_search_pos:
                        break
                    if pos > start_index and (pos - start_index) % 3 == 0:
                        feature_seq_start = pos - atg_pwm_pos
                        feature_seq_end = feature_seq_start + pwm_plus.length
                        if feature_seq_start > 0 and feature_seq_end <= len(oriented_seq):
                            feature_seq = oriented_seq[feature_seq_start:feature_seq_end]
                            score = pwm_plus.score_plus_minus(
                                feature_seq, pwm_minus, pwm_range=(
                                    atg_pwm_pos - pwm_range[0] - 1,
                                    atg_pwm_pos + 2 + pwm_range[1] - 1,
                                )
                            )
                            if score is not None:
                                score = round(score, 3)
                                alt_score_entries.append(
                                    f"{pos}_{self._translate_brief(oriented_seq[pos:pos + 15])}_{score:.3f}"
                                )
                                if score >= min_threshold and (
                                    best_alt_score is None or score > best_alt_score
                                ):
                                    best_alt_start = pos
                                    best_alt_score = score
                    search_pos = pos + 1

                if alt_score_entries:
                    gene_id = attrs.get("Parent", [""])[0]
                    scores_handle.write(f"{seqid}\t{gene_id}\t" + "\t".join(alt_score_entries) + "\n")

                if best_alt_start is None:
                    dst.write(line)
                    continue

                best_alt_start += 1
                if strand == "-":
                    new_start = len(oriented_seq) - best_alt_start + 1
                    parts[4] = str(new_start)
                else:
                    new_start = best_alt_start
                    parts[3] = str(new_start)

                attrs["start_revised"] = ["true"]
                attrs["start_revised_score"] = [f"{best_alt_score:.3f}"]
                parts[8] = self._format_attrs(attrs)
                dst.write("\t".join(parts) + "\n")
                num_revised += 1

        return num_revised

    def _load_transcripts(self) -> Dict[str, str]:
        transcripts = {}
        with open(self.transcripts_file) as handle:
            for record in SeqIO.parse(handle, "fasta"):
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
    
    def _write_final_gff3(self, src_gff3: Path, out_gff3: Path):
        """
        Write final GFF3 with full gene/mRNA/CDS hierarchy.
        The best_candidates.gff3 only has mRNA entries; cdna_alignment_orf_to_genome_orf.pl
        requires CDS children to build gene models, so we expand each mRNA into:
            gene  (parent)
            mRNA  (child of gene)
            CDS   (child of mRNA, same coords)
        """
        with open(src_gff3) as fh, open(out_gff3, 'w') as out:
            out.write("##gff-version 3\n")
            for line in fh:
                if line.startswith('#'):
                    continue
                line = line.rstrip('\n')
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) != 9 or parts[2] != 'mRNA':
                    continue

                seqid, source, _, start, end, score, strand, phase, attrs = parts

                # Parse ID and Parent from attributes
                mrna_id, parent_id = '', ''
                for attr in attrs.split(';'):
                    if attr.startswith('ID='):
                        mrna_id = attr[3:]
                    elif attr.startswith('Parent='):
                        parent_id = attr[7:]

                # gene line
                out.write(f"{seqid}\t{source}\tgene\t{start}\t{end}\t.\t{strand}\t.\t"
                          f"ID={parent_id};Name=ORF_{parent_id}\n")
                # mRNA line (unchanged)
                out.write(line + '\n')
                # exon line (required by GFF3_utils2 for gene model building)
                out.write(f"{seqid}\t{source}\texon\t{start}\t{end}\t.\t{strand}\t.\t"
                          f"ID=exon.{mrna_id};Parent={mrna_id}\n")
                # CDS line
                out.write(f"{seqid}\t{source}\tCDS\t{start}\t{end}\t.\t{strand}\t0\t"
                          f"ID=cds.{mrna_id};Parent={mrna_id}\n")

    def _gff3_to_bed(self, gff3_file: Path, bed_file: Path):
        """Convert GFF3 to BED format"""
        from pytransdecoder.core.gff3_parser import parse_gff3
        
        with open(bed_file, 'w') as out:
            for feature in parse_gff3(gff3_file):
                # BED format: chrom, chromStart, chromEnd, name, score, strand
                # Note: BED uses 0-based start, 1-based end (half-open)
                chrom = feature['seqid']
                start = feature['start'] - 1  # Convert to 0-based
                end = feature['end']
                name = feature['attributes'].get('ID', [''])[0]
                score = feature['score']
                strand = feature['strand']
                
                # BED12 format with additional fields
                # For now, use simple BED6 format
                out.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")

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
        com_name = f"ORF type:{orf_type} ({strand})"

        score = attrs.get('score', [''])[0]
        if score:
            com_name += f",score={score}"

        blast_info = attrs.get('blast', [''])[0]
        if blast_info:
            com_name += blast_info if blast_info.startswith(',') else f",{blast_info}"

        pfam_info = attrs.get('pfam', [''])[0]
        if pfam_info:
            com_name += pfam_info if pfam_info.startswith(',') else f",{pfam_info}"

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
            for record in SeqIO.parse(f, "fasta"):
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
                    # Remove stop codon if present
                    if protein_seq.endswith('*'):
                        protein_seq = protein_seq[:-1]
                    output_seq = protein_seq
                else:
                    # CDS sequence
                    output_seq = orf_seq
                
                # Write FASTA record
                out.write(f">{header}\n")
                
                # Write sequence (60 chars per line)
                for i in range(0, len(output_seq), 60):
                    out.write(output_seq[i:i+60] + '\n')
