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
from pathlib import Path
from typing import Optional, Dict, Set, List, Tuple
import click
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
            overlaps = False
            for selected_candidate in selected:
                sel_orf = selected_candidate['orf']
                if sel_orf['strand'] != strand:
                    continue
                
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
        
        if checkpoint.exists():
            logger.info("Start codon refinement already done (checkpoint exists)")
            return gff3_file
        
        # TODO: Implement PWM-based start codon refinement
        # For now, just pass through the input file
        logger.info("Start codon refinement not yet fully implemented - skipping")
        
        checkpoint.write_text("OK")
        return gff3_file
    
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
        
        # Copy GFF3 to final location
        logger.info(f"Creating {final_gff3}")
        import shutil
        shutil.copy(gff3_file, final_gff3)
        
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
                
                # Generate output sequence
                if seq_type == 'pep':
                    # Translate to protein
                    protein_seq = translator.translate(orf_seq, frame=0)
                    # Remove stop codon if present
                    if protein_seq.endswith('*'):
                        protein_seq = protein_seq[:-1]
                    output_seq = protein_seq
                else:
                    # CDS sequence
                    output_seq = orf_seq
                
                # Write FASTA record
                # Include additional info from attributes
                gene_info = feature['attributes'].get('gene', [''])[0]
                type_info = feature['attributes'].get('type', [''])[0]
                
                header = f">{orf_id}"
                if gene_info:
                    header += f" gene={gene_info}"
                if type_info:
                    header += f" type={type_info}"
                
                out.write(f"{header}\n")
                
                # Write sequence (60 chars per line)
                for i in range(0, len(output_seq), 60):
                    out.write(output_seq[i:i+60] + '\n')


@click.command()
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
def main(
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
    TransDecoder.Predict - Final coding region predictions
    
    This command takes the output from TransDecoder.LongOrfs and identifies
    the most likely coding regions using:
    - Hexamer composition scoring
    - Homology support (BLAST and Pfam)
    - ORF length
    - Start codon refinement
    """
    # Setup logging
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # Run prediction
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


if __name__ == '__main__':
    main()
