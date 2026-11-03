# Changelog

## v0.4.0 (2026-11-03)
Faster LongOrfs, with identical output: every output file is byte-identical to v0.3.0 (and v0.2.1), at any `--cpu`.
- **`longorfs --cpu N`** (and `pyTransdecoder --cpu N` now also applies to LongOrfs): ORF extraction in N
  worker processes, chunks of transcripts written back in input order; output files identical at any N.
  The input is read once (base counts are summed per chunk, exactly). With `--cpu 1` the input is streamed
  as before; with N > 1 it is held in memory. Inputs with repeated accessions fall back to one process
  (their `.pN` numbering continues across copies).
- LongOrfs with `--cpu 8`: median 3.4x faster than v0.3.0 on our benchmark inputs (a 429k-transcript Trinity
  assembly: 62 s -> 15 s).

## v0.3.0 (2026-11-02)
Faster, with identical predictions: every output file is byte-identical to v0.2.1, at any `--cpu`.
- **`--cpu N`** (`predict` and `pyTransdecoder`): candidate and decoy ORFs are scored in N worker processes.
  Work is split and reassembled in input order, and each decoy shuffle is seeded per transcript, so results
  do not depend on N.
- **Faster single-process steps**, each computing the same values in the same order:
  - NumPy Markov scoring (additions in the reference order; any non-standard score table uses the original loop);
  - codon-memoised translation built from Biopython's own results;
  - a one-pass codon scan and bisect stop lookup in the ORF finder;
  - faster k-mer counting, redundancy filtering, base counts, reverse complement and decoy shuffles;
  - transcripts and final models read once per run.
- Typical wall time (LongOrfs + Predict), on our benchmark inputs: 1.75-3.45x faster than v0.2.1 with
  `--cpu 1`, 5-10x faster with `--cpu 8` (e.g. a 429k-transcript GM12878 Trinity assembly: 30 min -> 3.7 min).
- New dependency: NumPy (>= 1.21).
- Biopython's "partial codon" warning is no longer printed during translation (output unchanged).

## v0.2.1 (2026-11-02)
- Packaging only, for the first release on PyPI (`pip install pytransdecoder`): author, summary and project URL;
  README links made absolute so they work on PyPI. No code changes; predictions are identical to v0.2.0.

## v0.2.0 (2026-11-02)

### New in PyTransDecoder only
- **Homology-assisted start refinement.** When Predict is given `--retain-blastp-hits` and/or
  `--retain-pfam-hits`, the hits decide whether an ORF that runs off the transcript's 5' end keeps
  its open end or is trimmed to its first in-frame ATG:
  - trim when the best BLASTP hit's subject N-terminus aligns within 5 aa of the ATG;
  - keep it 5' partial when the alignment starts >= 5 aa upstream of the ATG (with more of the
    homolog upstream), or when a Pfam domain envelope starts upstream of the ATG.

  The ab initio coding-vs-UTR rule decides where no rule fires or the calls conflict. On by default
  whenever hits are supplied; `--no-homology-start-assist` turns it off. Each ORF's evidence and
  decision are written to `start_refinement.homology.tsv`; ORFs trimmed by homology carry
  `start_homology=R2` in the GFF3. See `docs/start_codon_refinement.md`.
- Without homology hits, or with `--no-homology-start-assist`, the output is identical to the
  previous commit (3a793ff). TransDecoder (Perl) does not have this mode, so Perl and Python agree
  on ab initio runs but not on homology runs.

### Shared with TransDecoder (Perl) v7.0.0
- ORF acceptance by a target-decoy q-value from dinucleotide-shuffled transcripts:
  - `--max-qvalue`, default 0.05; `--precise` = 0.01, `--sensitive` = 0.2;
  - `--decoy-seed`, `--decoy-shuffles`.
- Presets for small or mostly partial inputs: `--preset`, `build-preset`, and bundled clade presets.
- Start refinement as two competing scores, coding vs UTR, with the Kozak score in bits:
  - `--start-contrast-alpha`, `--start-kozak-weight`;
  - `--no-refine-starts` turns refinement off.
- One ORF per transcript by default (`--retain-multiple-orfs` reports all).
- `genome-map` writes genome-coordinate ORFs (GFF3 and BED12).
- GC-dependent long-ORF retention from a regenerated simulation.
- Default minimum protein length of 50 aa.
- Deterministic tie-breaks.
