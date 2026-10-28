# Changelog

## v0.2.0 (unreleased)

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
