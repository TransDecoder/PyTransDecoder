"""
Basic tests for PyTransDecoder core modules
"""

import subprocess
from pathlib import Path

import pytest
from pytransdecoder.core.translator import Translator
from pytransdecoder.core.sequence import reverse_complement
from pytransdecoder.core.orf_finder import ORFFinder
from pytransdecoder.core.models import ORF
from pytransdecoder.core.pwm import build_pwm, PWM
from pytransdecoder import cli
from pytransdecoder.predict import TransDecoderPredict


def test_translator_standard_code():
    """Test basic translation"""
    translator = Translator("Standard")
    
    # ATG = Met, GCA = Ala, TAA = stop
    sequence = "ATGGCATAA"
    protein = translator.translate(sequence)
    
    # Should be "MA*" (Met-Ala-stop)
    assert protein.startswith("MA")
    assert "*" in protein


def test_translator_stop_codons():
    """Test stop codon identification"""
    translator = Translator("Standard")
    stops = translator.get_stop_codons()
    
    assert "TAA" in stops
    assert "TAG" in stops
    assert "TGA" in stops
    assert len(stops) == 3


def test_translator_start_codons():
    """Test start codon identification"""
    translator = Translator("Standard")
    
    # Default: only ATG
    starts = translator.get_start_codons(allow_non_met=False)
    assert starts == ["ATG"]
    
    # With alternatives: includes ATG and others
    starts_alt = translator.get_start_codons(allow_non_met=True)
    assert "ATG" in starts_alt


def test_reverse_complement():
    """Test reverse complement"""
    seq = "ATGC"
    rc = reverse_complement(seq)
    assert rc == "GCAT"
    
    seq2 = "AAAA"
    rc2 = reverse_complement(seq2)
    assert rc2 == "TTTT"


def test_orf_finder_basic():
    """Test basic ORF finding"""
    finder = ORFFinder(min_protein_length=1)  # Allow short ORFs for testing
    
    # Simple ORF: ATG...TAA
    sequence = "ATGGCAGCAGCATAA"  # M A A A *
    orfs = finder.find_all_orfs(sequence, "test_seq")
    
    # Should find at least one ORF
    assert len(orfs) > 0
    
    # Check the longest ORF
    longest = orfs[0]
    assert longest.transcript_id == "test_seq"
    assert "M" in longest.protein
    assert longest.strand in ['+', '-']


def test_orf_finder_complete_only():
    """Test complete ORFs only mode"""
    finder = ORFFinder(
        min_protein_length=1,
        complete_orfs_only=True
    )
    
    # Complete ORF
    complete_seq = "ATGGCATAA"  # M A *
    orfs = finder.find_all_orfs(complete_seq, "test")
    
    # Should find ORFs
    assert len(orfs) > 0
    
    # All should be complete
    for orf in orfs:
        if orf.strand == '+':  # Only check forward strand
            assert orf.orf_type == "complete"


def test_orf_finder_3prime_partial_coordinates_match_trimmed_cds():
    finder = ORFFinder(
        min_protein_length=1,
        forward_strand=True,
        reverse_strand=False,
    )

    # No in-frame stop codon after the start. The trailing 2 nt should be
    # trimmed from the CDS and must not remain in the reported coordinates.
    sequence = "ATGAAACCCGGGTTTTT"
    orfs = finder.find_all_orfs(sequence, "test")

    partial = next(orf for orf in orfs if orf.strand == '+' and orf.orf_type == "3prime_partial")

    assert partial.sequence == "ATGAAACCCGGGTTT"
    assert partial.start == 1
    assert partial.end == 15
    assert partial.length == len(partial.sequence)


def test_orf_finder_min_length():
    """Test minimum length filtering"""
    # Very short ORF
    sequence = "ATGGCATAA"  # 9 nt = 3 aa
    
    # Should find with min_protein_length=1
    finder1 = ORFFinder(min_protein_length=1)
    orfs1 = finder1.find_all_orfs(sequence, "test")
    assert len(orfs1) > 0
    
    # Should not find with min_protein_length=10
    finder2 = ORFFinder(min_protein_length=10)
    orfs2 = finder2.find_all_orfs(sequence, "test")
    # May find reverse complement ORFs, so check lengths
    for orf in orfs2:
        assert len(orf.protein) >= 10


def test_orf_model_to_gff3():
    """Test GFF3 output generation"""
    orf = ORF(
        transcript_id="TRINITY_DN1000_c0_g1_i1",
        gene_id="GENE.TEST~~TEST.p1",
        model_id="TEST.p1",
        start=100,
        end=399,
        strand='+',
        sequence="ATG" * 100,  # 300 nt
        protein="M" * 100,
        orf_type="complete",
        length=300
    )
    
    gff3 = orf.to_gff3()
    
    # Check format
    lines = gff3.split('\n')
    assert len(lines) == 3  # gene, mRNA, CDS
    
    # Check that all lines have 9 fields
    for line in lines:
        fields = line.split('\t')
        assert len(fields) == 9
        
    # Check feature types
    assert "gene" in lines[0]
    assert "mRNA" in lines[1]
    assert "CDS" in lines[2]


def test_orf_to_fasta():
    """Test FASTA output generation"""
    orf = ORF(
        transcript_id="TEST",
        gene_id="GENE.TEST~~TEST.p1",
        model_id="TEST.p1",
        start=1,
        end=99,
        strand='+',
        sequence="ATGAAATAA",
        protein="MK*",
        orf_type="complete",
        length=9
    )
    
    # Test CDS FASTA
    cds_fasta = orf.to_fasta_cds()
    assert cds_fasta.startswith(">TEST.p1")
    assert "ATGAAATAA" in cds_fasta
    
    # Test protein FASTA
    pep_fasta = orf.to_fasta_protein()
    assert pep_fasta.startswith(">TEST.p1")
    assert "MK*" in pep_fasta


if __name__ == '__main__':
    pytest.main([__file__, '-v'])


def test_pwm_roundtrip_and_scoring(tmp_path):
    positive = ["CCCCCATGAAAA", "CCCCCATGAAAT", "CCCCCATGAAAG"]
    negative = ["TTTTTATGCCCC", "GGGGGATGTTTT", "AAAACATGCCCC"]

    pwm_plus = build_pwm(positive)
    pwm_minus = build_pwm(negative)

    positive_score = pwm_plus.score_plus_minus(positive[0], pwm_minus)
    negative_score = pwm_plus.score_plus_minus(negative[0], pwm_minus)

    assert positive_score is not None
    assert negative_score is not None
    assert positive_score > negative_score

    pwm_file = tmp_path / "test.pwm"
    pwm_plus.write(pwm_file)
    loaded = PWM.load(pwm_file)
    loaded_score = loaded.score_plus_minus(positive[0], pwm_minus)

    assert loaded_score == pytest.approx(positive_score)


def _predictor_with_one_orf(tmp_path, transcript_seq, mrna_line):
    transcripts_fasta = tmp_path / "transcripts.fa"
    transcripts_fasta.write_text(">tx1\n" + transcript_seq + "\n")
    predictor = TransDecoderPredict(transcripts_file=transcripts_fasta, output_dir=tmp_path)
    predictor.workdir.mkdir(parents=True, exist_ok=True)
    gff3_file = tmp_path / "best_candidates.gff3"
    gff3_file.write_text("##gff-version 3\n" + mrna_line + "\n")
    return predictor, gff3_file


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("segment_contrast,trimmed", [(-5.0, True), (50.0, False)])
def test_open_5prime_end_trimmed_unless_upstream_segment_is_coding(
        tmp_path, monkeypatch, strand, segment_contrast, trimmed):
    """An ORF open at its 5' end starts at its first in-frame ATG when the segment before the
    ATG has little frame contrast (UTR), and stays 5'-partial when it looks coding."""
    from pytransdecoder.core.sequence import reverse_complement
    orf = "GCAGCAGCA" + "ATG" + "AAACCCGGGTTTAAACCC" + "TAA"   # open end, ATG at offset 9
    tx = "CC" + orf + "GG"
    start, end = 3, 2 + len(orf)
    if strand == "-":
        tx = reverse_complement(tx)
        start, end = len(tx) - end + 1, len(tx) - start + 1
    predictor, gff3_file = _predictor_with_one_orf(
        tmp_path, tx,
        f"tx1\ttransdecoder\tmRNA\t{start}\t{end}\t.\t{strand}\t.\tID=tx1.p1;Parent=GENE.tx1~~tx1.p1;score=10.0")
    monkeypatch.setattr(predictor, "_frame_contrast",
                        lambda seq, scores: segment_contrast if len(seq) == 9 else 0.0)

    revised = tmp_path / "revised.gff3"
    n = predictor._revise_start_sites(gff3_file, revised, {"tx1": tx}, {}, coding_contrast=1.0, kozak=None)

    fields = revised.read_text().splitlines()[1].split("\t")
    assert n == (1 if trimmed else 0)
    new_start, new_end = int(fields[3]), int(fields[4])
    if trimmed:
        assert (new_start, new_end) == ((start + 9, end) if strand == "+" else (start, end - 9))
        cds = tx[new_start - 1:new_end]
        assert (cds if strand == "+" else reverse_complement(cds)).startswith("ATG")
        assert "start_revised=true" in fields[8]
    else:
        assert (new_start, new_end) == (start, end)


def test_frame_contrast_matches_orf_selection_frames(tmp_path):
    """The contrast is the in-frame Markov score minus the better of the scores of the sequence
    shifted by one and two bases, as in ORF selection."""
    predictor, _ = _predictor_with_one_orf(tmp_path, "ATGAAATAA", "")
    scores = {"A-0": 1.0, "A-1": -1.0, "A-2": -2.0}
    seq = "AAAAAA"
    s0 = predictor._score_cds_via_markov(seq, 5, scores)
    s1 = predictor._score_cds_via_markov(seq[1:], 5, scores)
    s2 = predictor._score_cds_via_markov(seq[2:], 5, scores)
    assert predictor._frame_contrast(seq, scores) == pytest.approx(s0 - max(s1, s2))


def test_long_orf_auto_accept_only_for_standard_stop_codons(tmp_path):
    """The GC-dependent length table assumes TAA/TAG/TGA; codes with other stop sets (ciliate:
    TGA only) retain no ORF on length alone."""
    transcripts_fasta = tmp_path / "transcripts.fa"
    transcripts_fasta.write_text(">tx1\n" + "ATGAAACCCGGGTTTTAA" * 10 + "\n")
    std = TransDecoderPredict(transcripts_file=transcripts_fasta, output_dir=tmp_path)
    cil = TransDecoderPredict(transcripts_file=transcripts_fasta, output_dir=tmp_path,
                              genetic_code="Tetrahymena")
    assert std._min_length_auto_accept() < 3000
    assert cil._min_length_auto_accept() == float("inf")
    strict = TransDecoderPredict(transcripts_file=transcripts_fasta, output_dir=tmp_path,
                                 genetic_code="Tetrahymena", retain_long_orfs_mode="strict",
                                 retain_long_orfs_length=900)
    assert strict._min_length_auto_accept() == 900


@pytest.mark.parametrize("single_best_only", [True, False])
def test_single_best_orf_is_top_by_coding_score_not_length(tmp_path, single_best_only):
    """One ORF per transcript by default, chosen by coding score before length (a longer
    low-scoring ORF, e.g. antisense, must not win); --retain-multiple-orfs keeps both."""
    transcripts_fasta = tmp_path / "transcripts.fa"
    transcripts_fasta.write_text(">tx1\n" + "ACGT" * 300 + "\n")
    predictor = TransDecoderPredict(transcripts_file=transcripts_fasta, output_dir=tmp_path,
                                    single_best_only=single_best_only)
    predictor.workdir.mkdir(parents=True, exist_ok=True)
    predictor.checkpoints_dir.mkdir(exist_ok=True)
    predictor.gff3_file.write_text(
        "tx1\ttransdecoder\tmRNA\t1\t600\t.\t+\t.\tID=tx1.p1;Parent=GENE.tx1~~tx1.p1\n"
        "tx1\ttransdecoder\tmRNA\t700\t990\t.\t-\t.\tID=tx1.p2;Parent=GENE.tx1~~tx1.p2\n")
    scores = tmp_path / "cds.scores"
    scores.write_text("#acc\tMarkov_order\tseq_length\tscore_1\tscore_2\tscore_3\tscore_4\tscore_5\tscore_6\n"
                      "tx1.p1\t5\t600\t5.0\t1.0\t1.0\t0\t0\t0\n"
                      "tx1.p2\t5\t291\t40.0\t1.0\t1.0\t0\t0\t0\n")
    out = predictor._select_best_orfs(scores)
    ids = [l.split("ID=")[1].split(";")[0] for l in out.read_text().splitlines() if "\tmRNA\t" in l]
    assert ids == (["tx1.p2"] if single_best_only else ["tx1.p2", "tx1.p1"])


def test_cli_reports_single_orf_by_default():
    import argparse
    from pytransdecoder.cli import create_predict_parser
    sub = argparse.ArgumentParser().add_subparsers()
    p = create_predict_parser(sub)
    assert p.parse_args(["-t", "x.fa"]).single_best_only is True
    assert p.parse_args(["-t", "x.fa", "--retain-multiple-orfs"]).single_best_only is False
    assert p.parse_args(["-t", "x.fa", "--single_best_only"]).single_best_only is True


def test_kozak_model_learns_start_context(tmp_path):
    """Training starts with A at -3 and G at +4 vs internal ATGs without them give positive
    log-odds there; the model is written to start_refinement.kozak.tsv."""
    tx, records = [], []
    for k in range(60):
        cds = "ATG" + "GCC" + ("CTGATGCAG" * 4) + "TAA"     # internal in-frame ATGs: TGATGC context
        seq = "CCCCCCACCATGG"[:-4] + cds + "CCCC"           # ...ACC|ATG|G
        tx.append(f">t{k}\n{seq}\n")
        start = seq.index(cds) + 1
        records.append(f">t{k}.p1 t{k}:{start}-{start + len(cds) - 1}(+)\n{cds}\n")
    transcripts_fasta = tmp_path / "transcripts.fa"
    transcripts_fasta.write_text("".join(tx))
    training = tmp_path / "train.cds"
    training.write_text("".join(records))
    predictor = TransDecoderPredict(transcripts_file=transcripts_fasta, output_dir=tmp_path)
    predictor.workdir.mkdir(parents=True, exist_ok=True)
    model = predictor._train_kozak_model(training, predictor._load_transcripts())
    assert model[3]["A"] > 0 > model[3]["C"]       # -3: A at starts, C at internal ATGs
    assert model[6]["G"] > 0                       # +4
    assert (predictor.workdir / "start_refinement.kozak.tsv").exists()
    seq = "CCCACCATGGCC" + "TGATGCAG"
    assert predictor._kozak_score(model, seq, 6) > predictor._kozak_score(model, seq, seq.index("ATG", 9))
    # too few training starts: no model
    training.write_text("".join(records[:10]))
    assert predictor._train_kozak_model(training, predictor._load_transcripts()) is None


@pytest.mark.parametrize("kozak_score,trimmed", [(2.0, True), (-2.0, False)])
def test_kozak_context_shifts_keep_or_trim(tmp_path, monkeypatch, kozak_score, trimmed):
    """With contrast between the two cutoffs, a strong Kozak ATG is taken as the start and a
    weak one leaves the ORF 5'-partial."""
    orf = "GCAGCAGCA" + "ATG" + "AAACCCGGGTTTAAACCC" + "TAA"
    tx = "CC" + orf + "GG"
    predictor, gff3_file = _predictor_with_one_orf(
        tmp_path, tx,
        f"tx1\ttransdecoder\tmRNA\t3\t{2 + len(orf)}\t.\t+\t.\tID=tx1.p1;Parent=GENE.tx1~~tx1.p1;score=10.0")
    # alpha 0.4, weight 0.15: cutoffs 0.7 (kozak +2) and 0.1 (kozak -2); segment ratio 0.4
    monkeypatch.setattr(predictor, "_frame_contrast", lambda seq, scores: 0.4 * len(seq))
    monkeypatch.setattr(predictor, "_kozak_score", lambda model, seq, i: kozak_score)
    revised = tmp_path / "revised.gff3"
    n = predictor._revise_start_sites(gff3_file, revised, {"tx1": tx}, {}, 1.0, kozak=[{}])
    assert n == (1 if trimmed else 0)


def test_final_pep_header_includes_orf_metadata(tmp_path):
    transcripts_fasta = tmp_path / "transcripts.fa"
    transcripts_fasta.write_text(">tx1\nATGAAATAA\n")

    predictor = TransDecoderPredict(transcripts_file=transcripts_fasta, output_dir=tmp_path)

    gff3_file = tmp_path / "best_candidates.gff3"
    gff3_file.write_text(
        "##gff-version 3\n"
        "tx1\ttransdecoder\tmRNA\t1\t9\t.\t+\t.\t"
        "ID=tx1.p1;Parent=GENE.tx1~~tx1.p1;score=10.0;blast=blast:sp|P12345|TEST|1e-5|50\n"
    )

    pep_file = tmp_path / "out.pep"
    predictor._gff3_to_proteins(gff3_file, pep_file, seq_type="pep")

    lines = pep_file.read_text().splitlines()
    assert lines[0] == (
        ">tx1.p1 GENE.tx1~~tx1.p1 ORF type:complete (+),score=10.0,"
        "blast:sp|P12345|TEST|1e-5|50 len:2 tx1:1-9(+)"
    )
    assert lines[1] == "MK*"


def test_pipeline_runs_pfam_search_and_passes_domtblout(tmp_path, monkeypatch):
    transcripts_fasta = tmp_path / "transcripts.fa"
    transcripts_fasta.write_text(">tx1\nATGAAATAA\n")
    pfam_db = tmp_path / "Pfam-A.hmm"
    pfam_db.write_text("HMMER3/f\n")

    recorded_cmds = []
    predict_kwargs = {}

    def fake_run_longorfs(**kwargs):
        output_dir = Path(kwargs["output_dir"]) if kwargs["output_dir"] else Path(".")
        workdir = output_dir / f"{Path(kwargs['transcripts_file']).name}.transdecoder_dir"
        workdir.mkdir(parents=True, exist_ok=True)
        (workdir / "longest_orfs.pep").write_text(">tx1.p1\nMK\n")

    class DummyPredict:
        def __init__(self, **kwargs):
            predict_kwargs.update(kwargs)

        def run(self):
            return None

    def fake_subprocess_run(cmd, *args, **kwargs):
        recorded_cmds.append(cmd)
        if cmd[0] == "hmmpress":
            for ext in (".h3f", ".h3i", ".h3m", ".h3p"):
                Path(str(pfam_db) + ext).write_text("")
        elif cmd[0] == "hmmsearch":
            domtblout = Path(cmd[2])
            domtblout.write_text("# hmmer domtblout\n")
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(cli, "run_longorfs", fake_run_longorfs)
    monkeypatch.setattr("pytransdecoder.predict.TransDecoderPredict", DummyPredict)
    monkeypatch.setattr("subprocess.run", fake_subprocess_run)

    parser = cli.create_pipeline_parser()
    args = parser.parse_args([
        "-t", str(transcripts_fasta),
        "-O", str(tmp_path),
        "--pfam-search-db", str(pfam_db),
    ])

    cli.pipeline_cmd(args)

    assert ["hmmpress", "-f", str(pfam_db)] in recorded_cmds
    assert [
        "hmmsearch",
        "--domtblout",
        str(tmp_path / "transcripts.fa.transdecoder_dir" / "pfam.domtblout"),
        str(pfam_db),
        str(tmp_path / "transcripts.fa.transdecoder_dir" / "longest_orfs.pep"),
    ] in recorded_cmds
    assert predict_kwargs["retain_pfam_hits"] == str(
        tmp_path / "transcripts.fa.transdecoder_dir" / "pfam.domtblout"
    )


def test_score_start_window_bounds_and_normalization():
    """The start-refinement window must cover exactly `up` bases 5' of the ATG, the ATG,
    and `down` bases 3' of it, and be normalized by (up + down), matching training."""
    from pytransdecoder.core.pwm import build_pwm, score_start_window

    atg = 20
    base = "C" * atg + "ATG" + "C" * 10  # 33-nt feature, ATG at index 20

    def with_base(seq, idx, b):
        return seq[:idx] + b + seq[idx + 1:]

    # positives carry an A just outside the window (idx 18) and a T at the last
    # downstream position inside it (idx 23); negatives carry neither.
    plus = build_pwm([with_base(with_base(base, 18, "A"), 23, "T")] * 20)
    minus = build_pwm([base] * 20)

    probe_outside = with_base(base, 18, "A")
    probe_inside = with_base(base, 23, "T")

    # up=1 -> window starts at idx 19, so idx 18 must not contribute
    assert abs(score_start_window(plus, minus, probe_outside, atg, 1, 1) - 
               score_start_window(plus, minus, base, atg, 1, 1)) < 1e-9
    # down=1 -> window ends at idx 23, so idx 23 must contribute
    assert score_start_window(plus, minus, probe_inside, atg, 1, 1) > \
        score_start_window(plus, minus, base, atg, 1, 1)
    # normalization by up + down
    raw = plus.score_plus_minus(probe_inside, minus, pwm_range=(atg - 3, atg + 2 + 4))
    assert abs(score_start_window(plus, minus, probe_inside, atg, 3, 4) - raw / 7) < 1e-9


def test_base_freqs_roundtrip_matches_perl_format(tmp_path):
    """base_freqs.dat must be base/count/ratio (as compute_base_probs.pl writes it) and be
    read back as the Markov background model, not silently replaced by 0.25."""
    from pytransdecoder.core.base_freqs import compute_base_counts, write_base_frequencies, read_base_frequencies

    fa = tmp_path / "t.fa"
    fa.write_text(">a\nAAAAAAACCG\n")
    counts = compute_base_counts(str(fa))
    assert counts == {"A": 7, "T": 7, "C": 3, "G": 3}
    out = tmp_path / "base_freqs.dat"
    write_base_frequencies(counts, str(out))
    assert out.read_text().splitlines()[0] == "A\t7\t0.350"
    assert read_base_frequencies(str(out)) == {"A": 0.35, "C": 0.15, "G": 0.15, "T": 0.35}


def test_orf_finder_min_length_counts_implicit_stop_like_perl():
    """LongOrfs counts the implicit 3' stop as a codon, so a 3'-partial ORF of 99 codons
    passes -m 100 exactly as a complete ORF of 99 codons + stop does."""
    partial = "ATG" + "GCT" * 98                     # 297 nt, runs off the 3' end
    complete = "ATG" + "GCT" * 98 + "TAA"            # 300 nt incl. stop
    too_short = "ATG" + "GCT" * 97                   # 294 nt, runs off the 3' end

    finder = ORFFinder(min_protein_length=100, reverse_strand=False)

    def frame0_types(seq):  # the GCT repeat is also open in frames 1-2; look only at the ATG frame
        return [o.orf_type for o in finder.find_all_orfs(seq, "x") if o.start == 1]

    assert frame0_types(partial) == ["3prime_partial"]
    assert frame0_types(complete) == ["complete"]
    assert frame0_types(too_short) == []


# Translation tables from the Perl TransDecoder's Nuc_translator.pm (codons in TCAG order).
PERL_GENETIC_CODES = {
    'Universal': "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    'Mitochondrial-Vertebrates': "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
    'Mitochondrial-Yeast': "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    'Mitochondrial-Protozoan': "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    'Mitochondrial-Invertebrates': "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
    'Ciliate': 'FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    'Mitochondrial-Echinoderm': "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    'Euplotid': "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    'Candida': "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    'Mitochondrial-Ascidian': "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
    'Mitochondrial-Chlorophycean': "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    'Mitochondrial-Trematode': "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    'Mitochondrial-Scenedesmus_obliquus': "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    'Mitochondrial-Thraustochytrium': "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    'Mitochondrial-Pterobranchia': "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
    'SR1_Gracilibacteria': "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    'Pachysolen_tannophilus': "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    'Mesodinium': "FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    'Peritrich': "FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
}


@pytest.mark.parametrize("name", sorted(PERL_GENETIC_CODES))
def test_genetic_codes_match_perl(name):
    bases = "TCAG"
    codons = [a + b + c for a in bases for b in bases for c in bases]
    t = Translator(name)
    expected = PERL_GENETIC_CODES[name]
    assert "".join(t.translate(c) for c in codons) == expected
    assert set(t.get_stop_codons()) == {c for c, aa in zip(codons, expected) if aa == "*"}


def test_markov_terminal_stop_follows_genetic_code(tmp_path):
    """The terminal stop excluded from Markov scoring must be a stop of the selected code:
    AGA ends a vertebrate mitochondrial ORF; TAA is a glutamine codon under Tetrahymena."""
    transcripts_fasta = tmp_path / "t.fa"
    transcripts_fasta.write_text(">t\nATG\n")
    scores = {"AGA-2": 5.0, "TAA-2": 5.0}  # trinucleotide ending at position 2

    universal = TransDecoderPredict(transcripts_file=transcripts_fasta, output_dir=tmp_path)
    mito = TransDecoderPredict(transcripts_file=transcripts_fasta, output_dir=tmp_path,
                               genetic_code="Mitochondrial-Vertebrates")
    tetra = TransDecoderPredict(transcripts_file=transcripts_fasta, output_dir=tmp_path,
                                genetic_code="Tetrahymena")

    # 3-nt sequences: position 0 starts the terminal triplet, so a stop scores nothing
    assert universal._score_cds_via_markov("AGA", 2, scores) != 0.0
    assert mito._score_cds_via_markov("AGA", 2, scores) == 0.0
    assert universal._score_cds_via_markov("TAA", 2, scores) == 0.0
    assert tetra._score_cds_via_markov("TAA", 2, scores) != 0.0
