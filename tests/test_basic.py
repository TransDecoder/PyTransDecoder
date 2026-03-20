"""
Basic tests for PyTransDecoder core modules
"""

import pytest
from pytransdecoder.core.translator import Translator
from pytransdecoder.core.sequence import reverse_complement
from pytransdecoder.core.orf_finder import ORFFinder
from pytransdecoder.core.models import ORF
from pytransdecoder.core.pwm import build_pwm, PWM
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


def test_refine_start_sites_updates_plus_strand(tmp_path):
    transcript_seq = "CCCCCCCCCCCCCCCCCCCCCCTGAAACCCGGGTTTAAA"
    transcript_seq = transcript_seq[:21] + "ATG" + transcript_seq[24:]
    transcripts_fasta = tmp_path / "transcripts.fa"
    transcripts_fasta.write_text(">tx1\n" + transcript_seq + "\n")

    predictor = TransDecoderPredict(transcripts_file=transcripts_fasta, output_dir=tmp_path)
    predictor.workdir.mkdir(parents=True, exist_ok=True)

    gff3_file = tmp_path / "best_candidates.gff3"
    gff3_file.write_text(
        "##gff-version 3\n"
        "tx1\ttransdecoder\tmRNA\t1\t34\t.\t+\t.\t"
        "ID=tx1.p1;Parent=GENE.tx1~~tx1.p1;Name=tx1.p1;score=10.0\n"
    )

    revised_gff3 = tmp_path / "best_candidates.gff3.revised"
    pwm_plus = build_pwm(["CCCCCCCCCCCCCCCCCCCCCATGAAACCCGGGTT"])
    pwm_minus = build_pwm(
        [
            "TTTTTTTTTTTTTTTTTTTTTATGCCCCGGGAAAA",
            "GGGGGGGGGGGGGGGGGGGGGATGTTTTCCCAAAA",
            "AAAACCCCAAAACCCCAAAACCATGCCCCAAAACC",
        ]
    )

    num_revised = predictor._revise_start_sites(
        gff3_file,
        revised_gff3,
        {"tx1": transcript_seq},
        pwm_plus,
        pwm_minus,
        (5, 5),
        min_threshold=0.0,
    )

    assert num_revised == 1
    revised_line = revised_gff3.read_text().splitlines()[1]
    fields = revised_line.split("\t")
    assert fields[3] == "22"
    assert "start_revised=true" in fields[8]


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
    assert lines[1] == "MK"
