"""
Basic tests for PyTransDecoder core modules
"""

import pytest
from pytransdecoder.core.translator import Translator
from pytransdecoder.core.sequence import reverse_complement
from pytransdecoder.core.orf_finder import ORFFinder
from pytransdecoder.core.models import ORF


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
