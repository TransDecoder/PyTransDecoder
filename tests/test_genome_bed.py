"""genome-map writes BED12 next to the genome GFF3 (as Perl TransDecoder's genome mode)."""

from pytransdecoder.core.gene_model import CodingGene, CodingModel
from pytransdecoder.genome import run


def test_genome_map_writes_bed(tmp_path):
    # tx1 (30 nt) = chr1:101-115 + chr1:201-215 (+ strand); ORF at transcript 4-27
    (tmp_path / "tx.fa").write_text(">tx1\n" + "A" * 30 + "\n")
    (tmp_path / "tx.gtf").write_text(
        'chr1\tx\ttranscript\t101\t215\t.\t+\t.\tgene_id "g1"; transcript_id "tx1";\n'
        'chr1\tx\texon\t101\t115\t.\t+\t.\tgene_id "g1"; transcript_id "tx1";\n'
        'chr1\tx\texon\t201\t215\t.\t+\t.\tgene_id "g1"; transcript_id "tx1";\n')
    orf = CodingModel("tx1", "+", [(1, 30)], 4, 27, "tx1.p1", name="ORF type:complete (+),score=9.9")
    (tmp_path / "orfs.gff3").write_text(
        CodingGene(gene_id="GENE.tx1~~tx1.p1", models=[orf]).gff3_text())

    out_gff3 = tmp_path / "orfs.genome.gff3"
    n_done, n_total = run(tmp_path / "orfs.gff3", tmp_path / "tx.gtf", tmp_path / "tx.fa", out_gff3)
    assert (n_done, n_total) == (1, 1)
    fields = (tmp_path / "orfs.genome.bed").read_text().rstrip("\n").split("\t")
    assert fields[:3] == ["chr1", "100", "215"]
    assert fields[3].startswith("ID=tx1.p1;")
    assert fields[5:] == ["+", "103", "212", "0", "2", "15,15", "0,100"]
