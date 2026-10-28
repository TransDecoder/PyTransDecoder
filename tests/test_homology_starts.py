"""Homology-assisted start refinement (core/homology_starts.py and Predict._revise_start_sites)."""
from pathlib import Path

from pytransdecoder.core.homology_starts import homology_start_call
from pytransdecoder.predict import TransDecoderPredict

# --- rules -------------------------------------------------------------------------------------
# first in-frame ATG at peptide position 21 in all cases


def test_r2_trim_subject_nterm_at_atg():
    assert homology_start_call(21, blast=(21, 1)) == (["R2"], "trim")
    assert homology_start_call(21, blast=(23, 3)) == (["R2"], "trim")      # sstart 3, maps to 21
    assert homology_start_call(21, blast=(26, 1)) == (["R2"], "trim")      # +5 aa
    assert homology_start_call(21, blast=(16, 1)) == (["R2"], "trim")      # -5 aa
    assert homology_start_call(21, blast=(27, 1)) == ([], None)            # +6 aa: too far
    assert homology_start_call(21, blast=(25, 4)) == ([], None)            # sstart 4: not the N-terminus


def test_r1_keep_alignment_upstream():
    assert homology_start_call(21, blast=(16, 40)) == (["R1"], "keep")     # 5 aa upstream, subject continues
    assert homology_start_call(21, blast=(17, 40)) == ([], None)           # only 4 aa upstream
    assert homology_start_call(21, blast=(1, 1)) == ([], None)             # sstart 1: subject ends there too
    assert homology_start_call(21, blast=(1, 2)) == (["R1"], "keep")


def test_r3_keep_pfam_envelope_upstream():
    assert homology_start_call(21, pfam_env_from=20) == (["R3"], "keep")
    assert homology_start_call(21, pfam_env_from=21) == ([], None)
    assert homology_start_call(21, blast=None, pfam_env_from=None) == ([], None)


def test_conflict_and_combination():
    # R2 (subject start at the ATG) and R3 (domain envelope upstream) disagree -> conflict
    assert homology_start_call(21, blast=(21, 1), pfam_env_from=3) == (["R2", "R3"], "conflict")
    # R1 and R3 agree
    assert homology_start_call(21, blast=(2, 30), pfam_env_from=3) == (["R1", "R3"], "keep")
    # R1 and R2 cannot both fire on one hit (R1 needs sstart > 1 and qstart <= a - 5, so the subject
    # start maps to <= a - 6): at the boundary only R1 fires
    assert homology_start_call(21, blast=(16, 2)) == (["R1"], "keep")


# --- integration with the start refinement ------------------------------------------------------
UTR = "CCGCCA" * 10  # 60 nt, no ATG or stop in frame 0
CDS = "ATG" + "GCT" * 60 + "TAA"
TX = UTR + CDS


def _run(tmp_path: Path, blast_lines=None, pfam_lines=None, assist=True):
    fa = tmp_path / "t.fa"
    fa.write_text(">tx1\n" + TX + "\n")
    blast = pfam = None
    if blast_lines is not None:
        blast = tmp_path / "blastp.outfmt6"
        blast.write_text("".join(l + "\n" for l in blast_lines))
    if pfam_lines is not None:
        pfam = tmp_path / "pfam.domtblout"
        pfam.write_text("# comment\n" + "".join(l + "\n" for l in pfam_lines))
    p = TransDecoderPredict(transcripts_file=fa, output_dir=tmp_path, retain_blastp_hits=blast,
                            retain_pfam_hits=pfam, homology_start_assist=assist)
    p.workdir.mkdir(parents=True, exist_ok=True)
    gff = tmp_path / "in.gff3"
    gff.write_text(f"tx1\ttransdecoder\tmRNA\t1\t{len(TX)}\t.\t+\t.\tID=tx1.p1;Parent=g1\n")
    out = tmp_path / "out.gff3"
    homology = p._homology_start_evidence()
    # empty hexamer scores: coding score 0 < UTR score 0.4 * 60 -> the ab initio rule trims
    n = p._revise_start_sites(gff, out, {"tx1": TX}, {}, 1.0, None, homology)
    line = out.read_text().strip()
    return n, line, p.workdir / "start_refinement.homology.tsv"


def _blast(qstart, sstart):
    return f"tx1.p1\tsp|X|Y\t90\t50\t1\t0\t{qstart}\t{qstart + 49}\t{sstart}\t{sstart + 49}\t1e-20\t100"


def _pfam(env_from):
    f = ["tx1.p1", "-", "81", "Dom", "PF00001.1", "100", "1e-10", "40", "0", "1", "1", "1e-10", "1e-10",
         "40", "0", "1", "50", str(env_from), "60", str(env_from), "60", "0.9", "desc"]
    return " ".join(f)


def test_ab_initio_without_hits(tmp_path):
    n, line, table = _run(tmp_path)
    assert n == 1 and line.split("\t")[3] == "61" and "start_homology" not in line   # trimmed to the ATG
    assert not table.exists()


def test_r1_keeps_open(tmp_path):
    n, line, table = _run(tmp_path, blast_lines=[_blast(2, 30)])   # ATG at peptide position 21
    assert n == 0 and line.split("\t")[3] == "1"
    rows = table.read_text().splitlines()
    assert rows[1].split("\t")[8:] == ["R1", "0.0", "24.0", "trim", "keep", "homology"]


def test_conflict_falls_back_to_ab_initio(tmp_path):
    n, line, table = _run(tmp_path, blast_lines=[_blast(21, 1)], pfam_lines=[_pfam(3)])
    assert n == 1 and "start_homology" not in line
    assert table.read_text().splitlines()[1].endswith("R2,R3\t0.0\t24.0\ttrim\ttrim\tab_initio(conflict)")


def test_r2_trim_recorded_and_opt_out(tmp_path):
    n, line, _ = _run(tmp_path, blast_lines=[_blast(21, 1)])
    assert n == 1 and "start_homology=R2" in line
    sub = tmp_path / "optout"
    sub.mkdir()
    n, line, table = _run(sub, blast_lines=[_blast(2, 30)], assist=False)
    assert n == 1 and not table.exists()      # opt-out: ab initio trims despite R1


def test_blast_only_ignores_missing_pfam(tmp_path):
    n, line, table = _run(tmp_path, pfam_lines=[_pfam(3)])   # Pfam only: R3 keeps open
    assert n == 0
    assert table.read_text().splitlines()[1].split("\t")[5:9] == [".", ".", "3", "R3"]
