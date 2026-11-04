"""External programs (gffread, DIAMOND/BLAST+, HMMER, samtools) are checked before any work starts."""

import stat
import subprocess
import sys
from pathlib import Path

import pytest

from pytransdecoder.external_tools import missing_programs

TX = ">tx1\n" + "ATG" + "GCT" * 120 + "TAA" + "\n"


def _fake_programs(bin_dir: Path, names):
    bin_dir.mkdir(exist_ok=True)
    for n in names:
        p = bin_dir / n
        p.write_text("#!/bin/sh\nexit 0\n")
        p.chmod(p.stat().st_mode | stat.S_IXUSR)


def _pipeline(tmp_path, path_dir, *args):
    """Run pyTransdecoder with PATH limited to path_dir (the Python interpreter is called by full path)."""
    code = "import sys; from pytransdecoder.cli import pipeline_main; sys.argv[0] = 'pyTransdecoder'; pipeline_main()"
    env = {"PATH": str(path_dir), "PYTHONPATH": str(Path(__file__).resolve().parents[1]), "HOME": str(tmp_path)}
    return subprocess.run([sys.executable, "-c", code, *args], cwd=tmp_path, env=env,
                          capture_output=True, text=True)


def test_missing_programs(tmp_path, monkeypatch):
    _fake_programs(tmp_path / "bin", ["diamond"])
    monkeypatch.setenv("PATH", str(tmp_path / "bin"))
    needed = [("diamond", "a"), ("hmmsearch", "b"), ("hmmsearch", "c"), ("samtools", "d")]
    assert missing_programs(needed) == [("hmmsearch", "b"), ("samtools", "d")]


@pytest.mark.parametrize("args, missing", [
    (["--blast-search-pep", "db.pep"], ["diamond"]),
    (["--blast-search-pep", "db.pep", "--blast-tool", "blastp"], ["makeblastdb", "blastp"]),
    (["--pfam-search-db", "Pfam-A.hmm"], ["hmmsearch", "hmmpress"]),
    (["--blast-search-pep", "db.pep", "--pfam-search-db", "Pfam-A.hmm"], ["diamond", "hmmsearch", "hmmpress"]),
])
def test_pipeline_fails_before_longorfs(tmp_path, args, missing):
    (tmp_path / "t.fa").write_text(TX)
    (tmp_path / "db.pep").write_text(">p\nMA\n")
    (tmp_path / "Pfam-A.hmm").write_text("")
    empty = tmp_path / "empty_bin"
    empty.mkdir()
    r = _pipeline(tmp_path, empty, "-t", "t.fa", *args)
    assert r.returncode == 1
    assert "required program(s) not found on PATH" in r.stderr
    for m in missing:
        assert f"  {m}  (needed for" in r.stderr
    assert "Traceback" not in r.stderr
    assert not (tmp_path / "t.fa.transdecoder_dir").exists()   # nothing ran


def test_pressed_pfam_db_does_not_need_hmmpress(tmp_path):
    (tmp_path / "t.fa").write_text(TX)
    for ext in ("", ".h3f", ".h3i", ".h3m", ".h3p"):
        (tmp_path / f"Pfam-A.hmm{ext}").write_text("")
    empty = tmp_path / "empty_bin"
    empty.mkdir()
    r = _pipeline(tmp_path, empty, "-t", "t.fa", "--pfam-search-db", "Pfam-A.hmm")
    assert r.returncode == 1
    assert "  hmmsearch  (needed for" in r.stderr and "hmmpress" not in r.stderr


def test_genome_mode_needs_gffread(tmp_path):
    (tmp_path / "g.fa").write_text(">chr1\nACGT\n")
    (tmp_path / "a.gtf").write_text("")
    empty = tmp_path / "empty_bin"
    empty.mkdir()
    r = _pipeline(tmp_path, empty, "--genome", "g.fa", "--gtf", "a.gtf")
    assert r.returncode == 1 and "  gffread  (needed for --genome/--gtf)" in r.stderr


def test_bam_structures_need_samtools(tmp_path):
    (tmp_path / "t.fa").write_text(TX)
    (tmp_path / "s.bam").write_bytes(b"")
    empty = tmp_path / "empty_bin"
    empty.mkdir()
    r = _pipeline(tmp_path, empty, "-t", "t.fa", "--structures", "s.bam")
    assert r.returncode == 1 and "  samtools  (needed for" in r.stderr
    assert not (tmp_path / "t.fa.transdecoder_dir").exists()


def test_programs_present_pass_the_check(tmp_path):
    # fake diamond is on PATH: the check passes and the run proceeds to LongOrfs (the fake search then
    # writes no hits file, which is not this test's concern)
    (tmp_path / "t.fa").write_text(TX)
    (tmp_path / "db.pep").write_text(">p\nMA\n")
    _fake_programs(tmp_path / "bin", ["diamond"])
    r = _pipeline(tmp_path, tmp_path / "bin", "-t", "t.fa", "--blast-search-pep", "db.pep")
    assert "required program(s) not found" not in r.stderr
    assert (tmp_path / "t.fa.transdecoder_dir" / "longest_orfs.pep").exists()


def test_genome_map_bam_needs_samtools(tmp_path):
    for f in ("orfs.gff3", "t.fa"):
        (tmp_path / f).write_text("")
    (tmp_path / "s.bam").write_bytes(b"")
    empty = tmp_path / "empty_bin"
    empty.mkdir()
    code = "import sys; from pytransdecoder.cli import main; sys.argv[0] = 'pytransdecoder'; main()"
    env = {"PATH": str(empty), "PYTHONPATH": str(Path(__file__).resolve().parents[1]), "HOME": str(tmp_path)}
    r = subprocess.run([sys.executable, "-c", code, "genome-map", "--orfs", "orfs.gff3", "--structures", "s.bam",
                        "-t", "t.fa"], cwd=tmp_path, env=env, capture_output=True, text=True)
    assert r.returncode == 1 and "  samtools  (needed for" in r.stderr
