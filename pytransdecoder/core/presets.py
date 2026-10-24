"""
Presets: models and ORF acceptance thresholds learned from a large reference transcriptome, for inputs
too small (or too partial) to learn from (experimental).

A preset directory holds
  preset.json            name, source, genetic code, GC %, long-ORF auto-accept length (nt or null),
                         decoy strata, q levels
  hexamer.scores         the reference run's Markov coding model
  training_orfs.cds      the reference run's training ORFs (fixes the median per-base coding contrast)
  kozak.tsv              the reference run's start-context model (optional)
  thresholds.tsv         per decoy stratum and q level, the smallest frame contrast accepted in the
                         reference run's target-decoy q-values (decoy_qvalues.tsv)
With a preset, Predict scores the input with the preset model and keeps a candidate ORF (in-frame score
> 0, frame contrast > 0) when its contrast reaches the preset threshold of its stratum at the chosen q.
"""
import json
import math
import shutil
from collections import defaultdict
from pathlib import Path
from typing import Dict, Optional, Tuple

Q_LEVELS = (0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5)


def build_preset(run_workdir: Path, out_dir: Path, name: str, gc_pct: float,
                 auto_accept_nt: Optional[float], genetic_code: str, strata: str, source: str = "") -> Path:
    """Write a preset from a finished decoy-rule Predict run (its <input>.transdecoder_dir)."""
    run_workdir, out_dir = Path(run_workdir), Path(out_dir)
    qfile = run_workdir / "decoy_qvalues.tsv"
    if not qfile.exists():
        raise FileNotFoundError(f"{qfile} not found: build presets from a --short-orf-rule decoy run")
    out_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy(run_workdir / "hexamer.scores", out_dir / "hexamer.scores")
    shutil.copy(run_workdir / "__checkpoints_TDpredict" / "top_training_orfs.cds", out_dir / "training_orfs.cds")
    kozak = run_workdir / "start_refinement.kozak.tsv"
    if kozak.exists():
        shutil.copy(kozak, out_dir / "kozak.tsv")
    by_stratum = defaultdict(list)
    with open(qfile) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            _, stratum, contrast, q = line.rstrip("\n").split("\t")
            by_stratum[stratum].append((float(q), float(contrast)))
    with open(out_dir / "thresholds.tsv", "w") as fh:
        fh.write("#stratum\tq\tmin_contrast\tn_targets\n")
        for stratum in sorted(by_stratum, key=lambda s: tuple(_key(s))):
            rows = by_stratum[stratum]
            for level in Q_LEVELS:
                ok = [c for q, c in rows if q <= level]
                fh.write(f"{stratum}\t{level:g}\t{min(ok) if ok else 'inf'}\t{len(rows)}\n")
    meta = {"name": name, "source": source, "genetic_code": genetic_code, "gc_pct": round(gc_pct, 2),
            "auto_accept_nt": None if auto_accept_nt in (None, float("inf")) else auto_accept_nt,
            "strata": strata, "q_levels": list(Q_LEVELS)}
    (out_dir / "preset.json").write_text(json.dumps(meta, indent=1) + "\n")
    return out_dir


def _key(stratum: str) -> Tuple:
    out = []
    for x in stratum.split(","):
        out.append(x == "True" if x in ("True", "False") else int(x))
    return tuple(out)


class Preset:
    def __init__(self, path):
        self.path = Path(path)
        self.meta = json.loads((self.path / "preset.json").read_text())
        self.hexamer_scores = self.path / "hexamer.scores"
        self.training_orfs = self.path / "training_orfs.cds"
        self.kozak = self.path / "kozak.tsv" if (self.path / "kozak.tsv").exists() else None
        self.strata = self.meta["strata"]
        self._thr: Dict[float, Dict[Tuple, float]] = defaultdict(dict)
        with open(self.path / "thresholds.tsv") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                stratum, q, c, _ = line.rstrip("\n").split("\t")
                self._thr[float(q)][_key(stratum)] = float(c)

    @property
    def auto_accept_nt(self) -> float:
        a = self.meta.get("auto_accept_nt")
        return float("inf") if a is None else float(a)

    def level(self, q: float) -> float:
        """Largest stored q level <= q (the smallest level if q is below all of them)."""
        levels = sorted(self._thr)
        ok = [x for x in levels if x <= q + 1e-12]
        return ok[-1] if ok else levels[0]

    def threshold(self, stratum: Tuple, q: float) -> float:
        """Minimum frame contrast at q for a stratum; strata absent from the reference use the nearest
        ORF-length bin with the same remaining keys, then any transcript-length bin of the same type."""
        table = self._thr[self.level(q)]
        if stratum in table:
            return table[stratum]
        rest = stratum[1:]
        same = [k for k in table if k[1:] == rest]
        if not same and len(stratum) == 3:
            same = [k for k in table if k[2] == stratum[2]]
        if not same:
            same = list(table)
        best = min(same, key=lambda k: (abs(k[0] - stratum[0]), abs(k[1] - stratum[1]) if len(k) > 1 else 0))
        return table[best]
