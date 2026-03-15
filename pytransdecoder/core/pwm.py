"""
Position weight matrix utilities for start codon refinement.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


BASE_ORDER = ("G", "A", "T", "C")
PSEUDOCOUNT = 0.1


@dataclass
class PWM:
    """Minimal PWM implementation matching the Perl PWM.pm behavior."""

    pos_freqs: List[Dict[str, int]] = field(default_factory=list)
    pos_probs: List[Dict[str, float]] = field(default_factory=list)
    built: bool = False

    def add_feature(self, feature_seq: str) -> None:
        feature_seq = feature_seq.upper()
        if any(base not in BASE_ORDER for base in feature_seq):
            return

        if self.pos_freqs and len(feature_seq) != self.length:
            raise ValueError(
                f"PWM length {self.length} does not match feature length {len(feature_seq)}"
            )

        while len(self.pos_freqs) < len(feature_seq):
            self.pos_freqs.append({})

        for index, base in enumerate(feature_seq):
            self.pos_freqs[index][base] = self.pos_freqs[index].get(base, 0) + 1

        self.built = False

    def remove_feature(self, feature_seq: str) -> None:
        feature_seq = feature_seq.upper()
        if len(feature_seq) != self.length:
            raise ValueError(
                f"PWM length {self.length} does not match feature length {len(feature_seq)}"
            )

        for index, base in enumerate(feature_seq):
            current = self.pos_freqs[index].get(base, 0)
            self.pos_freqs[index][base] = current - 1

        self.built = False

    @property
    def length(self) -> int:
        if self.pos_probs:
            return len(self.pos_probs)
        return len(self.pos_freqs)

    def build(self) -> None:
        self.pos_probs = []
        for pos_freq in self.pos_freqs:
            total = sum(pos_freq.values())
            pos_probs = {}
            for base in BASE_ORDER:
                count = pos_freq.get(base, 0)
                pos_probs[base] = (count + PSEUDOCOUNT) / (total + 4 * PSEUDOCOUNT)
            self.pos_probs.append(pos_probs)

        self.built = True

    def write(self, filename: Path) -> None:
        if not self.built:
            self.build()

        with open(filename, "w") as handle:
            handle.write("pos\tG\tA\tT\tC\n")
            for index, pos_probs in enumerate(self.pos_probs):
                handle.write(
                    f"{index}\t{pos_probs['G']:.6f}\t{pos_probs['A']:.6f}\t"
                    f"{pos_probs['T']:.6f}\t{pos_probs['C']:.6f}\n"
                )

    @classmethod
    def load(cls, filename: Path) -> "PWM":
        pwm = cls()
        with open(filename) as handle:
            header = next(handle, None)
            if header is None:
                raise ValueError(f"PWM file is empty: {filename}")
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) != 5:
                    continue
                pwm.pos_probs.append(
                    {
                        "G": float(parts[1]),
                        "A": float(parts[2]),
                        "T": float(parts[3]),
                        "C": float(parts[4]),
                    }
                )
        pwm.built = True
        return pwm

    def score_plus_minus(
        self,
        target_sequence: str,
        pwm_minus: "PWM",
        pwm_range: Optional[Tuple[int, int]] = None,
    ) -> Optional[float]:
        target_sequence = target_sequence.upper()
        if any(base not in BASE_ORDER for base in target_sequence):
            return None
        if not self.built or not pwm_minus.built:
            raise ValueError("Both PWMs must be built before scoring")
        if len(target_sequence) != self.length or self.length != pwm_minus.length:
            raise ValueError("Target sequence and PWMs must have matching lengths")

        start, end = (0, self.length - 1) if pwm_range is None else pwm_range
        motif_score = 0.0
        for index in range(start, end + 1):
            base = target_sequence[index]
            pos_prob = self.pos_probs[index].get(base)
            neg_prob = pwm_minus.pos_probs[index].get(base)
            if not pos_prob or not neg_prob:
                return None
            motif_score += math.log(pos_prob / neg_prob)

        return motif_score


def build_pwm(features: Iterable[str]) -> PWM:
    pwm = PWM()
    for feature in features:
        pwm.add_feature(feature)
    pwm.build()
    return pwm


def trapezoid_auc(points: Sequence[Tuple[float, float]]) -> float:
    ordered = sorted(points, key=lambda item: (item[0], item[1]))
    ordered = [(0.0, 0.0)] + ordered + [(1.0, 1.0)]
    total_auc = 0.0
    prev_x, prev_y = ordered[0]
    for curr_x, curr_y in ordered[1:]:
        total_auc += ((curr_x - prev_x) * min(prev_y, curr_y)) + (
            0.5 * (curr_x - prev_x) * abs(curr_y - prev_y)
        )
        prev_x, prev_y = curr_x, curr_y
    return total_auc
