"""
Start-codon (Kozak) context model used by start codon refinement.

Positions -6..-1 and +4..+5 around an ATG (A of ATG = +1; the ATG itself is not scored).
Log2-odds (bits, +1 pseudocounts) of training ORF start contexts vs in-frame ATGs inside the
training ORFs: +1 bit = twice as start-like as an average internal ATG. Same model and file formats as Perl TransDecoder (PerlLib/Kozak_model.pm,
util/train_kozak_model.pl).
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Dict, List, Optional, Sequence

UPSTREAM = 6
DOWNSTREAM = 2
BASES = ("A", "C", "G", "T")
MIN_TRAINING_STARTS = 50

Counts = List[Dict[str, int]]
Model = List[Dict[str, float]]


def position_labels() -> List[str]:
    return [str(k) for k in range(-UPSTREAM, 0)] + [f"+{k}" for k in range(4, 4 + DOWNSTREAM)]


def window(sequence: str, atg_index: int) -> List[Optional[str]]:
    """Bases at -6..-1 and +4..+5 around the ATG at atg_index (None outside the sequence)."""
    positions = list(range(atg_index - UPSTREAM, atg_index)) + \
        list(range(atg_index + 3, atg_index + 3 + DOWNSTREAM))
    return [sequence[j] if 0 <= j < len(sequence) and sequence[j] in "ACGT" else None
            for j in positions]


def new_counts() -> Counts:
    return [dict.fromkeys(BASES, 0) for _ in range(UPSTREAM + DOWNSTREAM)]


def add_site(counts: Counts, sequence: str, atg_index: int) -> None:
    for k, base in enumerate(window(sequence, atg_index)):
        if base:
            counts[k][base] += 1


def log_odds_model(pos_counts: Counts, neg_counts: Counts) -> Model:
    model = []
    for pc, nc in zip(pos_counts, neg_counts):
        pt, nt = sum(pc.values()) + 4, sum(nc.values()) + 4
        model.append({b: math.log2(((pc[b] + 1) / pt) / ((nc[b] + 1) / nt)) for b in BASES})
    return model


def score(model: Model, sequence: str, atg_index: int) -> float:
    return sum(model[k][b] for k, b in enumerate(window(sequence, atg_index)) if b)


def write_model(model: Model, num_starts: int, filename: Path) -> None:
    with open(filename, "w") as out:
        out.write("#position\tA\tC\tG\tT  (log2-odds in bits, start vs internal in-frame ATG; "
                  f"{num_starts} training starts)\n")
        for label, row in zip(position_labels(), model):
            out.write(label + "\t" + "\t".join(f"{row[b]:.3f}" for b in BASES) + "\n")


def read_model(filename: Path) -> Model:
    model = []
    with open(filename) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 5:
                model.append({b: float(v) for b, v in zip(BASES, parts[1:])})
    if len(model) != UPSTREAM + DOWNSTREAM:
        raise ValueError(f"{filename} does not have {UPSTREAM + DOWNSTREAM} positions")
    return model


def write_counts(pos_counts: Counts, neg_counts: Counts, filename: Path) -> None:
    with open(filename, "w") as out:
        out.write("\t".join(["#position"] + [f"start_{b}" for b in BASES] +
                            [f"internal_{b}" for b in BASES]) + "\n")
        for label, pc, nc in zip(position_labels(), pos_counts, neg_counts):
            out.write("\t".join([label] + [str(pc[b]) for b in BASES] + [str(nc[b]) for b in BASES]) + "\n")


def _nice_ceiling(value: float) -> float:
    for step in (0.25, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10):
        if value <= step:
            return step
    return int(value) + 1


def _fmt(v: float) -> str:
    return f"{v:g}"


def write_logo_svg(pos_counts: Counts, model: Model, num_starts: int, filename: Path) -> None:
    """
    Two panels over -6..-1, ATG, +4..+5 (the ATG itself is labeled, not drawn):
      top:    information content (bits) of the training start contexts
      bottom: log2-odds vs internal in-frame ATGs, positive above and negative below the axis
    Both y scales fit the data (at least 0.5 bits / 1 log2 unit).
    """
    color = {"A": "#109648", "C": "#255C99", "G": "#F7B32B", "T": "#D62839"}
    glyph_w = {"A": 72, "C": 72, "G": 78, "T": 61}   # bold Helvetica/Arial at 100px; cap height ~72px
    labels = position_labels()
    cols = labels[:UPSTREAM] + ["A", "T", "G"] + labels[UPSTREAM:]
    colw, left, top_h, bot_h, gap = 40, 56, 150, 220, 56

    ic_stack, up_stack, dn_stack = [], [], []
    max_ic, max_lo = 0.5, 1.0
    for k in range(UPSTREAM + DOWNSTREAM):
        n = sum(pos_counts[k].values())
        f = {b: (pos_counts[k][b] + 0.5) / (n + 2) for b in BASES}
        ic = 2 + sum(f[b] * math.log2(f[b]) for b in BASES)
        ic_stack.append([(b, f[b] * ic) for b in sorted(BASES, key=lambda b: f[b])])
        max_ic = max(max_ic, ic)
        lo = model[k]
        up = [(b, lo[b]) for b in sorted(BASES, key=lambda b: lo[b]) if lo[b] > 0]
        dn = [(b, -lo[b]) for b in sorted(BASES, key=lambda b: -lo[b]) if lo[b] < 0]
        up_stack.append(up)
        dn_stack.append(dn)
        max_lo = max(max_lo, sum(h for _, h in up), sum(h for _, h in dn))
    ic_top, lo_top = _nice_ceiling(max_ic), _nice_ceiling(max_lo)

    top_y0 = 40 + top_h
    bot_axis = top_y0 + gap + bot_h / 2
    height = top_y0 + gap + bot_h + 40
    width = left + len(cols) * colw + 20
    ic_px, lo_px = top_h / ic_top, (bot_h / 2) / lo_top

    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{_fmt(height)}" '
           'font-family="Helvetica,Arial,sans-serif">',
           '<rect width="100%" height="100%" fill="white"/>',
           f'<text x="{left}" y="22" font-size="14">Start-codon context of {num_starts} training ORF starts (bits)</text>',
           f'<text x="{left}" y="{_fmt(top_y0 + gap - 14)}" font-size="14">Log2-odds vs in-frame ATGs inside the training ORFs</text>']

    def letter(b: str, x: float, y_base: float, h: float) -> None:
        if h < 0.5:
            return
        svg.append(f'<text x="0" y="0" font-size="100" font-weight="bold" fill="{color[b]}" '
                   f'transform="translate({x + 2:.1f},{y_base:.1f}) scale({(colw - 4) / glyph_w[b]:.4f},{h / 72:.4f})">{b}</text>')

    for c, label in enumerate(cols):
        x = left + c * colw
        is_atg = UPSTREAM <= c < UPSTREAM + 3
        style = ' font-weight="bold" fill="#888888"' if is_atg else ""
        svg.append(f'<text x="{x + colw / 2:.1f}" y="{height - 14:.1f}" font-size="{16 if is_atg else 12}" '
                   f'text-anchor="middle"{style}>{label}</text>')
        if is_atg:
            continue
        k = c if c < UPSTREAM else c - 3
        y = top_y0
        for b, v in ic_stack[k]:
            letter(b, x, y, v * ic_px)
            y -= v * ic_px
        y = bot_axis
        for b, v in up_stack[k]:
            letter(b, x, y, v * lo_px)
            y -= v * lo_px
        y = bot_axis
        for b, v in dn_stack[k]:
            y += v * lo_px
            letter(b, x, y, v * lo_px)

    xr = left + len(cols) * colw
    svg.append(f'<line x1="{left}" y1="{_fmt(top_y0)}" x2="{xr}" y2="{_fmt(top_y0)}" stroke="black"/>')
    svg.append(f'<line x1="{left}" y1="{_fmt(bot_axis)}" x2="{xr}" y2="{_fmt(bot_axis)}" stroke="black"/>')
    svg.append(f'<line x1="{left - 4}" y1="{top_y0 - top_h:.1f}" x2="{left - 4}" y2="{top_y0:.1f}" stroke="black"/>')
    svg.append(f'<line x1="{left - 4}" y1="{bot_axis - bot_h / 2:.1f}" x2="{left - 4}" y2="{bot_axis + bot_h / 2:.1f}" stroke="black"/>')
    for text, y in [(0, top_y0), (_fmt(ic_top / 2), top_y0 - top_h / 2), (_fmt(ic_top), top_y0 - top_h),
                    (f"+{_fmt(lo_top)}", bot_axis - bot_h / 2), (0, bot_axis), (f"-{_fmt(lo_top)}", bot_axis + bot_h / 2)]:
        svg.append(f'<text x="{left - 8}" y="{y + 4:.1f}" font-size="11" text-anchor="end">{text}</text>')
    svg.append("</svg>")

    with open(filename, "w") as out:
        out.write("\n".join(svg) + "\n")
