"""
Homology-assisted start refinement (used only when BLASTP and/or Pfam hits are supplied to Predict).

For an ORF that is open at its 5' end and has a first in-frame ATG at peptide position a (1-based
position in the longest-ORF peptide, whose position 1 is the first whole codon at the open end):

  R2 (BLASTP, trim):  the best hit's subject N-terminus is aligned (sstart <= R2_MAX_SSTART) and maps
                      to within R2_MAX_OFFSET aa of the ATG: |qstart - (sstart - 1) - a| <= R2_MAX_OFFSET
  R1 (BLASTP, keep):  the best hit's alignment starts >= R1_MIN_UPSTREAM aa upstream of the ATG
                      (qstart <= a - R1_MIN_UPSTREAM) and the subject continues further N-terminal
                      (sstart > 1)
  R3 (Pfam, keep):    a Pfam domain envelope starts >= R3_MIN_UPSTREAM aa upstream of the ATG
                      (env_from <= a - R3_MIN_UPSTREAM)

R2 calls trim-to-ATG; R1 or R3 call keep-5'-open. When both a trim and a keep call fire, or none
fires, the ab initio coding-vs-UTR rule decides. Thresholds: TransDecoderPaper evidence E33 step 3
(fixed before evaluation on simulated Trinity assemblies and the cross-species benchmark).
"""
from typing import List, Optional, Tuple

R2_MAX_SSTART = 3
R2_MAX_OFFSET = 5
R1_MIN_UPSTREAM = 5
R3_MIN_UPSTREAM = 1


def homology_start_call(atg_aa: int, blast: Optional[Tuple[int, int]] = None,
                        pfam_env_from: Optional[int] = None) -> Tuple[List[str], Optional[str]]:
    """
    Args:
        atg_aa: 1-based peptide position of the first in-frame ATG
        blast: (qstart, sstart) of the ORF's best BLASTP hit, or None (no hit / no BLASTP input)
        pfam_env_from: smallest envelope start of the ORF's Pfam domains, or None
    Returns:
        (rules fired, call) with call 'trim', 'keep', 'conflict' or None (no rule fired)
    """
    fired = []
    if blast is not None:
        qstart, sstart = blast
        if sstart <= R2_MAX_SSTART and abs(qstart - (sstart - 1) - atg_aa) <= R2_MAX_OFFSET:
            fired.append("R2")
        if qstart <= atg_aa - R1_MIN_UPSTREAM and sstart > 1:
            fired.append("R1")
    if pfam_env_from is not None and pfam_env_from <= atg_aa - R3_MIN_UPSTREAM:
        fired.append("R3")
    trim = "R2" in fired
    keep = "R1" in fired or "R3" in fired
    if trim and keep:
        return fired, "conflict"
    if trim:
        return fired, "trim"
    if keep:
        return fired, "keep"
    return fired, None
