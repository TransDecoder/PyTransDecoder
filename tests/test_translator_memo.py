"""Memoised codon translation equals whole-sequence Biopython translation."""
import random
import warnings

from pytransdecoder.core.translator import Translator


def test_memo_equals_biopython():
    rng = random.Random(9)
    warnings.simplefilter("ignore")
    for code in ("Standard", "Tetrahymena", "mitochondrial_vertebrates"):
        t = Translator(code)
        for alphabet in ("ACGT", "ACGTN", "ACGTNRYKMSWBDHV", "acgtn"):
            for _ in range(400):
                s = "".join(rng.choice(alphabet) for _ in range(rng.randint(0, 120)))
                for fr in (0, 1, 2):
                    assert t.translate(s, fr) == t._translate_reference(s, fr), (code, s, fr)


def test_invalid_letters_raise_like_reference():
    t = Translator("Standard")
    for s in ("ATGXXX", "AT-GCC"):
        try:
            ref = t._translate_reference(s)
        except Exception as e:  # noqa: BLE001
            ref = type(e)
        try:
            got = t.translate(s)
        except Exception as e:  # noqa: BLE001
            got = type(e)
        assert got == ref
