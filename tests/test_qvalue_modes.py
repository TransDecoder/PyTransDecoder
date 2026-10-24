import argparse

import pytest

from pytransdecoder import cli
from pytransdecoder.predict import DEFAULT_DECOY_Q, PRECISE_DECOY_Q, SENSITIVE_DECOY_Q


def _cutoff(*argv):
    p = argparse.ArgumentParser()
    cli._add_qvalue_args(p)
    return cli._qvalue_cutoff(p.parse_args(list(argv)))


def test_named_and_numeric_cutoffs():
    assert _cutoff() == DEFAULT_DECOY_Q == 0.05
    assert _cutoff("--precise") == PRECISE_DECOY_Q
    assert _cutoff("--sensitive") == SENSITIVE_DECOY_Q
    assert _cutoff("--max-qvalue", "0.1") == 0.1
    assert _cutoff("--decoy-q", "0.3") == 0.3


def test_modes_exclusive():
    with pytest.raises(SystemExit):
        _cutoff("--precise", "--max-qvalue", "0.1")
