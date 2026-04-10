from bioalgos.alignment.gotoh_global import gotoh_global
from bioalgos.alignment.gotoh_local import gotoh_local
from bioalgos.alignment.gotoh_freeshift import gotoh_freeshift

def test_gotoh_global():
    aligned1, aligned2, score, checked_score, matrices, stats = gotoh_global(
        "AAAAAA", "AAA"
    )

    assert score == checked_score
    assert stats["gaps"] == 3

def test_gotoh_local():
    aligned1, aligned2, score, checked_score, matrices, stats = gotoh_local(
        "TTACGTAA", "GGACGTGG"
    )

    assert aligned1 == "ACGT"
    assert aligned2 == "ACGT"

def test_gotoh_freeshift():
    aligned1, aligned2, score, checked_score, matrices, stats = gotoh_freeshift(
        "TTACGTAA", "ACGT"
    )

    assert score == checked_score
    assert stats["matches"] == 4