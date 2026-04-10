from bioalgos.alignment.smith_waterman import smith_waterman

def test_sw_local_alignment():
    aligned1, aligned2, score, checked_score, matrices, stats = smith_waterman(
        "TTACGTAA", "GGACGTGG"
    )

    assert aligned1 == "ACGT"
    assert aligned2 == "ACGT"
    assert score == checked_score