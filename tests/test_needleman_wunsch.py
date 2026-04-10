from bioalgos.alignment.needleman_wunsch import needleman_wunsch

def test_nw_simple_match():
    aligned1, aligned2, score, checked_score, matrices, stats = needleman_wunsch("A", "A")

    assert aligned1 == "A"
    assert aligned2 == "A"
    assert score == checked_score
    assert stats["matches"] == 1