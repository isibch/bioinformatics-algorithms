from .needleman_wunsch import needleman_wunsch
from .smith_waterman import smith_waterman


def align(seq1, seq2, method="nw", match_score=1, mismatch_score=-1, gap_cost=1):
    if method in {"needleman-wunsch", "nw"}:
        return needleman_wunsch(
            seq1,
            seq2,
            match_score=match_score,
            mismatch_score=mismatch_score,
            gap_cost=gap_cost,
        )

    if method in {"smith-waterman", "sw"}:
        return smith_waterman(
            seq1,
            seq2,
            match_score=match_score,
            mismatch_score=mismatch_score,
            gap_cost=gap_cost,
        )

    raise ValueError(
        "Unknown method. Supported methods are: 'nw', 'needleman-wunsch', 'sw', 'smith-waterman'."
    )