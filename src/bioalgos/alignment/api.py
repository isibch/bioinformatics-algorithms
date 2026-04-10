from .needleman_wunsch import needleman_wunsch
from .smith_waterman import smith_waterman
from .gotoh_global import gotoh_global
from .gotoh_local import gotoh_local
from .gotoh_freeshift import gotoh_freeshift


def align(
    seq1,
    seq2,
    method="nw",
    match_score=1,
    mismatch_score=-1,
    gap_cost=1,
    gap_open=-5,
    gap_extend=-1,
):
    if method in {"nw", "needleman-wunsch"}:
        return needleman_wunsch(
            seq1,
            seq2,
            match_score=match_score,
            mismatch_score=mismatch_score,
            gap_cost=gap_cost,
        )

    if method in {"sw", "smith-waterman"}:
        return smith_waterman(
            seq1,
            seq2,
            match_score=match_score,
            mismatch_score=mismatch_score,
            gap_cost=gap_cost,
        )

    if method in {"gotoh-global", "gg"}:
        return gotoh_global(
            seq1,
            seq2,
            match_score=match_score,
            mismatch_score=mismatch_score,
            gap_open=gap_open,
            gap_extend=gap_extend,
        )

    if method in {"gotoh-local", "gl"}:
        return gotoh_local(
            seq1,
            seq2,
            match_score=match_score,
            mismatch_score=mismatch_score,
            gap_open=gap_open,
            gap_extend=gap_extend,
        )

    if method in {"gotoh-freeshift", "gf"}:
        return gotoh_freeshift(
            seq1,
            seq2,
            match_score=match_score,
            mismatch_score=mismatch_score,
            gap_open=gap_open,
            gap_extend=gap_extend,
        )

    raise ValueError(
        "Unknown method. Use one of: nw, sw, gotoh-global, gotoh-local, gotoh-freeshift."
    )