import argparse

from bioalgos.alignment import align
from bioalgos.alignment.visualization import format_alignment


def main():
    parser = argparse.ArgumentParser(description="Run pairwise sequence alignment.")

    parser.add_argument("seq1", help="First input sequence")
    parser.add_argument("seq2", help="Second input sequence")

    parser.add_argument(
        "--method",
        default="nw",
        choices=["nw", "sw", "gotoh-global", "gotoh-local", "gotoh-freeshift"],
        help="Alignment method",
    )

    parser.add_argument("--match", type=int, default=1, help="Match score")
    parser.add_argument("--mismatch", type=int, default=-1, help="Mismatch score")
    parser.add_argument("--gap", type=int, default=1, help="Linear gap cost")
    parser.add_argument("--gap-open", type=int, default=-5, help="Affine gap open score")
    parser.add_argument("--gap-extend", type=int, default=-1, help="Affine gap extend score")

    parser.add_argument(
        "--chunk-size",
        type=int,
        default=60,
        help="Line width for alignment visualization",
    )

    args = parser.parse_args()

    result = align(
        args.seq1,
        args.seq2,
        method=args.method,
        match_score=args.match,
        mismatch_score=args.mismatch,
        gap_cost=args.gap,
        gap_open=args.gap_open,
        gap_extend=args.gap_extend,
    )

    # Unified return format
    aligned_seq1, aligned_seq2, score, checked_score, matrices, stats = result

    print(f"Method: {args.method}")
    print(f"Score: {score}")
    print(f"Checked score: {checked_score}")
    print()

    print(format_alignment(aligned_seq1, aligned_seq2, chunk_size=args.chunk_size))
    print()

    print("Stats:")
    for key, value in stats.items():
        print(f"  {key}: {value}")


if __name__ == "__main__":
    main()