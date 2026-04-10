import time
from statistics import mean

from bioalgos.alignment import (
    needleman_wunsch,
    smith_waterman,
    gotoh_global,
    gotoh_local,
    gotoh_freeshift,
)


METHODS = {
    "nw": lambda s1, s2: needleman_wunsch(
        s1,
        s2,
        match_score=1,
        mismatch_score=-1,
        gap_cost=1,
    ),
    "sw": lambda s1, s2: smith_waterman(
        s1,
        s2,
        match_score=1,
        mismatch_score=-1,
        gap_cost=1,
    ),
    "gotoh-global": lambda s1, s2: gotoh_global(
        s1,
        s2,
        match_score=1,
        mismatch_score=-1,
        gap_open=-5,
        gap_extend=-1,
    ),
    "gotoh-local": lambda s1, s2: gotoh_local(
        s1,
        s2,
        match_score=1,
        mismatch_score=-1,
        gap_open=-5,
        gap_extend=-1,
    ),
    "gotoh-freeshift": lambda s1, s2: gotoh_freeshift(
        s1,
        s2,
        match_score=1,
        mismatch_score=-1,
        gap_open=-5,
        gap_extend=-1,
    ),
}


TEST_CASES = [
    ("GATTACA", "GCATGCU"),
    ("TTACGTAA", "GGACGTGG"),
    ("TTACGTAA", "ACGT"),
    ("ACCCCCGT", "ACGT"),
    ("AAAAAA", "AAA"),
]


def benchmark_method(func, seq1, seq2, repeats=20):
    runtimes = []
    score = None
    checked_score = None

    for _ in range(repeats):
        start = time.perf_counter()
        result = func(seq1, seq2)
        end = time.perf_counter()

        runtimes.append(end - start)

        _, _, score, checked_score, _, _ = result

    return score, checked_score, mean(runtimes)


def main():
    print("Benchmarking alignment methods\n")

    for seq1, seq2 in TEST_CASES:
        print(f"Sequences: {seq1} vs {seq2}")
        print(f"{'method':<18} {'score':<8} {'checked':<8} {'avg_time_ms':<12}")
        print("-" * 52)

        for method_name, func in METHODS.items():
            score, checked_score, avg_time = benchmark_method(func, seq1, seq2)
            print(f"{method_name:<18} {score:<8} {checked_score:<8} {avg_time * 1000:<12.4f}")

        print()


if __name__ == "__main__":
    main()