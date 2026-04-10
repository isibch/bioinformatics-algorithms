from .utils import score_pair, check_alignment_score, compute_alignment_stats

def needleman_wunsch(seq1, seq2, match_score=1, mismatch_score=-1, gap_cost=1):
    # Initialize the scoring and traceback matrices
    n = len(seq1)
    m = len(seq2)

    # DP matrix storing alignment scores
    score_matrix = [[0] * (m + 1) for _ in range(n + 1)]
    # Traceback matrix storing the direction of the optimal move
    traceback_matrix = [[None] * (m + 1) for _ in range(n + 1)]

    # Initialize first row and column with gap penalties
    for i in range(1, n + 1):
        score_matrix[i][0] = -gap_cost * i
        traceback_matrix[i][0] = "up"

    for j in range(1, m + 1):
        score_matrix[0][j] = -gap_cost * j
        traceback_matrix[0][j] = "left"

    traceback_matrix[0][0] = "done"

    # Fill the DP matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):

            diagonal_score = score_matrix[i - 1][j - 1] + score_pair(
                seq1[i - 1], seq2[j - 1], match_score, mismatch_score
            )

            up_score = score_matrix[i - 1][j] - gap_cost
            left_score = score_matrix[i][j - 1] - gap_cost

            best_score = max(diagonal_score, up_score, left_score)
            score_matrix[i][j] = best_score

            if best_score == diagonal_score:
                traceback_matrix[i][j] = "diagonal"
            elif best_score == up_score:
                traceback_matrix[i][j] = "up"
            else:
                traceback_matrix[i][j] = "left"

    # Traceback from bottom-right to reconstruct global alignment
    aligned_seq1 = []
    aligned_seq2 = []

    i, j = n, m

    while i > 0 or j > 0:

        direction = traceback_matrix[i][j]

        if direction == "diagonal":
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1

        elif direction == "up":
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append("-")
            i -= 1

        else:
            aligned_seq1.append("-")
            aligned_seq2.append(seq2[j - 1])
            j -= 1

    # Reverse because traceback builds alignment backwards
    aligned_seq1.reverse()
    aligned_seq2.reverse()

    aligned_seq1 = "".join(aligned_seq1)
    aligned_seq2 = "".join(aligned_seq2)

    alignment_score = check_alignment_score(
        aligned_seq1, aligned_seq2, match_score, mismatch_score, gap_cost
    )

    stats = compute_alignment_stats(aligned_seq1, aligned_seq2)

    return (
        aligned_seq1,
        aligned_seq2,
        score_matrix[n][m],
        alignment_score,
        {
            "score_matrix": score_matrix,
            "traceback_matrix": traceback_matrix,
        },
        stats,
    )