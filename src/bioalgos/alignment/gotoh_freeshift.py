from .utils import score_pair, compute_alignment_stats, check_affine_freeshift_score


def gotoh_freeshift(seq1, seq2, match_score=1, mismatch_score=-1, gap_open=-5, gap_extend=-1):
    n = len(seq1)
    m = len(seq2)

    # Negative infinity is used to disable invalid transitions
    neg_inf = float("-inf")

    # Three DP matrices are used in Gotoh:
    #
    # M  -> alignment ends with a match/mismatch
    # Ix -> alignment ends with a gap in seq2 (vertical move)
    # Iy -> alignment ends with a gap in seq1 (horizontal move)
    #
    # Using separate matrices allows us to distinguish between opening a gap and extending an existing gap.
    M = [[0] * (m + 1) for _ in range(n + 1)]
    Ix = [[neg_inf] * (m + 1) for _ in range(n + 1)]
    Iy = [[neg_inf] * (m + 1) for _ in range(n + 1)]

    # Traceback matrices for each state
    tb_M = [[None] * (m + 1) for _ in range(n + 1)]
    tb_Ix = [[None] * (m + 1) for _ in range(n + 1)]
    tb_Iy = [[None] * (m + 1) for _ in range(n + 1)]

    # First row and column are initialized with 0 because leading gaps are free in freeshift alignment.
    for i in range(n + 1):
        M[i][0] = 0
        tb_M[i][0] = "done"

    for j in range(m + 1):
        M[0][j] = 0
        tb_M[0][j] = "done"

    # Fill the DP matrices
    for i in range(1, n + 1):
        for j in range(1, m + 1):

            # Score for aligning the two current characters
            pair_score = score_pair(seq1[i - 1], seq2[j - 1], match_score, mismatch_score)

            # Match/Mismatch state
            # Alignment continues with a character pair
            prev_M = M[i - 1][j - 1]
            prev_Ix = Ix[i - 1][j - 1]
            prev_Iy = Iy[i - 1][j - 1]

            best_prev = max(prev_M, prev_Ix, prev_Iy)
            M[i][j] = best_prev + pair_score

            # Store where we came from for traceback
            if best_prev == prev_M:
                tb_M[i][j] = "M"
            elif best_prev == prev_Ix:
                tb_M[i][j] = "Ix"
            else:
                tb_M[i][j] = "Iy"

            # Gap in seq2 (vertical gap)
            # Either open a new gap or extend an existing one
            open_gap_x = M[i - 1][j] + gap_open
            extend_gap_x = Ix[i - 1][j] + gap_extend

            Ix[i][j] = max(open_gap_x, extend_gap_x)

            if Ix[i][j] == open_gap_x:
                tb_Ix[i][j] = "M"
            else:
                tb_Ix[i][j] = "Ix"

            # Gap in seq1 (horizontal gap)
            open_gap_y = M[i][j - 1] + gap_open
            extend_gap_y = Iy[i][j - 1] + gap_extend

            Iy[i][j] = max(open_gap_y, extend_gap_y)

            if Iy[i][j] == open_gap_y:
                tb_Iy[i][j] = "M"
            else:
                tb_Iy[i][j] = "Iy"

    # In freeshift alignment, trailing gaps are also free.
    # Therefore the best endpoint can lie anywhere in the last
    # row or the last column.
    final_score = neg_inf
    end_i, end_j, state = n, m, "M"

    for j in range(m + 1):
        if M[n][j] > final_score:
            final_score = M[n][j]
            end_i, end_j, state = n, j, "M"
        if Ix[n][j] > final_score:
            final_score = Ix[n][j]
            end_i, end_j, state = n, j, "Ix"
        if Iy[n][j] > final_score:
            final_score = Iy[n][j]
            end_i, end_j, state = n, j, "Iy"

    for i in range(n + 1):
        if M[i][m] > final_score:
            final_score = M[i][m]
            end_i, end_j, state = i, m, "M"
        if Ix[i][m] > final_score:
            final_score = Ix[i][m]
            end_i, end_j, state = i, m, "Ix"
        if Iy[i][m] > final_score:
            final_score = Iy[i][m]
            end_i, end_j, state = i, m, "Iy"

    # Traceback starts from the best endpoint on the border
    aligned_seq1 = []
    aligned_seq2 = []

    i, j = end_i, end_j

    # Add trailing free gaps after the chosen endpoint
    while i < n:
        aligned_seq1.append(seq1[i])
        aligned_seq2.append("-")
        i += 1

    while j < m:
        aligned_seq1.append("-")
        aligned_seq2.append(seq2[j])
        j += 1

    i, j = end_i, end_j

    while i > 0 or j > 0:

        if state == "M":
            if tb_M[i][j] == "done":
                break

            prev_state = tb_M[i][j]

            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])

            i -= 1
            j -= 1
            state = prev_state

        elif state == "Ix":
            prev_state = tb_Ix[i][j]

            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append("-")

            i -= 1
            state = prev_state

        elif state == "Iy":
            prev_state = tb_Iy[i][j]

            aligned_seq1.append("-")
            aligned_seq2.append(seq2[j - 1])

            j -= 1
            state = prev_state

        else:
            break

    # Add leading free gaps before the aligned region
    while i > 0:
        aligned_seq1.append(seq1[i - 1])
        aligned_seq2.append("-")
        i -= 1

    while j > 0:
        aligned_seq1.append("-")
        aligned_seq2.append(seq2[j - 1])
        j -= 1

    # Reverse because traceback builds alignment backwards
    aligned_seq1.reverse()
    aligned_seq2.reverse()

    aligned_seq1 = "".join(aligned_seq1)
    aligned_seq2 = "".join(aligned_seq2)

    # Recompute the alignment score independently to verify that traceback is consistent
    alignment_score = check_affine_freeshift_score(
        aligned_seq1,
        aligned_seq2,
        match_score,
        mismatch_score,
        gap_open,
        gap_extend
    )

    stats = compute_alignment_stats(aligned_seq1, aligned_seq2)

    return aligned_seq1, aligned_seq2, final_score, alignment_score, M, Ix, Iy, stats