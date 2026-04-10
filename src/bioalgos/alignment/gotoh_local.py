from .utils import score_pair, compute_alignment_stats, check_affine_alignment_score


def gotoh_local(seq1, seq2, match_score=1, mismatch_score=-1, gap_open=-5, gap_extend=-1):
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

    # Variables for tracking the best local alignment endpoint
    max_score = 0
    max_position = (0, 0)
    max_state = "M"

    # Initialize first row and column
    # For local alignment, the alignment may start anywhere,
    # so the first row and column allow a restart with score 0.
    for i in range(n + 1):
        tb_M[i][0] = "done"

    for j in range(m + 1):
        tb_M[0][j] = "done"

    # Fill the DP matrices
    for i in range(1, n + 1):
        for j in range(1, m + 1):

            # Score for aligning the two current characters
            pair_score = score_pair(seq1[i - 1], seq2[j - 1], match_score, mismatch_score)

            # Match/Mismatch state
            # Alignment continues with a character pair
            from_M = M[i - 1][j - 1] + pair_score
            from_Ix = Ix[i - 1][j - 1] + pair_score
            from_Iy = Iy[i - 1][j - 1] + pair_score

            # Local alignment allows restarting the alignment at score 0
            M[i][j] = max(0, from_M, from_Ix, from_Iy)

            # Store where we came from for traceback
            if M[i][j] == 0:
                tb_M[i][j] = "done"
            elif M[i][j] == from_M:
                tb_M[i][j] = "M"
            elif M[i][j] == from_Ix:
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

            # Track the best scoring cell across all matrices
            if M[i][j] > max_score:
                max_score = M[i][j]
                max_position = (i, j)
                max_state = "M"

            if Ix[i][j] > max_score:
                max_score = Ix[i][j]
                max_position = (i, j)
                max_state = "Ix"

            if Iy[i][j] > max_score:
                max_score = Iy[i][j]
                max_position = (i, j)
                max_state = "Iy"

    # Traceback starts from the best local endpoint
    aligned_seq1 = []
    aligned_seq2 = []

    i, j = max_position
    state = max_state

    while i > 0 or j > 0:

        if state == "M":
            # Stop when the alignment reaches score 0
            if M[i][j] == 0 or tb_M[i][j] == "done":
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

    # Reverse because traceback builds alignment backwards
    aligned_seq1.reverse()
    aligned_seq2.reverse()

    aligned_seq1 = "".join(aligned_seq1)
    aligned_seq2 = "".join(aligned_seq2)

    # Recompute the alignment score independently to verify that traceback is consistent
    alignment_score = check_affine_alignment_score(
        aligned_seq1,
        aligned_seq2,
        match_score,
        mismatch_score,
        gap_open,
        gap_extend
    )

    stats = compute_alignment_stats(aligned_seq1, aligned_seq2)

    return (
        aligned_seq1,
        aligned_seq2,
        max_score,
        alignment_score,
        {
            "M": M,
            "Ix": Ix,
            "Iy": Iy,
            "tb_M": tb_M,
            "tb_Ix": tb_Ix,
            "tb_Iy": tb_Iy,
        },
        stats,
    )