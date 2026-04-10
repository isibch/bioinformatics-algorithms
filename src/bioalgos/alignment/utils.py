# Function to calculate the score of a pair of characters
def score_pair(a, b, match_score, mismatch_score):
    if a == b:
        return match_score
    return mismatch_score
    
# Function to calculate the score of the aligned sequences
def check_alignment_score(aligned_seq1, aligned_seq2, match_score, mismatch_score, gap_cost):
    total_score = 0

    for a, b in zip(aligned_seq1, aligned_seq2):
        if a == "-" or b == "-":
            total_score -= gap_cost
        else:
            total_score += score_pair(a, b, match_score, mismatch_score)

    return total_score
    
# Function to calculate alignment statistics
def compute_alignment_stats(aligned_seq1, aligned_seq2):
    matches = 0
    mismatches = 0
    gaps = 0

    for a, b in zip(aligned_seq1, aligned_seq2):
        if a == "-" or b == "-":
            gaps += 1
        elif a == b:
            matches += 1
        else:
            mismatches += 1

    alignment_length = len(aligned_seq1)
    identity = matches / alignment_length if alignment_length > 0 else 0

    return {
        "matches": matches,
        "mismatches": mismatches,
        "gaps": gaps,
        "alignment_length": alignment_length,
        "identity": identity,
    }

# Function to print the score matrix
def print_score_matrix(seq1, seq2, matrix):

    header = [" "] + [" "] + list(seq2)
    print("   ".join(header))

    for i, row in enumerate(matrix):

        if i == 0:
            label = " "
        else:
            label = seq1[i - 1]

        row_values = [str(x) for x in row]

        print(label + "  " + "  ".join(row_values))

# Function to check the score of an affine gap alignment
def check_affine_alignment_score(
    aligned_seq1,
    aligned_seq2,
    match_score,
    mismatch_score,
    gap_open,
    gap_extend,
):
    total_score = 0
    gap_in_seq1 = False
    gap_in_seq2 = False

    for a, b in zip(aligned_seq1, aligned_seq2):
        if a == "-":
            if not gap_in_seq1:
                total_score += gap_open
            else:
                total_score += gap_extend

            gap_in_seq1 = True
            gap_in_seq2 = False

        elif b == "-":
            if not gap_in_seq2:
                total_score += gap_open
            else:
                total_score += gap_extend

            gap_in_seq2 = True
            gap_in_seq1 = False

        else:
            total_score += score_pair(a, b, match_score, mismatch_score)
            gap_in_seq1 = False
            gap_in_seq2 = False

    return total_score

# Function to check the score of an affine gap free-shift alignment (ignoring terminal gaps)
def check_affine_freeshift_score(
    aligned_seq1,
    aligned_seq2,
    match_score,
    mismatch_score,
    gap_open,
    gap_extend,
):
    # Ignore leading terminal gaps
    left = 0
    while left < len(aligned_seq1) and (
        aligned_seq1[left] == "-" or aligned_seq2[left] == "-"
    ):
        left += 1

    # Ignore trailing terminal gaps
    right = len(aligned_seq1) - 1
    while right >= left and (
        aligned_seq1[right] == "-" or aligned_seq2[right] == "-"
    ):
        right -= 1

    # If nothing remains, score is 0
    if left > right:
        return 0

    total_score = 0
    gap_in_seq1 = False
    gap_in_seq2 = False

    for a, b in zip(aligned_seq1[left:right + 1], aligned_seq2[left:right + 1]):
        if a == "-":
            if not gap_in_seq1:
                total_score += gap_open
            else:
                total_score += gap_extend

            gap_in_seq1 = True
            gap_in_seq2 = False

        elif b == "-":
            if not gap_in_seq2:
                total_score += gap_open
            else:
                total_score += gap_extend

            gap_in_seq2 = True
            gap_in_seq1 = False

        else:
            total_score += score_pair(a, b, match_score, mismatch_score)
            gap_in_seq1 = False
            gap_in_seq2 = False

    return total_score