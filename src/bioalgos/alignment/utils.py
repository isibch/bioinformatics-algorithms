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