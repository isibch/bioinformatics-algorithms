def needleman_wunsch(seq1, seq2, match_score=1, mismatch_score=-1, gap_cost=1):
    def score(a, b):
        if a == b:
            return match_score
        return mismatch_score

    # Function to calculate the score of the aligned sequences
    def check_score(aligned_seq1, aligned_seq2):
        total_score = 0
        
        for a, b in zip(aligned_seq1, aligned_seq2):
            if a == "-" or b == "-":
                total_score -= gap_cost
            else:
                total_score += score(a, b)
            
        return total_score

    # Function to calculate alignment statistics
    def alignment_stats(aligned_seq1, aligned_seq2):
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

        if alignment_length > 0:
            identity = matches / alignment_length
        else:
            identity = 0

        return {
            "matches": matches,
            "mismatches": mismatches,
            "gaps": gaps,
            "alignment_length": alignment_length,
            "identity": identity
        }
        
    # Initialize the scoring matrix
    n = len(seq1)
    m = len(seq2)
    
    score_matrix = [[0] * (m + 1) for _ in range(n + 1)]
    traceback_matrix = [[None] * (m + 1) for _ in range(n + 1)]
    
    traceback_matrix[0][0] = "done"

    # Fill the first row and column with gap penalties
    for i in range(n + 1):
        score_matrix[i][0] = -gap_cost * i
        if i > 0:
            traceback_matrix[i][0] = "up"
    for j in range(m + 1):
        score_matrix[0][j] = -gap_cost * j
        if j > 0:
            traceback_matrix[0][j] = "left"

    # Fill the rest of the scoring matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = score_matrix[i - 1][j - 1] + score(seq1[i - 1], seq2[j - 1])
            delete = score_matrix[i - 1][j] - gap_cost
            insert = score_matrix[i][j - 1] - gap_cost
            
            best_score = max(match, delete, insert)
            score_matrix[i][j] = best_score
            
            if best_score == match:
                traceback_matrix[i][j] = "diagonal"
            elif best_score == delete:
                traceback_matrix[i][j] = "up"
            else:
                traceback_matrix[i][j] = "left"

    # Traceback to find the optimal alignment
    aligned_seq1 = []
    aligned_seq2 = []
    i, j = n, m

    while i > 0 or j > 0:
        current_direction = traceback_matrix[i][j]
        
        if current_direction == "diagonal":
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
            
        elif current_direction == "up":
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append("-")
            i -= 1
            
        else:
            aligned_seq1.append("-")
            aligned_seq2.append(seq2[j - 1])
            j -= 1

    # Add remaining gaps
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
    
    alignment_score = check_score(aligned_seq1, aligned_seq2)
    stats = alignment_stats(aligned_seq1, aligned_seq2)

    return aligned_seq1, aligned_seq2, score_matrix[n][m], score_matrix, alignment_score, stats

# Utility function to print the score matrix in a readable format
def print_score_matrix(seq1, seq2, matrix):

    header = [" "] + [" "] + list(seq2)
    print("   ".join(header))
    
    for i, row in enumerate(matrix):
        
        if i == 0:
            label = " "
        else:
            label = seq1[i-1]
        
        row_values = [str(x) for x in row]
        
        print(label + "  " + "  ".join(row_values))