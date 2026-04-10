def format_alignment(aligned_seq1, aligned_seq2, chunk_size=60):
    match_line = []

    for a, b in zip(aligned_seq1, aligned_seq2):
        if a == b and a != "-":
            match_line.append("|")
        elif a == "-" or b == "-":
            match_line.append(" ")
        else:
            match_line.append(".")

    match_line = "".join(match_line)

    blocks = []
    for start in range(0, len(aligned_seq1), chunk_size):
        end = start + chunk_size
        blocks.append(
            "\n".join(
                [
                    aligned_seq1[start:end],
                    match_line[start:end],
                    aligned_seq2[start:end],
                ]
            )
        )

    return "\n\n".join(blocks)