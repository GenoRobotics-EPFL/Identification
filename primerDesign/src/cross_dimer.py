

def reverse_dna(seq):
    return seq[::-1]

def anti_sense_dna(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in seq)

def extend_match_length(seq1, seq2, start, end, step):
    match_length = 0
    while 0 <= start < len(seq1) and 0 <= end < len(seq2) and seq1[start] == seq2[end]:
        match_length += 1
        start += step
        end += step
    return match_length

def dimer_look(seq1, seq2, n_base):
    n_base = max(1, min(n_base, 10))
    init_k = n_base + 1
    str_len = n_base + 2

    if len(seq1) < len(seq2):
        seq1, seq2 = seq2, reverse_dna(seq1)
    else:
        seq2 = reverse_dna(seq2)

    anti_sense_seq2 = anti_sense_dna(seq2)
    found_matches = {}

    for start in range(len(seq2) - init_k):
        pos = -1
        while True:
            pos = seq1.find(anti_sense_seq2[start:start + init_k], pos + 1)
            if pos == -1:
                break
            match_pos = len(seq2) + pos - start
            if found_matches.get(match_pos, 0) == 0:
                match_length = init_k
                match_length += extend_match_length(seq1, anti_sense_seq2, pos + init_k, start + init_k, 1)
                match_length += extend_match_length(seq1, anti_sense_seq2, pos - 1, start - 1, -1)

                if match_length > str_len:
                    found_matches[match_pos] = 1

    return any(found_matches.values())

"""
change the value n_base -> increasing it makes dimer search more permissive, decreasing more restrictive
"""
def has_dimer(primer1, primer2, n_base=6):
    return dimer_look(primer1, primer1, n_base) or dimer_look(primer1, primer2, n_base) or dimer_look(primer2, primer2, n_base)

