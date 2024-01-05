

def reverse_dna(seq):
    return seq[::-1]

def anti_sense_dna(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in seq)

def dimer_look(r1, r2, n_base):
    n_base = int(n_base)

    z1 = ""
    z2 = ""
    z3 = ""
    j = 0
    n = 0
    x = 0
    y = 0
    z = 0
    r = 0
    f = -1
    w = 0
    l = 0
    l1 = len(r1)
    l2 = len(r2)

    if n_base < 1:
        n_base = 1
    if n_base > 10:
        n_base = 10

    init_k = n_base + 1
    str_len = n_base + 2
    filter_size = n_base + 5

    if l1 < l2:
        z1 = r2
        z3 = reverse_dna(r1)
        l = l1
        l1 = l2
        l2 = l
    else:
        z1 = r1
        z3 = reverse_dna(r2)

    z2 = anti_sense_dna(z3)
    px = {}
    r1x = z1
    r2y = z2

    for y in range(l2 - init_k):
        x = -1
        while True:
            x = z1.find(z2[y:y + init_k], x + 1)
            if x == -1:
                break
            w = l2 + x - y
            z = y + l1 - x - 1
            j = 0
            l = init_k
            r = init_k

            if px.get(w, 0) == 0:
                while x + l <= l1 - 1:
                    if y + l > l2 - 1:
                        break
                    if r1x[x + l] == r2y[y + l]:
                        l += 1
                        r += 1
                    else:
                        if x + l + 1 > l1 - 1:
                            break
                        if y + l + 1 > l2 - 1:
                            break
                        if r1x[x + l + 1] != r2y[y + l + 1]:
                            break
                        l += 2
                        r += 1

                while x - j - 1 >= 0:
                    if y - j - 1 < 0:
                        break
                    if r1x[x - j - 1] == r2y[y - j - 1]:
                        l += 1
                        r += 1
                        j += 1
                    else:
                        if x - j - 2 < 0:
                            break
                        if y - j - 2 < 0:
                            break
                        if r1x[x - j - 2] != r2y[y - j - 2]:
                            break
                        l += 2
                        r += 1
                        j += 2

                if r > str_len:
                    f += 1
                    px[w] = 1

    return f > -1

"""
change the value n_base -> increasing it makes dimer search more permissive, decreasing more restrictive
"""
def has_dimer(primer1, primer2, n_base=6):
    return dimer_look(primer1, primer1, n_base) or dimer_look(primer1, primer2, n_base) or dimer_look(primer2, primer2, n_base)

