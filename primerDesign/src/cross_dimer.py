

def reverse_dna(seq):
    return seq[::-1]

def anti_sense_dna(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in seq)

def dimer_look(r1, r2, n_base):
    n_base, z1, z2, z3, f = int(n_base), "", "", "", -1

    if n_base < 1: n_base = 1
    if n_base > 10: n_base = 10

    init_k, str_len, filter_size = n_base + 1, n_base + 2, n_base + 5

    if len(r1) < len(r2): z1, z3, l = r2, reverse_dna(r1), len(r1)
    else: z1, z3, l = r1, reverse_dna(r2), len(r2)

    z2, px, r1x, r2y = anti_sense_dna(z3), {}, z1, anti_sense_dna(z3)

    for y in range(len(r2) - init_k):
        x, w, z, j, l, r = -1, 0, 0, 0, init_k, init_k

        while True:
            x = z1.find(z2[y:y + init_k], x + 1)
            if x == -1: break
            print(z2)
            print(x)
            print(y)
            w, z = z2 + x - y, y + z1 - x - 1
            j = 0

            if px.get(w, 0) == 0:
                while x + l <= l1 - 1:
                    if y + l > l2 - 1: break
                    if r1x[x + l] == r2y[y + l]: l, r = l + 1, r + 1
                    else:
                        if x + l + 1 > l1 - 1 or y + l + 1 > l2 - 1 or r1x[x + l + 1] != r2y[y + l + 1]: break
                        l, r = l + 2, r + 1

                while x - j - 1 >= 0:
                    if y - j - 1 < 0: break
                    if r1x[x - j - 1] == r2y[y - j - 1]: l, r, j = l + 1, r + 1, j + 1
                    else:
                        if x - j - 2 < 0 or y - j - 2 < 0 or r1x[x - j - 2] != r2y[y - j - 2]: break
                        l, r, j = l + 2, r + 1, j + 2

                if r > str_len:
                    f += 1
                    px[w] = 1

    return f > -1

def has_dimer(primer1, primer2, n_base):
    return dimer_look(primer1, primer1, n_base) or dimer_look(primer1, primer2, n_base) or dimer_look(primer2, primer2, n_base)

sequences = ["CTGTTGAAGTTTCATCTACAAATGGATAATACT", "CATACAATCGCTATTCATAATGGAAAGGAA", "GAATTCAGTTAACGACGAGATTTAGTATCC", 
"GACATAGGATGCCACTCTTTAAAATGAAAA", "CAAGAAGGTTGGTATTGCTCCTTTATTTTT", "CAACCTCAAAGAAAAGCTTTCTCTTTTATC"]

self_or_cross_binding = False
n_base = 4
for i in range(len(sequences)):
    self_or_cross_binding |= has_dimer(sequences[i], sequences[i], n_base)
    for j in range(i+1, len(sequences)):
        self_or_cross_binding |= has_dimer(sequences[i], sequences[j], n_base)

if (self_or_cross_binding):
    print("The primers provided have a self or cross binding")
else:
    print("The primers provided are free of any self or cross binding")

