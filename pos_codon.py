def pos_codon(list_codon, codon):
    codon_position = []
    for i in range(len(list_codon)):
        if list_codon[i] in codon:
            codon_position.append(i)
    return codon_position

