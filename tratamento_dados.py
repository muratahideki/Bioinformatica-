#criar um lista, separando em uma lista de 3 em 3
list_codon =  [sequencia[i:i+3] for i in range(0,len(sequencia)-len(sequencia)%3,3)]

codon_start = 'ATG'
codon_stop = 'TAA', 'TAG', 'TGA'


def pos_codon(list_codon, codon):
    codon_position = []
    for i in range(len(list_codon)):
        if list_codon[i] in codon:
            codon_position.append(i)
    return codon_position

pos_start = pos_codon(list_codon,codon_start)
pos_stop = pos_codon(list_codon, codon_stop)

#criar os pares de [start,stop]

pars = []
for s in pos_start:
    stops__after = [stop for stop in pos_stop if stop > s]
    if stops__after:
        pars.append([s, stops__after[0]])

#dimensão dessa seq de codons 
dimension = []
for p in range(len(pars)):
    mag = pars[p][1]-pars[p][0] + 1 
    dimension.append(mag)

#sequencia desejada 
def sequencia(start,stop):
    seq = list_codon[start,stop]
    return seq 

