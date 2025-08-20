# molde_dna = sequence
#duplication

def duplication(molde_dna):
    sequencia = {A:T, T:A, C:G, G:C}
    fita_duplicada = ""
    for base in fita_DNA:
        fita_duplicada += sequencia[base]
    return fita_duplicada

#transcription
# molde2rna = sequence 
def transcription(molde2rna):
    molde2rna.replace("A","U")
    molde2rna.replace("T","A")
    molde2rna.replace("C","G")
    molde2rna.replace("G","C")
