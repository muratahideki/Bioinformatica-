# molde_dna = sequence
#duplication

def duplication(molde2dna):
    sequencia = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    fita_duplicada = ""
    for base in molde2dna:
        fita_duplicada += sequencia[base]
    return fita_duplicada

#transcription
# molde2rna = sequence 

def transcription(molde2rna):
    sequencia = {'A':'U', 'T':'A', 'C':'G', 'G':'C' }
    fita_rna = ""
    for base in molde2rna:
        fita_rna += sequencia[base]
    return fita_rna
