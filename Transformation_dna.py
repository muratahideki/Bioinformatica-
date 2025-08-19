# molde2dna = sequence
#duplication
def duplication(molde2dna):
    molde2dna.replace("A","T")
    molde2dna.replace("T","A")
    molde2dna.replace("C","G")
    molde2dna.replace("G","C")

    return molde2dna

#transcription
# molde2rna = sequence 
def transcription(molde2rna):
    molde2rna.replace("A","U")
    molde2rna.replace("T","A")
    molde2rna.replace("C","G")
    molde2rna.replace("G","C")
