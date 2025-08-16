def duplication(molde2dna):
    molde2dna = molde2dna.replace("A", "t")
    molde2dna = molde2dna.replace("T", "a")
    molde2dna = molde2dna.replace("C", "g")
    molde2dna = molde2dna.replace("G", "c")
    # Troca para maiúsculas para o resultado final
    return molde2dna.upper()

def transcription(molde2rna):
    molde2rna = molde2rna.replace("A", "u")
    molde2rna = molde2rna.replace("T", "a")
    molde2rna = molde2rna.replace("C", "g")
    molde2rna = molde2rna.replace("G", "c")
    # Troca para maiúsculas para o resultado final
    return molde2rna.upper()
