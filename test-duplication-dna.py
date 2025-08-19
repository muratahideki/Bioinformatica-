from Transformation-dna.py import duplication 

def test_duplication_adenina():
  assert duplication("A") == "T"

def test_duplication_timina():
  assert duplication("T") == "A"

def test_duplication_citocina():
  assert duplication("C") == "G"

def test_duplication_guanina():
  assert duplication("G") ==  "C"
