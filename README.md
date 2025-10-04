# Bioinformatica-

# Objetivos<br>
Estudo de Oryza sativa 

# Métodos<br>

## 0. Pré configurando o ambiente

### 0.1 criando o ambiente virtual 

- Instalar o conda<br>
``` Bash
# Baixar o instalador
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Rodar o instalador
bash Miniconda3-latest-Linux-x86_64.sh
```

- O conda é como se fosse um virtual environment mais poderoso, fazendo a comparação com o virtual environment para python 

| Recurso                  | virtualenv / venv (Python)    | Conda                                                          |
| ------------------------ | ----------------------------- | -------------------------------------------------------------- |
| Linguagem                | Só Python                     | Python, R, C/C++, Java, etc.                                   |
| Pacotes                  | Apenas pacotes Python (pip)   | Qualquer pacote, de qualquer linguagem, com dependências       |
| Isolamento               | Sim, cria um ambiente isolado | Sim, cria um ambiente isolado completo                         |
| Gerenciamento de versões | Só Python e pip packages      | Python, pacotes e até bibliotecas do sistema (ex.: zlib, hdf5) |

- Depois criamos então um virtual environment especíico para trabalhar com bioinformática dentro do próprio conda, e onde podemos instalar os pacotes mais que desejamos

``` bash 
# Criar ambiente com Python e pacotes
conda create -n ngs_env python=3.11 trimmomatic fastqc multiqc

# Ativar o ambiente
conda activate ngs_env

# Instalar um pacote adicional depois
conda install -c bioconda samtools
```

## 1. Levantamento e escolha de dados públicos e pré processamento de dados:<br>
 - Localizar conjuntos de dados RNA-Seq públicos no NCBI GEO ou SRA relacionados a
plantas sob estresse.
 - Justificar a escolha do experimento.

Para realizar isso foi realizado o download do arquivo FASTA no [National Library of Medicine](https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&acc=SRR30542639&display=download) que foi refenciados pelo [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276295)<br>
No caso, foi usado o banco de dados do NCBI SRA de código de acesso [SRP302030](https://www.ncbi.nlm.nih.gov/sra/?term=SRP302030.)

1.1 FATQC
Foi baixado os arquivos Fastq no formato gz. e usado a flag:<br>
`fastqc` `nome do arquivo(pode ser mais de um)` `-o diretório de destino/`<br>

 - Isso vai gerar um arquivo html. 

## Resultados<br>
