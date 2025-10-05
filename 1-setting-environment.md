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

 Vou detalhar como foi criado o venv: 
 
 ``` bash
 conda create -n ngs_env -c bioconda -c conda-forge trimmomatic fastqc multiqc cutadapt samtools bwa
 ```
- como o arquivo é paired-end, e os dois sentidos estão misturados no mesmo arquivo preciso instalar sra-tools para usar o fasterqc
- o fasterqc vai separa o arquivo fasta em dois sentidos de leitura: o direto e o reverso 

```bash
conda install -c bioconda sra-tools
conda install -c bioconda hisat2
```
