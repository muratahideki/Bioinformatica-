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

 Vou detalhar como foi criado o venv: 
 
 ``` bash
 conda create -n ngs_env -c bioconda -c conda-forge trimmomatic fastqc multiqc cutadapt samtools bwa
 ```

## 1. Levantamento e escolha de dados públicos e pré processamento de dados:<br>
 - Localizar conjuntos de dados RNA-Seq públicos no NCBI GEO ou SRA relacionados a
plantas sob estresse.
 - Justificar a escolha do experimento.

Para realizar isso foi realizado o download do arquivo FASTA no [National Library of Medicine](https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&acc=SRR30542639&display=download) que foi refenciados pelo [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276295)<br>
No caso, foi usado o banco de dados do NCBI SRA de código de acesso [SRP302030](https://www.ncbi.nlm.nih.gov/sra/?term=SRP302030.)

### 1.1 FASTQC
Foi baixado os arquivos Fastq no formato gz. e usado a flag:<br>
`fastqc` `nome do arquivo(pode ser mais de um)` `-o diretório de destino/`<br>

 - Isso vai gerar um arquivo html. 

### 1.2 Aplicar Trimmomatic para remoção de adaptadores e filtragem de baixa qualidade.

- Foi usado o Trimmomatic
- As bibliotecas de RNA-Seq foram preparadas com o Illumina mRNA-Seq Sample Preparation v2 kit, e os fragmentos de cDNA foram ligados a adaptadores de sequenciamento da Illumina
- Isso significa que, para fazer o trimming dos dados (por exemplo, no Trimmomatic ou Cutadapt), você deve usar os adaptadores padrão da Illumina TruSeq (presentes no próprio pacote do Trimmomatic, no arquivo TruSeq3-PE.fa ou TruSeq3-SE.fa).

- algoritimo:
```bash
trimmomatic PE [opções] \
  <entrada_R1> <entrada_R2> \
  <saida_paired_R1> <saida_unpaired_R1> \
  <saida_paired_R2> <saida_unpaired_R2> \
  [adaptador]
  [filtros e parâmetros]
```
Adaptador:
- ILLUMINACLIP:<arquivo_fasta>:2:30:10 → remove adaptadores Illumina.
- <arquivo_fasta> = arquivo de adaptadores (TruSeq3-PE.fa para paired-end).
- 2:30:10 = parâmetros de corte (mismatch, palindrome clip, simple clip).

- proposição:
```bash
trimmomatic PE -threads 10 \
  aSRR13447971.fastq.gz  \
  aSRR13447971.fastq.gz aSRR13447971.fastq.gz \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
```
Filtro:
- LEADING:3 → corta bases com qualidade <3 no início da leitura
- TRAILING:3 → corta bases com qualidade <3 no final da leitura
- SLIDINGWINDOW:4:15 → verifica uma janela de 4 bases, corta se a média <15
- MINLEN:36 → descarta leituras menores que 36 bases após o corte



  
#### O que são adaptadores em RNA-Seq (ou qualquer NGS Illumina)?

Quando o cDNA (ou DNA) é preparado para sequenciamento, ele é fragmentado. Mas a máquina Illumina não consegue ler fragmentos de DNA soltos: é preciso dar a eles uma "etiqueta" para que a plataforma reconheça, amplifique e leia cada molécula.<br>
Essa “etiqueta” são os adaptadores.<br>
Eles são pequenas sequências sintéticas de DNA que são ligadas (coladas) nas extremidades dos fragmentos da sua amostra.

#### Funções principais dos adaptadores

- Reconhecimento na flow cell (chip de sequenciamento) → Os adaptadores têm regiões que se ligam quimicamente à superfície da flow cell da Illumina, fixando o DNA/cDNA para que ele seja amplificado e lido.
- Primers para PCR e sequenciamento → Dentro dos adaptadores existem sequências conhecidas que servem de ponto de ancoragem para os primers usados tanto na amplificação da biblioteca quanto no sequenciamento em si.
- Índices (barcodes) → Muitas vezes os adaptadores carregam sequências de “códigos de barras” (index sequences), que permitem misturar várias amostras no mesmo sequenciamento (multiplexing) e depois separar cada uma computacionalmente.

#### Por que remover adaptadores dos reads?

Durante o sequenciamento, se um fragmento de DNA/cDNA for muito curto, a máquina pode acabar lendo além do fragmento e invadindo a região do adaptador.
Isso gera sequências artificiais (não biológicas) no final dos reads.

Se você não remover isso:<br>
- Vai atrapalhar o alinhamento (o read pode não mapear corretamente).
- Pode gerar falsos positivos em análises downstream.

Por isso ferramentas como Trimmomatic ou Cutadapt usam esses arquivos (TruSeq3-PE.fa etc.) que contêm as sequências conhecidas dos adaptadores Illumina, para localizar e cortar esses trechos indesejados.
  
