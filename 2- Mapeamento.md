## 2- Mapeamento 

Faz mais sentido escolher o método com Histat2, pelos seguintes motivos:

| Característica           | HISAT2                         | STAR                                  |
| ------------------------ | ------------------------------ | ------------------------------------- |
| Velocidade               | Razoável                       | Muito rápido                          |
| Uso de RAM               | Baixo (4–8 GB)                 | Alto (30–60 GB)                       |
| Sensibilidade a splicing | Boa                            | Excelente                             |
| Detecção de isoformas    | Limitada                       | Muito boa                             |
| Leituras longas          | Bom                            | Excelente                             |
| Quando usar              | PC limitado / análises simples | Projetos grandes / alta sensibilidade |

### Preparação<br>

- precisa criar arquivos Fasta, a partir dos arquivos FASTQ
```bash
seqtk seq -a SRR13447971_1.fastq > SRR13447971_1.fasta
seqtk seq -a SRR13447971_2.fastq > SRR13447971_2.fasta
``` 
- precisa indexar o arquivo Fasta
```bash 
hisat2-build SRR13447971_1.fasta genome_index
hisat2-build SRR13447971_2.fasta genome_index
```

### Alinhar com HISAT2

```bash
hisat2 -x SRR1344797_1_trimmed.fastq.gz -1 SRR1344797_2_trimmed.fastq.gz -2 sample_R2.fastq -S output.sam -p
```
- -S output.sam -p: arquivo de saída 
