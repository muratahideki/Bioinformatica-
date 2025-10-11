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

- precisa indexar o arquivo Fasta do genoma de referência 
```bash 
hisat2-build SofficinarumxspontaneumR570_771_v2.0.softmasked.fa sugarcane_index
```
- Isso vai gerar alguns arquivos de referência:

```
sugarcane_index.1.ht2
sugarcane_index.2.ht2
...
sugarcane_index.8.ht2
```
- (Opcional) Usar a anotação no alinhamento

Se quiser usar locais conhecidos de splicing, você pode gerar um arquivo de splice sites a partir do GFF3:

```
hisat2_extract_splice_sites.py SofficinarumxspontaneumR570_771_v2.1.gene.gff3 > splicesites.txt
```

### Alinhar com HISAT2

```bash
hisat2 -x sugarcane_index \
  --known-splicesite-infile splicesites.txt \
  -1 SRR13447971_1.fastq.gz -2 SRR13447971_2.fastq.gz \
  -S output.sam -p 4
```
- -S output.sam -p: arquivo de saída 
