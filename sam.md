
-  o arquivo .sam é legível
-  o arquivo .bam é binário, mas é menos pesado
  
```bash
samtools view -S -b saida.sam > arquivo.bam
samtools sort arquivo.bam -o arquivo_sorted.bam
samtools index arquivo_sorted.bam
```
