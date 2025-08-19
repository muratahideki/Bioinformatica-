# Bioinformatica-

## Objetivos<br>
Estudo e prática da matéria bioinformatica

## Automação de testes unitários<br>
Foi criado triggers em um arquivo `yml`, para que seja ativado quando ocorre um `push` ou  `pull-request` na branch `testes-unitarios`<br>
No `jobs` foi criado uma `build` que rodar na última versão do Ubuntu. A próxima etapa é clonar o repositório para uma máquina virtual (o runner), onde ocorre todos os steps, por isso o `job` inicia com `uses: actions/checkout@v4` 

# Como preparar o arquivo de teste unitário<br>
1º Arrange/ Preparar
- Definir os dados que serão testados
- valores de entrada 
- vaolor esperado
  
2º Act/ Agir
- chama a função que vai ser testada

3º Assert/ Verificar
- Vertificar se o valor retornado pela função é o mesmo que o valor esperado
