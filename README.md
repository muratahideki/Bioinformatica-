# Bioinformatica-

## Objetivos<br>
Estudo e prática da matéria bioinformatica

## Automação de testes unitários<br>
Foi criado triggers em um arquivo `yml`, para que seja ativado quando ocorre um `push` ou  `pull-request` na branch `testes-unitarios`<br>
No `jobs` foi criado uma `build` que rodar na última versão do Ubuntu. A próxima etapa é clonar o repositório para uma máquina virtual (o runner), onde ocorre todos os steps, por isso o `job` inicia com `uses: actions/checkout@v4` 
