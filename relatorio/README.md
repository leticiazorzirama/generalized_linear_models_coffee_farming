# Relatório — PBL Cafeicultura (Grupo 2)

Relatório em LaTeX, padrão ABNT, limite de **8 páginas** (enunciado, seção 5).

**Você não precisa saber LaTeX para escrever neste relatório.** Comece pelo guia rápido
abaixo: são três comandos e um único arquivo para editar.

---

## Guia rápido

### 1. Abrir um terminal na pasta do projeto

| sistema | como fazer |
| --- | --- |
| **Windows** | abra a pasta `glm-cafeicultura` no Explorador de Arquivos, clique com o botão direito num espaço vazio e escolha **"Abrir no Terminal"** |
| **macOS** | Finder → botão direito na pasta → **Serviços** → **Novo Terminal na Pasta** |
| **Linux** | botão direito na pasta → **Abrir no Terminal** |

Se usa RStudio, é mais fácil: menu **Tools → Terminal → New Terminal** já abre no lugar
certo.

### 2. Preparar o computador — só na primeira vez

Copie e cole a linha do seu sistema:

**Windows**

```
.\relatorio\setup.ps1
```

**macOS ou Linux**

```
./relatorio/setup.sh
```

Demora alguns minutos e baixa cerca de 250 MB. Ele instala só o que faltar e não mexe no
que já existe. Quando terminar, imprime **AMBIENTE PRONTO**.

> **Windows reclamou de "execução de scripts está desabilitada"?**
> Rode esta linha uma vez e tente de novo:
> ```
> Set-ExecutionPolicy -Scope CurrentUser -ExecutionPolicy RemoteSigned
> ```

### 3. Gerar o PDF

```
Rscript relatorio/build.R
```

O arquivo sai em **`relatorio/main.pdf`**.

### 4. Escrever

Abra **`relatorio/conteudo.tex`** em qualquer editor (RStudio, VS Code, Bloco de Notas).
**É o único arquivo que você precisa editar.**

Escreva → salve → rode o comando do passo 3 → abra o PDF.

---

## Os avisos: como saber o que ainda falta

Existem dois tipos de aviso no projeto. Entender os dois é como a gente controla o que
está pronto e o que não está.

### Aviso 1 — as marcas vermelhas dentro do PDF

Abrindo o `main.pdf`, você vê trechos assim, em vermelho e negrito:

> **[PENDENTE: Seção em consolidação: o modelo final da equipe ainda está sendo revisado.]**

São lembretes que **nós mesmos** escrevemos no texto, marcando o que ainda não está
pronto. Aparecem no PDF de propósito, para ninguém esquecer.

Para escrever uma nova, dentro do `conteudo.tex`:

```latex
\pendente{falta a tabela de efeitos marginais}
```

**Todas têm que sumir antes da entrega.**

### As 8 pendências de hoje

| seção do relatório | o que falta |
| --- | --- |
| Desafio A — abertura | o modelo final da equipe ainda está em revisão |
| Desafio A — *offset* | o intervalo de confiança do γ não contém 1 (ver abaixo) |
| Desafio A — diagnóstico | resíduos e figura do envelope simulado |
| Desafio C — abertura | o modelo logístico ainda não foi ajustado |
| Desafio C — predição | matriz de confusão, curva ROC, AUC e ponto de corte |
| Decisão | consolidar depois que A e C fecharem |
| Considerações finais | redigir |

> **Sobre o γ do Desafio A:** deixando o coeficiente do esforço amostral livre em vez de
> fixá-lo em 1, o valor estimado é **0,898**, com intervalo **[0,816; 0,980]** — que
> **não contém 1**. É a hipótese que justifica usar o *offset*, e sob esse conjunto de
> covariáveis ela é rejeitada. Com os fatores categóricos no modelo, o intervalo passa a
> conter 1. Não é fatal, é questão de qual especificação reportar — mas é o item que vale
> **15% da nota**, e a pergunta vai vir na apresentação.

### Aviso 2 — o conferidor

```
Rscript relatorio/conferir.R
```

Ele revisa o relatório antes de você gastar uma compilação. A saída vem em blocos
numerados:

| bloco | o que significa | quebra o PDF? |
| --- | --- | --- |
| **1. MACROS** | você usou um apelido de número que não existe | **sim** |
| **3. FIGURAS** | citou uma imagem que não está na pasta | **sim** |
| **4. CITAÇÕES** | citou um artigo que não está no `refs.bib` | **sim** — sai `[?]` |
| **5. NÚMEROS DIGITADOS** | há número escrito à mão no texto | não, mas envelhece |
| **6. PENDÊNCIAS** | quantas marcas vermelhas restam | não |

No final ele resume: **SEM ERROS BLOQUEANTES** ou a lista do que precisa ser corrigido.

O bloco 5 é o mais traiçoeiro e é o assunto da próxima seção.

---

## Por que não digitamos números no texto

Os números do relatório saem dos modelos, não dos nossos dedos.

Ao escrever no `conteudo.tex`, em vez do número use o **apelido** dele:

```latex
% assim NÃO:
O parâmetro de precisão estimado foi 20,20.

% assim SIM:
O parâmetro de precisão estimado foi \bPhi.
```

Na hora de gerar o PDF, o `\bPhi` vira `20,20` sozinho.

**Por que isso importa:** o Desafio A ainda vai mudar e o C nem começou. Se um modelo
mudar e o número estiver digitado, o texto fica errado e **ninguém percebe** — até a
banca perceber. Com o apelido, o texto se corrige sozinho na próxima compilação.

### Como descobrir o apelido de um número

```
Rscript relatorio/numeros.R
```

Ele lista os **87 apelidos** disponíveis com o valor atual de cada um. Os nomes seguem um
padrão:

| começa com | é de |
| --- | --- |
| `\base…` | a base de dados (nº de talhões, esforço amostral) |
| `\a…` | Desafio A — broca |
| `\b…` | Desafio B — ferrugem |
| `\c…` | Desafio C — exportação |

Exemplos: `\baseN` (340 talhões), `\aTheta` (o θ da binomial negativa), `\bPhi` (a
precisão da regressão beta), `\cPrevalencia` (47,06%).

### Preciso de um número que não tem apelido

Avise quem estiver cuidando do `numeros.R` — são poucas linhas para criar um novo. Ou
digite o número e aceite que ele vai aparecer no bloco 5 do conferidor.

### Não é obrigatório

Se digitar `20,20` direto, o PDF sai igual. O apelido é uma garantia disponível, não uma
regra imposta. Mas todo número digitado é um número que alguém vai ter que conferir na
mão depois.

---

## Como escrever no `conteudo.tex`

Não precisa saber LaTeX. Estes são praticamente todos os comandos usados no relatório:

```latex
\section{Título de seção}          % título numerado
\subsection{Subtítulo}

Parágrafo normal, escrito direto.
Uma linha em branco separa parágrafos.

\textbf{negrito}    \emph{itálico}    \var{nome_da_variavel}

\cite{ferrari2004}        % vira (FERRARI; CRIBARI-NETO, 2004)
\textcite{ferrari2004}    % vira Ferrari e Cribari-Neto (2004)

\pendente{o que ainda falta aqui}
```

> **Cuidado com estes caracteres:** `%`, `_`, `&`, `#`, `$` têm significado especial no
> LaTeX. Para escrevê-los literalmente, ponha uma barra invertida antes: `\%`, `\_`,
> `\&`. Essa é a causa mais comum de erro de compilação para quem está começando.

### Colocar uma figura

```latex
\begin{figure}[H]
    \centering
    \includegraphics[width=0.9\textwidth]{nome_do_arquivo.png}
    \caption{Legenda explicando o que a figura mostra.}
    \label{fig:apelido}
\end{figure}
```

A imagem precisa estar em `challenge_a/figuras/`, `challenge_b/figuras/` ou
`challenge_c/figuras/`. O `build.R` copia todas sozinho — use só o nome do arquivo, sem
caminho.

---

## Como o relatório é montado

```
relatorio/
├── conteudo.tex  O TEXTO                    <- e aqui que voce escreve
├── main.tex      capa, margens, pacotes     <- raramente muda
├── numeros.tex   GERADO - nao editar        <- reescrito a cada compilacao
├── numeros.R     o gerador dos numeros      <- so mexe quem trocar de modelo
├── refs.bib      a bibliografia
├── figs/         GERADO - copia das figuras dos desafios
├── setup.R       prepara o computador (os tres sistemas)
├── setup.sh      atalho macOS/Linux - acha ou instala o R
├── setup.ps1     atalho Windows - idem
├── conferir.R    a verificacao antes de entregar
├── build.R       gera o PDF (os tres sistemas)
├── build.sh      atalho macOS/Linux
└── build.ps1     atalho Windows
```

O `build.R` faz três coisas, nesta ordem: recalcula os números a partir dos modelos,
copia as figuras dos desafios, e gera o PDF.

> **Dois itens são gerados automaticamente e não devem ser editados:** o `numeros.tex` e
> a pasta `figs/`. Qualquer alteração neles se perde na próxima compilação.

---

## Trocou o modelo final de algum desafio?

Abra o `numeros.R` e edite **apenas** o bloco marcado `MODELOS FINAIS`, no começo do
arquivo. São três fórmulas, uma por desafio, com comentário dizendo qual é qual.

**Nada nos scripts `challenge_a.R`, `challenge_b.R` ou `challenge_c.R` precisa mudar.** O
`numeros.R` é independente: ele lê o CSV e reajusta os modelos por conta própria.

---

## Se algo der errado

**`Rscript não é reconhecido como comando`**
O R não está instalado, ou está instalado mas fora do PATH do sistema. Rode o
`setup.ps1` (Windows) ou `setup.sh` (macOS/Linux): eles procuram o R em todos os lugares
comuns e, se não acharem, oferecem instalar.

**`! LaTeX Error: File 'alguma-coisa.sty' not found`**
Falta um pacote do LaTeX. Rode `Rscript relatorio/setup.R` de novo — ele instala o que
faltar.

**`! Undefined control sequence` seguido de `\algumaCoisa`**
Você usou um apelido de número que não existe, ou errou a digitação. Rode
`Rscript relatorio/conferir.R`: o bloco 1 diz exatamente qual.

**`Package biblatex Error: Style 'abnt' not found`**
```
Rscript -e "tinytex::tlmgr_install(c('biblatex','biblatex-abnt','biber'))"
```

**A bibliografia sai vazia, ou aparecem `[?]` no lugar das citações**
Faltou a passada do `biber`. Apague os arquivos temporários e compile de novo — no
Windows, pelo Explorador: entre em `relatorio/` e apague tudo que se chama `main.` com
extensão diferente de `.tex` e `.pdf` (`main.aux`, `main.bcf`, `main.bbl`, etc.).

**O PDF não muda depois que eu edito**
Confira se salvou o `conteudo.tex`. Se salvou, apague os temporários como acima.

**Passou de 8 páginas**
O limite é do enunciado. Corte texto — **não** mexa nas margens (`geometry`) nem no
espaçamento (`setspace`): os dois são do padrão ABNT, e alterá-los é descumprir a
formatação exigida.

**Windows: `permission denied` no meio da compilação**
Se a pasta está no OneDrive ou Dropbox, a sincronização pode travar um arquivo. Pause a
sincronização e tente de novo.

---

## Notas técnicas

Para quem quiser entender as escolhas:

- **Classe `article` com miolo ABNT**, e não a classe `abntex2` completa. A `abntex2` é
  feita para monografia e traz folha de rosto, sumário e listas obrigatórias — num limite
  de 8 páginas isso consumiria metade do espaço.
- **Citação por `biblatex-abnt`, e não `abntex2cite`.** O `abntex2cite` está sem
  manutenção e quebra com o kernel LaTeX atual (colisão de `\abntnextkey` com o `expl3`).
  Foi o primeiro erro que apareceu ao montar isto. O `biblatex-abnt` é mantido e produz o
  mesmo padrão. O backend é o `biber`.
- **TinyTeX em vez de MiKTeX ou TeX Live.** Mesmo comando nos três sistemas, não pede
  senha de administrador, e instala pacotes LaTeX faltantes sozinho durante a compilação.
- **`build.R` escrito em R, e não em shell.** Todo mundo do grupo já tem R, e o mesmo
  arquivo roda igual nos três sistemas. O `build.sh` e o `build.ps1` só chamam ele.
- Se você já tem MiKTeX, TeX Live ou MacTeX instalado, não precisa do TinyTeX — o
  `build.R` detecta `latexmk` ou `pdflatex` no PATH e usa o que encontrar.
