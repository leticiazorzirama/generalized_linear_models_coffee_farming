# GRUPO 2 — Manejo Fitossanitário e Qualidade na Cafeicultura

**Base de dados:** `dados/grupo2_cafeicultura.csv` — 340 talhões · 19 variáveis
**Área:** Agronomia / Fitopatologia

---

## 1. Contexto do problema

Uma cooperativa de cafeicultores de Minas Gerais reúne 340 talhões monitorados
ao longo de uma safra. A cooperativa enfrenta três problemas simultâneos que
consomem seu orçamento técnico e definem sua margem:

- a **broca-do-café** (*Hypothenemus hampei*), monitorada por armadilhas
  instaladas em cada talhão;
- a **ferrugem alaranjada** (*Hemileia vastatrix*), avaliada pela severidade —
  a fração da área foliar lesionada;
- a **classificação do lote** para o mercado de exportação, que define o preço
  recebido.

O corpo técnico precisa decidir onde concentrar o esforço de manejo na próxima
safra e quais recomendações fazer aos cooperados quanto a cultivar, densidade de
plantio, adubação e irrigação. As três respostas têm naturezas estatísticas
completamente distintas e exigem, cada uma, uma família diferente de MLG.

---

## 2. Dicionário de variáveis

> Observação: os níveis dos fatores estão gravados **sem acentuação** no CSV
> (`Nao`, `Organico`, `Catuai`, `Sul de Minas`, etc.).

### Identificação e delineamento

| Variável | Tipo | Descrição |
|---|---|---|
| `id_talhao` | texto | Identificador do talhão |
| `regiao_produtora` | categórica (4) | Sul de Minas, Cerrado Mineiro, Mogiana, Matas de Minas |
| `cultivar` | categórica (4) | Catuai, Mundo Novo, Bourbon, Icatu |
| `manejo` | categórica (3) | Convencional, Integrado, Organico |
| `irrigacao` | categórica (2) | `Sim` / `Nao` |

### Covariáveis candidatas

| Variável | Tipo | Descrição |
|---|---|---|
| `altitude_m` | numérica | Altitude média do talhão (m) |
| `declividade_pct` | numérica | Declividade média (%) |
| `idade_lavoura_anos` | numérica | Idade da lavoura (anos) |
| `densidade_plantio` | numérica | Plantas por hectare |
| `adubacao_n_kg_ha` | numérica | Adubação nitrogenada (kg N/ha) |
| `precipitacao_safra_mm` | numérica | Precipitação acumulada na safra (mm) |
| `umidade_relativa_pct` | numérica | Umidade relativa média (%) |
| `ph_solo` | numérica | pH do solo |
| `materia_organica_pct` | numérica | Matéria orgânica do solo (%) |

### Esforço de amostragem (monitoramento da broca)

| Variável | Tipo | Descrição |
|---|---|---|
| `n_armadilhas` | contagem | Armadilhas instaladas no talhão (4 a 12) |
| `dias_exposicao` | contagem | Dias de exposição das armadilhas (15, 20 ou 30) |

### Variáveis-resposta

| Variável | Tipo | Descrição |
|---|---|---|
| `n_brocas_capturadas` | contagem | Adultos de broca capturados |
| `severidade_ferrugem` | proporção contínua (0,1) | Fração da área foliar lesionada |
| `padrao_exportacao` | binária | Lote atingiu o padrão de exportação (1) ou não (0) |

---

## 3. Os três desafios

### Desafio A — Modelo **Poisson** (ou extensão para contagens)

> **Pergunta do corpo técnico:** *Quais condições de talhão favorecem a
> infestação pela broca, e quanto se ganha em pressão de praga ao alterar cada
> uma delas?*

**O problema oculto do esforço amostral:** o número de brocas capturadas depende
diretamente de quantas armadilhas foram instaladas e por quantos dias ficaram
expostas. Um talhão com 12 armadilhas por 30 dias tem **seis vezes** mais
esforço de captura que um com 4 armadilhas por 15 dias. Comparar contagens
brutas entre talhões é, portanto, comparar coisas diferentes.

**Tarefas:**
1. Construam a medida de esforço amostral adequada a partir de `n_armadilhas` e
   `dias_exposicao` e incorporem-na ao modelo pelo mecanismo correto do MLG.
   Justifiquem por que esse termo entra com coeficiente **fixado em 1** e não
   como preditor livre. *(Verifiquem essa afirmação empiricamente: ajustem
   também o modelo com o log do esforço como preditor livre e comparem o
   coeficiente estimado com 1.)*
2. Verifiquem a equidispersão. Reportem a razão de Pearson sobre os graus de
   liberdade e conduzam um teste formal de superdispersão.
3. Caso haja superdispersão, comparem quasi-Poisson e binomial negativa,
   explicando **como cada uma modela a variância** e por que isso importa aqui.
4. Selecionem variáveis e interpretem os efeitos como razões de taxas de captura
   por unidade de esforço.

---

### Desafio B — Modelo **Beta**

> **Pergunta do corpo técnico:** *Que fatores agronômicos determinam a
> severidade da ferrugem, e existe uma cultivar comprovadamente mais resistente?*

`severidade_ferrugem` é uma **proporção contínua** no intervalo aberto (0, 1).
Não existe denominador: não se contaram "folhas lesionadas de um total de N",
mediu-se diretamente a fração de área foliar comprometida por análise de
imagem. Isso **não é** um problema binomial.

A distribuição é assimétrica, limitada nos dois extremos e sua variância depende
da média — necessariamente menor perto de 0 e de 1 do que perto de 0,5. Um
modelo linear pode predizer severidade negativa ou acima de 100%.

**Tarefas:**
1. Explicitem por que este caso **não** admite modelo binomial nem regressão
   linear sobre a proporção, e por que a transformação arco-seno-raiz
   (historicamente usada em fitopatologia) é uma solução inferior à regressão
   beta. Comparem as duas abordagens empiricamente.
2. Ajustem a regressão beta (`betareg::betareg`), justificando a ligação
   escolhida para a média.
3. **Modelo de dispersão variável:** investiguem se o parâmetro de precisão φ é
   constante. Ajustem um modelo com preditor também para φ
   (`y ~ x1 + x2 | z1`) e comparem por teste da razão de verossimilhanças.
4. Interpretem os efeitos. Como a ligação `logit` é não linear, o efeito de uma
   variável em pontos percentuais de severidade **depende do nível das demais** —
   apresentem efeitos marginais em cenários agronômicos concretos.
5. Diagnóstico: resíduos quantílicos padronizados, meio-normal com envelope
   simulado e identificação de talhões influentes.

> **Nota sobre os dados:** a severidade foi registrada com piso de 0,001
> (limite de detecção do método). Discutam o que fariam se houvesse zeros
> exatos na base — a regressão beta padrão não os admite.

---

### Desafio C — Modelo **Binomial**

> **Pergunta do corpo técnico:** *O que determina a chance de um lote atingir o
> padrão de exportação?*

`padrao_exportacao` é binária. Regressão linear sobre ela produziria
probabilidades preditas fora de [0, 1] e resíduos estruturalmente
heterocedásticos.

**Tarefas:**
1. Ajustem o modelo logístico e conduzam seleção de variáveis.
2. **Questão de raciocínio causal:** `severidade_ferrugem` é ela própria um
   desfecho (Desafio B), influenciada por cultivar, irrigação e densidade.
   Ajustem o modelo **com** e **sem** essa variável e comparem os coeficientes
   de `cultivar`. Discutam: ao incluir a severidade, o coeficiente da cultivar
   passa a medir o efeito *direto* (não mediado pela ferrugem), enquanto sem ela
   mede o efeito *total*. Qual dos dois responde melhor à pergunta do produtor
   que precisa **escolher a cultivar** antes do plantio?
3. Avaliem a capacidade preditiva: matriz de confusão, curva ROC e AUC,
   discutindo a escolha do ponto de corte em função do custo assimétrico de
   classificar erradamente um lote.
4. Interpretem em razões de chances com intervalos de confiança.

---

## 4. Requisitos comuns aos três modelos

**(a) Especificação** — Família e ligação justificadas pela natureza da resposta.

**(b) Seleção de variáveis** — Modelo amplo reduzido com critério explícito
(AIC/BIC, razão de verossimilhanças, análise de deviance). Documentem o caminho.
Há variáveis irrelevantes por construção na base; identificá-las é parte da
tarefa. Verifiquem multicolinearidade entre as variáveis de solo e clima.

**(c) Diagnóstico de resíduos** — Para cada modelo final:
- resíduos quantílicos aleatorizados contra o preditor linear;
- envelope simulado / meio-normal;
- verificação de dispersão;
- pontos influentes identificados pelo `id_talhao`;
- linearidade das covariáveis contínuas na escala do preditor linear.

**(d) Decisão** — Recomendação de manejo para a próxima safra, com efeitos
quantificados na escala original (número de brocas, pontos percentuais de
severidade, probabilidade de exportação).

---

## 5. Entrega

- Relatório de até 8 páginas · script reprodutível · apresentação de 12 minutos.

## 6. Critérios de avaliação

| Critério | Peso |
|---|---|
| Justificativa da família e ligação a partir da natureza do dado | 20% |
| Uso correto do *offset* e distinção binomial × beta | 15% |
| Processo de seleção de variáveis documentado e criticado | 20% |
| Diagnóstico de resíduos e tratamento da dispersão | 25% |
| Interpretação na escala original e recomendação final | 20% |
