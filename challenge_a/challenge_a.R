# =======================================================================================
# UNIVERSIDADE DO VALE DO ITAJAI - UNIVALI
# ESCOLA POLITECNICA
# PROGRAMA DE POS-GRADUACAO EM COMPUTACAO APLICADA - PPGCA
# MESTRADO EM COMPUTACAO APLICADA
# Disciplina: Modelagem estatistica
# Prof. Dr.: Rodrigo Sant'Ana
# Discentes: Andre Lucas Ribeiro, Leticia Zorzi Rama, Matheus Neis
# Itajai, Santa Catarina, Brasil

# =======================================================================================
# DESAFIO A
# =======================================================================================
# Descricao:
# Modelo Poisson (ou extensao para contagens)
# Quais condições de talhao favorecem a infestacao pela broca, 
# e quanto se ganha em pressao de praga ao alterar cada uma delas?
#
# Variavel resposta:
# `n_brocas_capturadas`
# 
# Variaveis candidatas:
# `altitude_m`, `declividade_pct`, `idade_lavoura_anos`, `densidade_plantio`, `adubacao_n_kg_ha`
# `precipitacao_safra_mm`, `umidade_relativa_pct`, `ph_solo`, `materia_organica_pct`
#
# O problema oculto do esforco amostral: 
# o numero de brocas capturadas depende diretamente de quantas armadilhas foram 
# instaladas e por quantos dias ficaram expostas. Um talhao com 12 armadilhas por 
# 30 dias tem seis vezes mais esforco de captura que um com 4 armadilhas por 15 dias. 
# Comparar contagens brutas entre talhoes e, portanto, comparar coisas diferentes.
#
# 1. Construam a medida de esforço amostral adequada a partir de `n_armadilhas` e
#    `dias_exposicao` e incorporem-na ao modelo pelo mecanismo correto do MLG.
#    Justifiquem por que esse termo entra com coeficiente fixado em 1 e nao
#    como preditor livre. (Verifiquem essa afirmação empiricamente: ajustem
#    também o modelo com o log do esforço como preditor livre e comparem o
#    coeficiente estimado com 1.)
# 2. Verifiquem a equidispersão. Reportem a razão de Pearson sobre os graus de
#    liberdade e conduzam um teste formal de superdispersão.
# 3. Caso haja superdispersão, comparem quasi-Poisson e binomial negativa,
#    explicando como cada uma modela a variância e por que isso importa aqui.
# 4. Selecionem variáveis e interpretem os efeitos como razões de taxas de captura
#    por unidade de esforço.
# 
# O script esta estruturado em:
# - Configurar ambiente de trabalho
# - Carregar e manipular a base de dados
# - Analise exploratoria multivariada da variavel resposta com cada uma das variaveis candidatas
# - Analise exploratoria global da variavel resposta com todas as variaveis candidatas
# - Modelo Poisson: 
#   - Modelos lineares generalizados nulos e completos sem e com medidas de esforco amostral, com e sem offset
#   - Ajuste de modelos lineares generalizados exaurindo as combinacoes entre as variaveis preditoras com offset de esforco amostral 
#   - Diagnostico e analise dos residuos 
#   - Verificacao da equidispersao
# - Modelo Binomial negativa
#   - Modelos lineares generalizados nulos e completo
#   - Ajuste de modelos lineares generalizados exaurindo as combinacoes entre as variaveis preditoras
#   - Diagnostico e analise dos residuos 
# - Parecer final (resposta as perguntas do desafio)
#
# # Base de dados: cafeicultura.csv
# =======================================================================================

# =======================================================================================
# CONFIGURAR AMBIENTE DE TRABALHO
# =======================================================================================

# Conferir e configurar caminhos, se necessario
# getwd()
# setwd("inserir caminho")
# getwd()

# Pacotes
# Primeiro pacote base - pacman - verifica se os pacotes
# listados ja estão instalados, caso contrario, ira os instalar
if(!require(pacman)) {
    print(paste0("Pacote ainda não instalado. Instalando..."))
    install.packages("pacman", dependencies = TRUE)
    library(pacman)
} else {
    print(paste0("Pacote já instalado e carregado"))
    library(pacman)
}

# Demais pacotes
p_load(ggplot2, dplyr, MuMIn, ggplot, ggcorrplot, patchwork, readxl, AER,
       DataExplorer, ggpubr, scatterplot3d, effects, car, hnp, statmod,
       datasets, DCluster, stargazer, olsrr, performance, report)

# Opcoes gerais
options(scipen = 10)

# Funcao auxiliar de impressao de texto
secao <- function(texto) {
barra <- paste(rep("=", 150), collapse = "")
  cat("\n", barra, "\n", texto, "\n", barra, "\n\n", sep = "")
}

# =======================================================================================
# BASE DE DADOS
# =======================================================================================
secao("CARREGAR BASE DE DADOS")

# Carregar a base de dados
caf <- read.csv("data/cafeicultura.csv", sep = ",")

# Selecionar as variaveis de interesse (resposta, candidatas e esforco amostral)
caf1 <- caf %>%
  dplyr::select(
    n_brocas_capturadas,
    adubacao_n_kg_ha,
    altitude_m, 
    declividade_pct,
    densidade_plantio,
    idade_lavoura_anos,
    materia_organica_pct,
    ph_solo,
    precipitacao_safra_mm,
    umidade_relativa_pct
  )

# =======================================================================================
# ANALISE EXPLORATORIA MULTIVARIADA 
# Variavel resposta com cada variavel candidata
# =======================================================================================
secao("ANÁLISE EXPLORATÓRIA MULTIVARIADA")

# Verificar as relacoes da variavel resposta com as variaveis candidatas
secao("Análise multivariada - relações da variável resposta com as variáveis candidatas")
plot_scatterplot(data = caf1, by = "n_brocas_capturadas",
                 geom_point_args = list(pch = 21, fill = "white",
                                        size = 2))

# Funcao para calcular as correlacoes de Pearson, Kendall e Spearman 
correlacao_calc <- function(var, name) {
  
  pearson <- cor(caf1$n_brocas_capturadas, caf1[[var]])
  kendall <- cor(caf1$n_brocas_capturadas, caf1[[var]], method = "kendall")
  spearman <- cor(caf1$n_brocas_capturadas, caf1[[var]], method = "spearman")
  
  secao(paste("Análise multivariada - relação de n_brocas_capturadas com", var))
  
  cat(name, "\nPearson:", pearson, ">", correlacao_aval(pearson))
  cat("\nKendall:", kendall, ">", correlacao_aval(kendall))
  cat("\nSpearman:", spearman, ">", correlacao_aval(spearman), "\n")

}

# Funcao para avaliar a correlacao
correlacao_aval <- function(r) {
  if(is.na(r)) return("não disponível")
  if(r == 0) return("nula")
  
  intensidade <- case_when(
    abs(r) < 0.30 ~ "fraca",
    abs(r) < 0.70 ~ "moderada",
    TRUE ~ "forte"
  )
  
  direcao <- ifelse(r > 0, "positiva", "negativa")
  
  paste(intensidade, direcao)
}

# Aplicando as funcoes de correlacao de 'n_brocas_capturadas' com cada variavel candidata
correlacao_calc("adubacao_n_kg_ha", "Adubação nitrogenada (kg/ha)")
correlacao_calc("altitude_m", "Altidude (m)")
correlacao_calc("declividade_pct", "Declividade (%)")
correlacao_calc("densidade_plantio", "Densidade (plantas/ha)")
correlacao_calc("idade_lavoura_anos", "Idade (anos)")
correlacao_calc("materia_organica_pct", "Matéria orgânica (%)")
correlacao_calc("precipitacao_safra_mm", "Precipitação (mm)")
correlacao_calc("ph_solo", "pH")
correlacao_calc("umidade_relativa_pct", "Umidade (%)")

# Parecer
cat("\nOs resultados das correlações indicam uma associação moderada positiva do número 
    de brocas capturadas com a umidade. As demais variáveis possuem associação fraca 
    com a variável resposta.")

# =======================================================================================
# ANALISE EXPLORATORIA GLOBAL 
# Variavel resposta com todas as variaveis candidatas
# =======================================================================================
secao("ANÁLISE EXPLORATÓRIA GLOBAL")

# Matriz de correlacao
mat <- cor(caf1)

# Arredondar
round(mat, digits = 3)

# Plotar a matriz
mat_plot <- ggcorrplot(mat, method = "square", lab = TRUE, type = "lower")
print(mat_plot)

# Calcular teste de hipotese para verificar quais associacoes sao significativas
p_mat <- cor_pmat(caf1)

# Marcar as relacoes nao significativas
mat_plot <- ggcorrplot(mat, method = "square", lab = TRUE, p.mat = p_mat,
                  type = "lower")
print(mat_plot)

# Manter apenas as relacoes significativas
mat_plot <- ggcorrplot(mat, method = "square", lab = TRUE, p.mat = p_mat,
                  insig = "blank", type = "lower")
print(mat_plot)

# Parecer
cat("\nA associação moderada e positiva do número de brocas capturadas com a umidade relativa
    se demonstrou significativa a partir do teste de hipóteses. Das demais associações avaliadas 
    como fracas, aquelas significativas são:\n
    - precipitação, matéria orgânica e densidade de plantio com direção positiva, e\n
    - adubação e altitude com direção negativa.\n
    A matriz de correlação também revelou uma associação moderada entre umidade e precipitação.")

# =======================================================================================
# MODELOS LINEARES GENERALIZADOS - POISSON
# Modelos:
# - mod00: nulo
# - mod01: completo SEM medida de esforco amostral
# - mod02: modelo nulo COM medida de esforco amostral como OFFSET
# - mod03: completo COM medida de esforco amostral como OFFSET 
# - mdo04: completo COM log da medida de esforco amostral como preditor livre
# - mod05: modelo com minimo AICc
# - mod06: modelo final
# =======================================================================================
secao("MODELOS LINEARES GENERALIZADOS - POISSON")

# n_brocas_capturadas e uma variavel de contagem o que significa dizer que e uma variavel discreta
# Poisson para modelar variaveis discretas

# --------------------------------------
# mod00: modelo poisson nulo
# --------------------------------------
# Para este primeiro caso e para fins de aprendizagem, 
# o significado de cada argumento sera comentado conforme constante em ?glm
#?glm

secao("Modelo nulo")
mod00 <- glm(
  formula = n_brocas_capturadas ~ 1, # descricao simbolica do modelo, a funcao a ser modelada
  family = poisson(link = "log"), # descricao da distribuicao do erro e da funcao de ligacao 
  data = caf1, # base de dados com as variaveis de interesse
  weights = NULL, # pesos para as variaveis preditoras
  # subset, # subset de observacoes a ser utilizadas no ajuste do modelo
  na.action = "na.fail", # retorna o objeto se nao ha dados ausentes ou sinaliza um erro se ha dados ausentes
  # start, # valores de partida para os parametros do preditor linear
  # etastart, # valores de partida para o preditor (qual a diferenca deste com o start?)
  offset = NULL, # termo a ser adicionado no preditor linear com um coeficiente estimado em 1
  # control, lista de parametros para controlar o ajuste
  # model, retorna o frame do modelo
  # method, # metodo para ser usado no ajuste, .fit e o padrao e faz o ajuste, caso seja .frame nao se faz o ajuste
  # x, y, indica se as variaveis resposta e preditoras devem ser incluidas nos valores retornados
  # singular.ok, (verificar o que significa)
  # constrasts, (verificar o que significa)
  # intercept, indica se o intercepto deve ser incluido no modelo nulo
  # object, objeto derivado da classe glm
  # type, tipo dos pesos a serem extraidos do modelo ajustado
)

# Verificando se intercepto do modelo nulo coincide com o log da media da variavel resposta
mod00$coefficients
log(mean(caf1$n_brocas_capturadas))

# --------------------------------------
# mod01: modelo poisson completo
# sem esforco amostral
# --------------------------------------
secao("Modelo completo SEM medida de esforço amostral")
mod01 <- glm(
  n_brocas_capturadas ~ adubacao_n_kg_ha
    + altitude_m
    + declividade_pct
    + densidade_plantio
    + idade_lavoura_anos
    + materia_organica_pct
    + ph_solo
    + precipitacao_safra_mm
    + umidade_relativa_pct, 
  family = poisson(link = "log"),
  data = caf1,
  na.action = "na.fail"
)

# Teste de verossimilhanca entre os modelos nulo e completo SEM esforco amostral
secao("Teste de verossimilhança entre os modelos nulo e completo SEM medida de esforço amostral")

# Calcular o teste qui-quadrado passo a passo 
# Obter o log da verossimilhanca de cada modelo
ll_null <- logLik(mod00)
ll_full <- logLik(mod01)

# Obter o dobro da diferenca
chi_sq_stat <- as.numeric(2 * (ll_full - ll_null))

# Obter a diferenca dos graus de liberdade
df_diff <- attr(ll_full, "df") - attr(ll_null, "df")

# Calcular o p-valor
p_value <- pchisq(chi_sq_stat, df = df_diff, lower.tail = FALSE)

# Imprimir os resultados
cat("Chi-Square Stat:", chi_sq_stat, "\nDF:", df_diff, "\np-value:", p_value, "\n")

# Verificar que os mesmos resultados sao retornados com anova utilizando o teste qui-quadrado
anova(mod00, mod01, test = "Chisq")

# Reducao percentual da deviancia nula
residuos_mod01 <- ((deviance(mod00) - deviance(mod01)) / deviance(mod00)) * 100

# Parecer
cat("\nEm comparação ao modelo nulo, o modelo completo SEM esforço amostral reduziu a deviância nula em ",
    residuos_mod01, "%, indicando que o modelo explicou esse percentual da deviância nula")

# Analise da deviancia das covariaveis do modelo completo SEM esforco amostral
secao("Analise da deviância das covariáveis do modelo completo SEM medida de esforço amostral")
anova(mod01, test = "Chisq")

# Parecer 
cat("\nPara o modelo completo SEM esforço amostral, em termos de redução de resíduos de deviância nula, as variáveis com 
    redução significativa foram todas, exceto pH solo.")

# Tabela de estimacao dos parametros do modelo completo SEM esforco amostral
secao("Estimação dos parâmetros do modelo completo SEM medida de esforço amostral")
summary(mod01)

# Parecer 
cat("\nPara o modelo completo SEM esforço amostral, em termos de efeito, as variáveis com efeito 
    significativo foram quase todas, exceto declividade, pH solo e precipitação.")

# --------------------------------------
# Construcao do esforco amostral
# --------------------------------------

# Relembrando o problema
# 1. Construam a medida de esforço amostral adequada a partir de `n_armadilhas` e
#    `dias_exposicao` e incorporem-na ao modelo pelo mecanismo correto do MLG.
#    Justifiquem por que esse termo entra com coeficiente fixado em 1 e não
#    como preditor livre. (Verifiquem essa afirmação empiricamente: ajustem
#    também o modelo com o log do esforço como preditor livre e comparem o
#    coeficiente estimado com 1.)

# Construir medida de esforco amostral
caf <- caf %>%
  mutate(
    esforco_amostral = n_armadilhas * dias_exposicao
  )

# Adicionar a medida de esforco amostral na base de dados que esta sendo manipulada
caf1$esforco_amostral <- caf$esforco_amostral

# --------------------------------------
# mod02: modelo nulo com esforco
# amostral passado como offset
# --------------------------------------
secao("Modelo nulo COM medida de esforço amostral passada como OFFSET")
mod02 <- glm(
  n_brocas_capturadas ~ 1, 
  family = poisson(link = "log"),
  data = caf1,
  na.action = "na.fail",
  offset = log(esforco_amostral)
)

# --------------------------------------
# mod03: modelo completo com esforco 
# amostral passado como offset
# --------------------------------------
secao("Modelo completo COM medida de esforço amostral passada como OFFSET")
mod03 <- glm(
  n_brocas_capturadas ~ adubacao_n_kg_ha
    + altitude_m
    + declividade_pct
    + densidade_plantio
    + idade_lavoura_anos
    + materia_organica_pct
    + ph_solo
    + precipitacao_safra_mm
    + umidade_relativa_pct, 
  family = poisson(link = "log"),
  data = caf1,
  na.action = "na.fail",
  offset = log(esforco_amostral)
)

# Teste de verossimilhanca entre os modelos nulo e completo COM esforco amostral passado como OFFSET
secao("Teste de verossimilhança entre os modelos nulo e completo COM medida de esforço amostral passada como OFFSET")
anova(mod02, mod03, test = "Chisq")

# Reducao percentual da deviancia nula
residuos_mod03 <- ((deviance(mod02) - deviance(mod03)) / deviance(mod02)) * 100

# Parecer
cat("\nEm comparação ao modelo nulo, o modelo completo COM esforço amostral passado como OFFSET reduziu 
    a deviância nula em ", residuos_mod03, "%, indicando que o modelo explicou esse percentual da deviância nula. 
    A redução da deviância entre os modelos nulo e completo SEM offset de esforço amostral foi de", residuos_mod01, "%. 
    Portanto, a incorporação do offset no mecanismo do modelo linear generalizado resultou em uma diferença de", 
    residuos_mod03 - residuos_mod01, "% na redução da deviância nula e no ganho de explicação da mesma.")

# Analise da deviancia das covariaveis do modelo completo COM esforco amostral passado como OFFSET
secao("Analise da deviância das covariáveis do modelo completo COM medida de esforço amostral")
anova(mod03, test = "Chisq")

# Parecer 
cat("\nPara o modelo completo COM esforço amostral passado como OFFSET, em termos de redução de resíduos de deviância nula, 
    as variáveis com redução significativa foram todas, exceto declividade e ph solo.")

# Tabela de estimacao dos parametros do modelo completo COM esforco amostral passado como OFFSET
secao("Estimação dos parâmetros do modelo completo COM medida de esforço amostral passada como OFFSET")
summary(mod03)

# Parecer 
cat("\nPara o modelo completo COM esforço amostral passado como OFFSET, em termos de efeito, as variáveis com 
    efeito significativo foram todas, exceto declividade, ph solo e precipitação.")

# --------------------------------------
# mod04: modelo completo com log do esforco 
# amostral passado como preditor livre
# --------------------------------------
secao("Modelo completo COM medida de esforço amostral passada como PREDITOR LIVRE")
mod04 <- glm(
  n_brocas_capturadas ~ adubacao_n_kg_ha
    + altitude_m
    + declividade_pct
    + densidade_plantio
    + idade_lavoura_anos
    + materia_organica_pct
    + ph_solo
    + precipitacao_safra_mm
    + umidade_relativa_pct
    + log(esforco_amostral), # log
  family = poisson(link = "log"),
  data = caf1,
  na.action = "na.fail"
)

# Teste de verossimilhanca entre os modelos nulo e completo COM esforco amostral passado como PREDITOR LIVRE
secao("Teste de verossimilhança entre os modelos nulo e completo COM medida de esforço amostral passada como PREDITOR LIVRE")
anova(mod02, mod04, test = "Chisq")

# Tabela de estimacao dos parametros do modelo completo COM esforco amostral passado como PREDITOR LIVRE
secao("Estimação dos parâmetros do modelo completo COM medida de esforço amostral passada como PREDITOR LIVRE")
summary(mod04)

# Parecer da comparacao do esforco amostral enquanto offset e enquanto preditor livre
cat("Ao observar o coeficiente estimado para o log do esforço amostral passado como preditor livre, observa-se que o coeficiente de é muito de 1. 
    Isso sugere ajustar um modelo com o coeficiente fixado em 1. Dessa forma, modelamos a taxa de captura de brocas mantendo, ao mesmo tempo, 
    a contagem como variável resposta para o modelo de Poisson. 
    Segundo Faraway (2006), esse tipo de abordagem é conhecido como modelo de taxa. Fixa-se o coeficiente em 1 utilizando um termo de offset. 
    Tal termo, presente no lado do preditor da equação do modelo, não possui parâmetro associado.")

# =======================================================================================
# MODELOS LINEARES GENERALIZADOS - POISSON
# Exaurindo as combinacoes entre as variaveis preditoras 
# =======================================================================================
secao("MODELOS LINEARES GENERALIZADOS - POISSON - Exaurindo as combinações entre as variáveis preditoras")

# Ajustar diferentes modelos de forma iterativa
tab01 <- dredge(mod03, extra = "R^2")

# Quantidade de modelos ajustados
dim(tab01)

# Identificar o melhor modelo com o criterio de informacao de akaike (AICc)
filter(tab01, AICc == min(AICc))

# Parecer
cat("\nConforme o AICc, o modelo que melhor se ajusta aos dados possui as seguintes covariáveis preditoras:
  - adubação, 
  - altitude,
  - declividade,
  - densidade de plantio,
  - idade da lavoura
  - matéria orgânica
  - ph solo,
  - precipitação, e
  - umidade."
)

# --------------------------------------
# mod05: modelo poisson com minimo AICc
# --------------------------------------
secao("Modelo poisson com minimo AICc")
mod05 <- glm(
  n_brocas_capturadas ~ adubacao_n_kg_ha 
    + altitude_m 
    + declividade_pct 
    + densidade_plantio 
    + idade_lavoura_anos 
    + materia_organica_pct 
    + ph_solo 
    + precipitacao_safra_mm 
    + umidade_relativa_pct, 
  family = poisson(link = "log"), 
  data = caf1, 
  na.action = "na.fail", 
  offset = log(esforco_amostral)
)

# OBS.: modelo poisson com minimo AICc (mod05) identifica-se com o modelo completo com medida de esforco amostral como offset (mod03)

# Teste de Verossimilhança entre os modelos nulo e poisson com minimo AICc
secao("Teste de verossimilhança entre os modelos poisson nulo e com minimo AICc")
anova(mod02, mod05, test = "Chisq")

# Reducao percentual da deviancia nula
residuos_mod05 <- ((deviance(mod02) - deviance(mod05)) / deviance(mod02)) * 100

# Parecer
cat("\nEm comparação ao modelo nulo, o modelo completo COM esforço amostral passado como OFFSET, com o AICc mínimo, reduziu 
    a deviância nula em ", residuos_mod05, "%, indicando que o modelo explicou esse percentual da deviância nula.")

# Analise da deviancia das covariaveis do modelo poisson com minimo AICc 
secao("Analise da deviância das covariáveis do modelo poisson com minimo AICc")
anova(mod05, test = "Chisq")

# Parecer 
cat("\nPara o modelo poisson com minimo AICc, em termos de redução de resíduos de deviância nula, as variáveis com redução significativa foram:
    - adubação, 
    - altitude,
    - densidade de plantio,
    - idade da lavoura,
    - matéria orgânica,
    - precipitação, e
    - umidade."
  )

# Faixa de coeficientes possiveis para os parametros do modelo poisson com minimo AICc
secao("Faixa de coeficientes possíveis para os parâmetros do modelo poisson com minimo AICc")
confint(mod05)

# Tabela de estimacao dos parametros
secao("Estimação dos parâmetros do modelo poisson com minimo AICc")
summary(mod05)

# Visualizar efeitos
plot(effects::allEffects(mod05))

# Parecer 
cat("\nPara o modelo poisson com minimo AICc, em termos de efeito, as variáveis com efeito significativo foram:
  - adubação, 
  - altitude, 
  - densidade de plantio, 
  - idade da lavoura, 
  - matéria orgânica, e 
  - umidade.
  Um modelo final será ajustado apenas com essas variáveis a fim de se observar se há um ajuste ainda melhor.")

# --------------------------------------
# mod06: modelo poisson final
# --------------------------------------
secao("Modelo poisson final")
mod06 <- glm(
  n_brocas_capturadas ~ adubacao_n_kg_ha 
    + altitude_m 
    + densidade_plantio 
    + idade_lavoura_anos 
    + materia_organica_pct 
    + umidade_relativa_pct, 
  family = poisson(link = "log"), 
  data = caf1, 
  na.action = "na.fail", 
  offset = log(esforco_amostral)
)

# Teste de Verossimilhança entre os modelos nulo e poisson final
secao("Teste de verossimilhança entre os modelos poisson nulo e final")
anova(mod02, mod06, test = "Chisq")

# Reducao percentual da deviancia nula
residuos_mod06 <- ((deviance(mod02) - deviance(mod06)) / deviance(mod02)) * 100

# Analise da deviancia das covariaveis do modelo poisson final 
secao("Analise da deviância das covariáveis do modelo poisson finalc")
anova(mod06, test = "Chisq")

# Faixa de coeficientes possiveis para os parametros do modelo poisson final
secao("Faixa de coeficientes possíveis para os parâmetros do modelo poisson final")
confint(mod06)

# Tabela de estimacao dos parametros
secao("Estimação dos parâmetros do modelo poisson final")
summary(mod06)

# Visualizar efeitos
plot(effects::allEffects(mod06))

# Parecer
cat("\nEm comparação ao modelo nulo, o modelo completo COM esforço amostral passado como OFFSET, com o AICc mínimo e com variaveis
    de efeitos significativo selecionadas, reduziu a deviância nula em ", residuos_mod06, "%, indicando que o modelo explicou esse 
    percentual da deviância nula. O modelo de AICc mínimo antes do teste de efeito significativo das variáveis reduziu a deviância em",
    residuos_mod05,". Apesar de uma redução maior, o modelo final com as variáveis sem efeito significativo removidas, resultou em um AICc 
    ligeiramente menor (2261.2) quando comparado aquele em que as mesmas foram mantidas (2263.1).")

# =======================================================================================
# DIAGNOSTICO E ANALISE DE RESIDUOS
# =======================================================================================
secao("DIAGNÓSTICO DOS RESÍDUOS")

<<<<<<< Updated upstream
plot(mod04) # (TO DO: verificar se podem diagnosticar glms?)
# Interpretacao:
# Residuals vs. Fitted: detecta se ha falta de ajuste e se a variancia e constante
# se os residuos mostram tendencia curvilinea, sinal de que pode haver relacoes nao-lineares

# Scale-Location: mostra se os residuos estao distribuidos igualmente para verificar
# se ha homocedasticidade (igual variancia) demonstrado por uma linha horizontal com pontos
# distribuidos de forma igual e aleatoria

# Q-Q plot
# Compara os residuos com observacoes ideais e mostra a distribuicao dos mesmos

# Residuals vs. Leverage
# Verifica se os valores extremos influenciam no modelo 

# Envelope simulado dos residuos do modelo poisson final
set.seed(20260922)   # envelope simulado: semente para o script ser reproduzivel
res01 <- hnp(mod04, plot.sim = FALSE)

# Transformar em dataframe
res01 <- data.frame(x = res01$x, median = res01$median,
                    lower = res01$lower, upper = res01$upper,
                    residuals = res01$residuals)

# Grafico do envelope simulado dos residuos
p00 <- ggplot(data = res01) +
    geom_ribbon(aes(x = x, ymin = lower, ymax = upper), alpha = 0.8) +
    geom_line(aes(x = x, y = median), colour = "white") +
    geom_point(aes(x = x, y = residuals), pch = 21, fill = "white",
               colour = "black", size = 5, alpha = 0.5) +
    labs(x = "Quantis teóricos", y = "Resíduos") +
    theme_gray(base_size = 18)
p00

# Exportar os residuos do modelo
res02 <- fortify(mod04)
res02$ID <- 1:nrow(res02)

# Grafico de dispersao dos residuos
p01 <- ggplot(data = res02, aes(x = ID, y = .stdresid)) +
    geom_point(pch = 21, fill = "white", colour = "black", size = 5,
               alpha = 0.8) +
    geom_hline(yintercept = 0, colour = "red") +
    labs(x = "Índice da amostra", y = "Resíduos padronizados") +
    theme_gray(base_size = 18)
p01

# Histograma dos residuos
p02 <- ggplot(data = res02, aes(x = .stdresid)) +
    geom_histogram(binwidth = 1, boundary = 1, closed = "right",
                   fill = "white", colour = "black") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 25)) +
    labs(x = "Resíduos padronizados", y = "Frequência") +
    theme_gray(base_size = 18)
p02  # TO DO Verificar 

# Visualizar os graficos lado a lado
(p01 | p02 | p00)
=======
par(mfrow = c(2, 2))
plot(mod06)
>>>>>>> Stashed changes

# =======================================================================================
# VERIFICACAO DA EQUIDISPERSAO DO MODELO POISSON
# =======================================================================================
secao("VERIFICAÇÃO DA EQUIDISPERSÃO")

#?dispersiontest

# Conforme a documentacao, poission assume que:
# a esperanca condicional E[y] = μ e a variancia VAR[y] = μ sao iguais

# --------------------------------------
# Calcular a esperanca e a variancia 
# e observar se sao iguais
# --------------------------------------

# Agrupar n_brocas_capturadas pela quantidade de vezes que foram capturadas
caf2 <- caf1 %>%
  group_by(n_brocas_capturadas) %>%
  summarize(count = n())

# Frequencia relativa
#####@> Frequencia relativa
caf2$n_brocas_capt_freq_relat <- caf2$count/sum(caf2$count)

# Media de brocas capturadas
n_brocas_capt_media <- sum(caf2$n_brocas_capturadas * caf2$n_brocas_capt_freq_relat)

# Variancia de brocas capturadas
n_brocas_capt_variancia <- sum((caf2$n_brocas_capturadas - n_brocas_capt_media)^2)/(sum(caf2$count)-1)

# Verificar igualdade entre media e variancia
n_brocas_capt_media == n_brocas_capt_variancia

# Divisao da variancia pela media
n_brocas_capt_variancia/n_brocas_capt_media

# TO DO Parecer
cat()

# --------------------------------------
# Verificar a frequencia esperada
# e observada com um modelo poisson
# --------------------------------------

# Frequencia esperada de acordo com um modelo poisson
n_brocas_capt_freq_pois <- dpois(x = caf2$n_brocas_capturadas, lambda = n_brocas_capt_media)

# Visualizar diferenca entre frequencia relativa esperada e observada no modelo poisson
barplot(t(cbind(caf2$n_brocas_capt_freq_relat, n_brocas_capt_freq_pois)),
        names.arg = caf2$n_brocas_capt_freq_relat,
        ylim = c(0, 1), xlab = "Número de brocas capturadas",
        ylab = "Frequência relativa", beside = TRUE,
        legend.text = c("Observado", "Esperado"),
        main = "Poisson")

# Grafico dos residuos no modelo poisson
plot(caf2$n_brocas_capturadas, n_brocas_capt_freq_pois - caf2$n_brocas_capt_freq_relat,
     xlim = c(0, 30), ylim = c(-0.6, 0.6),
     xlab = "Número de brocas capturadas", ylab = "Resíduos",
     axes = TRUE, main = "Resíduos Poisson", pch = 19)
abline(h = 0)

# Qui-quadrado para verificar se o modelo poisson se ajusta aos dados
qui_pois <- with(caf2,
                 sum((count - n_brocas_capt_freq_pois*sum(count))^2 /
                     (n_brocas_capt_freq_pois * sum(count))))
pchisq(q = qui_pois, df = nrow(caf2) - 1, lower.tail = FALSE)

# TO DO Parecer
cat("")

# --------------------------------------
# Teste formal de equidispersao
# --------------------------------------

# Razao do desvio residual e dos graus de liberdade dos residuos
deviance(mod06)
df.residual(mod06)
deviance(mod06)/df.residual(mod05)

# Testar sobre a hipotese alternativa da variancia ser uma funcao linear
dispersiontest(mod05, trafo = 1)

# Testar sobre a hipotese alternativa da variancia ser uma funcao quadratica
dispersiontest(mod05, trafo = 2)

# TO DO Parecer
cat("")

# --------------------------------------
# Verificar a frequencia esperada
# e observada com um modelo binomial
# negativo
# --------------------------------------

# Parametro de superdispersao k
# calculado pelo metodo dos momentos
# um parametro que entra na binomial negativa
k <- (n_brocas_capt_media ^ 2) / (n_brocas_capt_variancia - n_brocas_capt_media)

# Frequencia esperada de acordo com um modelo binomial negativa
n_brocas_capt_freq_binom <- dnbinom(x = caf2$n_brocas_capturadas, size = k, mu = n_brocas_capt_media)

# Visualizar diferenca entre frequencia relativa esperada e observada no modelo binomial negativa
barplot(t(cbind(caf2$n_brocas_capt_freq_relat, n_brocas_capt_freq_binom)),
        names.arg = caf2$n_brocas_capt_freq_relat,
        ylim = c(0, 1), xlab = "Número de brocas capturadas",
        ylab = "Frequência relativa", beside = TRUE,
        legend.text = c("Observado", "Esperado"),
        main = "Binomial Negativa")

# Qui-quadrado para verificar se o modelo binomial negativa se ajusta aos dados
plot(caf2$n_brocas_capturadas, n_brocas_capt_freq_binom - caf2$n_brocas_capt_freq_relat,
     xlim = c(0, 30), ylim = c(-0.6, 0.6),
     xlab = "Número de brocas capturadas", ylab = "Resíduos",
     axes = TRUE, main = "Resíduos Binomial Negativa", pch = 19)
abline(h = 0)

# Qui-quadrado para verificar se o modelo poisson se ajusta aos dados
qui_pois <- with(caf2,
                 sum((count - n_brocas_capt_freq_pois*sum(count))^2 /
                     (n_brocas_capt_freq_pois * sum(count))))
pchisq(q = qui_pois, df = nrow(caf2) - 1, lower.tail = FALSE)

# TO DO Parecer
cat("")

# =======================================================================================
# MODELOS LINEARES GENERALIZADOS - BINOMIAL NEGATIVA
# =======================================================================================
secao("MODELOS LINEARES GENERALIZADOS - BINOMIAL NEGATIVA")

# --------------------------------------
# Modelo nulo
# --------------------------------------
mod07 <- glm.nb(
  n_brocas_capturadas ~ 1,
  data = caf1,
  link = "log")

# --------------------------------------
# Modelo completo
# --------------------------------------
mod08 <- glm.nb(
  n_brocas_capturadas ~ adubacao_n_kg_ha
    + altitude_m
    + declividade_pct
    + densidade_plantio
    + idade_lavoura_anos
    + materia_organica_pct
    + ph_solo
    + precipitacao_safra_mm
    + umidade_relativa_pct
    + offset(log(esforco_amostral)), 
  data = caf1,
  na.action = "na.fail"
)

# Teste de verossimilhanca entre os modelos nulo e completo COM esforco amostral passado como OFFSET
secao("Teste de verossimilhança entre os modelos nulo e completo COM medida de esforço amostral passada como OFFSET")
anova(mod07, mod08, test = "Chisq")
residuos_mod08 <- ((deviance(mod07) - deviance(mod08)) / deviance(mod07)) * 100

# Parecer
cat("\nEm comparação ao modelo nulo, o modelo completo COM esforço amostral passado como OFFSET reduziu 
    os resíduos em ", residuos_mod08, "%.")

# Analise da deviancia das covariaveis do modelo completo COM esforco amostral passado como OFFSET
secao("Analise da deviância das covariáveis do modelo completo COM medida de esforço amostral")
anova(mod08, test = "Chisq")

# Parecer 
cat("\nPara o modelo completo COM esforço amostral passado como OFFSET, em termos de redução de resíduos de deviância nula, 
    as variáveis com redução significativa foram todas, exceto declividade, idade da lavoura e ph solo.")

# Tabela de estimacao dos parametros do modelo completo COM esforco amostral passado como OFFSET
secao("Estimação dos parâmetros do modelo completo COM medida de esforço amostral passada como OFFSET")
summary(mod08)

# Parecer 
cat("\nPara o modelo completo COM esforço amostral passado como OFFSET, em termos de efeito, as variáveis com 
    efeito significativo foram todas, exceto declividade, ph solo e precipitação.")

# =======================================================================================
# MODELOS LINEARES GENERALIZADOS - BINOMIAL NEGATIVA
# Exaurindo as combinacoes entre as variaveis preditoras 
# =======================================================================================
secao("MODELOS LINEARES GENERALIZADOS - BINOMIAL NEGATIVA - Exaurindo as combinações entre as variáveis preditoras")

# Ajustar diferentes modelos de forma iterativa
tab02 <- dredge(mod08, extra = "R^2")

# Quantidade de modelos ajustados
dim(tab02)

# Identificar o melhor modelo com o criterio de informacao de akaike (AICc)
filter(tab02, AICc == min(AICc))

# Parecer
cat("\nConforme o AICc, o modelo que melhor se ajusta aos dados possui as seguintes covariáveis preditoras:
  - adubação, 
  - altitude,
  - declividade,
  - densidade de plantio,
  - idade da lavoura
  - matéria orgânica
  - ph solo,
  - precipitação, e
  - umidade."
)

# --------------------------------------
# Modelo binomial negativa 
# selecionado
# --------------------------------------
secao("Modelo binomial negativa final")
mod09 <- glm.nb(
  n_brocas_capturadas ~ adubacao_n_kg_ha 
    + altitude_m 
    + declividade_pct 
    + densidade_plantio 
    + idade_lavoura_anos 
    + materia_organica_pct 
    + ph_solo 
    + precipitacao_safra_mm 
    + umidade_relativa_pct 
    + offset(log(esforco_amostral)), 
  data = caf1, 
  na.action = "na.fail",
  link = log,
  init.theta = 4.377708848
) 

# Teste de Verossimilhança entre os modelos binomial negativa nulo e modelo final
secao("Teste de verossimilhança entre os modelos binomial negativa nulo e final")
anova(mod07, mod09, test = "Chisq")
residuos_mod09 <- ((deviance(mod07) - deviance(mod09)) / deviance(mod07)) * 100

# Parecer
cat("\nEm comparação ao modelo nulo, o modelo binomial negativa reduziu os resíduos em ", residuos_mod09, "%.")

# Analise da deviancia das covariaveis do modelo binomial negativa final 
secao("Analise da deviância das covariáveis do modelo binomial negativa final")
anova(mod09, test = "Chisq")

# Parecer 
cat("\nPara o modelo binomial negativa final, em termos de redução de resíduos de deviância nula, as variáveis com redução significativa foram:
    - adubação, 
    - altitude,
    - densidade de plantio,
    - matéria orgânica,
    - precipitação, e
    - umidade."
  )

# Faixa de coeficientes possiveis para os parametros do modelo binomial negativa final
secao("Faixa de coeficientes possíveis para os parâmetros do modelo binomial negativa final")
confint(mod09)

# Tabela de estimacao dos parametros
secao("Estimação dos parâmetros do modelo binomial negativa final")
summary(mod09)

# Visualizar efeitos
plot(effects::allEffects(mod09))

# Parecer 
cat("\nPara o modelo binomial negativa final, em termos de efeito, as variáveis com efeito significativo foram:
  - adubação, 
  - altitude, 
  - densidade de plantio, 
  - idade da lavoura, 
  - matéria orgânica, e 
  - umidade.")

# --------------------------------------
# Modelo binomial negativa final 
# --------------------------------------

secao("Modelo binomial negativa final")
mod10 <- glm.nb(
  n_brocas_capturadas ~ adubacao_n_kg_ha 
    + altitude_m 
    + densidade_plantio 
    + idade_lavoura_anos 
    + materia_organica_pct 
    + umidade_relativa_pct 
    + offset(log(esforco_amostral)), 
  data = caf1, 
  na.action = "na.fail",
  link = log,
  init.theta = 4.377708848
) 

# Teste de Verossimilhança entre os modelos binomial negativa nulo e modelo final
secao("Teste de verossimilhança entre os modelos binomial negativa nulo e final")
anova(mod07, mod10, test = "Chisq")
residuos_mod10 <- ((deviance(mod07) - deviance(mod10)) / deviance(mod07)) * 100

# Parecer
cat("\nEm comparação ao modelo nulo, o modelo binomial negativa reduziu os resíduos em ", residuos_mod10, "%.")

# Analise da deviancia das covariaveis do modelo binomial negativa final 
secao("Analise da deviância das covariáveis do modelo binomial negativa final")
anova(mod10, test = "Chisq")

# Parecer 
cat("\nPara o modelo binomial negativa final, em termos de redução de resíduos de deviância nula, as variáveis com redução significativa foram:
    - adubação, 
    - altitude,
    - densidade de plantio,
    - idade da lavoura
    - matéria orgânica, e
    - umidade."
  )

# Faixa de coeficientes possiveis para os parametros do modelo binomial negativa final
secao("Faixa de coeficientes possíveis para os parâmetros do modelo binomial negativa final")
confint(mod10)

# Tabela de estimacao dos parametros
secao("Estimação dos parâmetros do modelo binomial negativa final")
summary(mod10)

# Visualizar efeitos
plot(effects::allEffects(mod10))

# Multicolinearidade entre as covariaveis (requisito comum, secao 4b do enunciado)
secao("Multicolinearidade do modelo binomial negativa final (VIF)")
vif_mod08 <- car::vif(mod08)
print(round(vif_mod08, 3))
cat("\nMaior VIF:", round(max(vif_mod08), 3),
    "- limiares de referencia: 5 (atencao) e 10 (grave).\n")
cat(if (max(vif_mod08) < 5)
    "Abaixo de 5: nao ha inflacao de variancia que comprometa a estimacao\ndos coeficientes nem a selecao de variaveis.\n"
    else "ATENCAO: VIF acima do limiar - revisar a selecao de variaveis.\n")

# =======================================================================================
# DIAGNOSTICO E ANALISE DE RESIDUOS
# =======================================================================================
secao("DIAGNÓSTICO DOS RESÍDUOS")

plot(mod10) 

<<<<<<< Updated upstream
# Envelope simulado dos residuos do modelo poisson final
set.seed(20260922)   # envelope simulado: semente para o script ser reproduzivel
res05 <- hnp(mod08, plot.sim = FALSE)
=======
# Envelope simulado dos residuos do modelo binomial negativa
res05 <- hnp(mod10, plot.sim = FALSE)
>>>>>>> Stashed changes

# Transformar em dataframe
res05 <- data.frame(x = res05$x, median = res05$median,
                    lower = res05$lower, upper = res05$upper,
                    residuals = res05$residuals)

# Grafico do envelope simulado dos residuos
p06 <- ggplot(data = res05) +
    geom_ribbon(aes(x = x, ymin = lower, ymax = upper), alpha = 0.8) +
    geom_line(aes(x = x, y = median), colour = "white") +
    geom_point(aes(x = x, y = residuals), pch = 21, fill = "white",
               colour = "black", size = 5, alpha = 0.5) +
    labs(x = "Quantis teóricos", y = "Resíduos") +
    theme_gray(base_size = 18)
p06

# Exportar os residuos do modelo
res06 <- fortify(mod10)
res06$ID <- 1:nrow(res06)

# Grafico de dispersao dos residuos
p07 <- ggplot(data = res06, aes(x = ID, y = .stdresid)) +
    geom_point(pch = 21, fill = "white", colour = "black", size = 5,
               alpha = 0.8) +
    geom_hline(yintercept = 0, colour = "red") +
    labs(x = "Índice da amostra", y = "Resíduos padronizados") +
    theme_gray(base_size = 18)
p07

# Histograma dos residuos
p08 <- ggplot(data = res06, aes(x = .stdresid)) +
    geom_histogram(binwidth = 1, boundary = 1, closed = "right",
                   fill = "white", colour = "black") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 25)) +
    labs(x = "Resíduos padronizados", y = "Frequência") +
    theme_gray(base_size = 18)
p08 # TO DO Verificar 

# Visualizar os graficos lado a lado
(p08 | p07 | p06)

# ============================================================
# PARECERES
# ============================================================
secao("PARECERES")

# Pergunta do corpo tecnico:
# Quais condicoes de talhao favorecem a infestacao pela broca, 
# e quanto se ganha em pressao de praga ao alterar cada uma delas?

# Como a binomial negativa tem o log como funcao de ligacao, para se saber quanto se ganha em pressao
# de praga para cada condicao, deve-se fazer o inverso do log que e o exponencial e aplica-lo ao 
# coeficiente de cada variavel que entrou no modelo final

# --------------------------------------
# Parecer automatizado 
# --------------------------------------
secao("Parecer automatizado")

parecer_modelo <- function(modelo) {

  # Coeficientes e intervalos de confianca
  coeficientes <- coef(modelo)
  ic <- suppressMessages(confint(modelo))

  # Resumo do modelo
  resumo <- summary(modelo)$coefficients

  # Remover intercepto e offset
  vars <- names(coeficientes)
  vars <- vars[vars != "(Intercept)" &
               !grepl("offset", vars, fixed = TRUE)]

  # Legendas das variaveis
  nomes <- c(
    adubacao_n_kg_ha = "Adubação nitrogenada (kg/ha)",
    altitude_m = "Altitude (m)",
    declividade_pct = "Declividade",
    densidade_plantio = "Densidade de plantio (plantas/ha)",
    idade_lavoura_anos = "Idade da lavoura (anos)",
    materia_organica_pct = "Matéria orgânica do solo (%)",
    ph_solo = "pH do solo",
    precipitacao_safra_mm = "Precipitação da safra (mm)",
    umidade_relativa_pct = "Umidade relativa (%)"
  )

  # Calcular efeitos
  resultados <- data.frame(
    variavel = vars,
    nome = nomes[vars],
    coeficiente = coeficientes[vars],
    razao_taxa = exp(coeficientes[vars]),
    percentual = (exp(coeficientes[vars]) - 1) * 100,
    IC95_inferior = exp(ic[vars, 1]),
    IC95_superior = exp(ic[vars, 2]),
    p_valor = resumo[vars, 4],
    significativo = resumo[vars, 4] < 0.05
  )

  # Imprimir valores numericos
  print(resultados, row.names = FALSE)

  cat("\n\nPARECER\n\n")

  for (i in seq_len(nrow(resultados))) {

    r <- resultados[i, ]

    direcao <- ifelse(r$percentual > 0, "aumenta", "reduz")
    percentual <- abs(r$percentual)

    cat(
      sprintf(
        "- %s: cada aumento de uma unidade %s a pressão de praga em aproximadamente %.2f%% (razão de taxa = %.4f; IC95%%: %.4f–%.4f; p = %.4f).\n",
        r$nome,
        direcao,
        percentual,
        r$razao_taxa,
        r$IC95_inferior,
        r$IC95_superior,
        r$p_valor
      )
    )
  }

  cat("\n")

  # Variaveis nao significativas
  nao_significativos <- resultados[
    !resultados$significativo, ]

  if (nrow(nao_significativos) > 0) {

    cat(
      "Variáveis sem efeito estatisticamente significativo ",
      "sobre a taxa de captura (p >= 0,05):\n",
      sep = ""
    )

    cat(
      paste(nao_significativos$nome, collapse = ", "),
      ".\n\n",
      sep = ""
    )
  }

  # Significant variables
  significativos <- resultados[
    resultados$significativo, ]

  if (nrow(significativos) > 0) {

    positivos <- significativos[
      significativos$percentual > 0, ]

    negativos <- significativos[
      significativos$percentual < 0, ]

    cat("Síntese dos efeitos estatisticamente significativos:\n")

    if (nrow(positivos) > 0) {
      cat(
        "- Efeitos positivos: ",
        paste(positivos$nome, collapse = ", "),
        ".\n",
        sep = ""
      )
    }

    if (nrow(negativos) > 0) {
      cat(
        "- Efeitos negativos: ",
        paste(negativos$nome, collapse = ", "),
        ".\n",
        sep = ""
      )
    }
  }

  # Retornar a tabela numerica
  invisible(resultados)
}

resultados_parecer <- parecer_modelo(mod10)

# =======================================================================================
# PARECER FINAL
# =======================================================================================
secao("PARECER FINAL")

cat("Conforme o modelo binomial negativa final (mod10), as condições de talhão que favorecem 
a infestação pela broca e o quanto se ganha ou se perde em pressão de praga ao alterar cada uma delas, 
mantendo as demais variáveis e o esforço amostral constantes, são:

  - Adubação nitrogenada: cada aumento de 1 kg/ha reduz a pressão de praga em 
    aproximadamente 0,08% (razão de taxa = 0.9992). Ou seja, a adubação nitrogenada 
    atua como fator de proteção, e não de risco.

  - Altitude: cada aumento de 1 metro reduz a taxa de captura em aproximadamente 0.11% 
    (razão de taxa = 0.9989). Talhões em altitudes maiores tendem a apresentar menor 
    pressão de praga.

  - Densidade de plantio: cada planta adicional por hectare aumenta a pressão de praga 
    em aproximadamente 0.02% (razão de taxa = 1.0002). O efeito é estatisticamente 
    significativo, porém de magnitude muito pequena por unidade.

  - Idade da lavoura: cada ano adicional de idade aumenta a pressão de praga em 
    aproximadamente 1.73% (razão de taxa = 1.0173).

  - Matéria orgânica do solo: cada aumento de 1 ponto percentual de matéria orgânica 
    aumenta a pressão de praga em aproximadamente 14.04% (razão de taxa = 1.1403). Este 
    é um dos efeitos de maior magnitude entre os preditores significativos.

  - Umidade relativa: cada aumento de 1 ponto percentual de umidade aumenta a pressão de 
    praga em aproximadamente 5.00% (razão de taxa = 1.0497), consistente com a 
    associação moderada e positiva identificada na análise exploratória.

  - Declividade, pH do solo e precipitação da safra não apresentaram efeito 
    estatisticamente significativo sobre a taxa de captura (intervalos de confiança 
    de 95% contendo o valor 1).

Em síntese, matéria orgânica e umidade relativa são as condições de talhão que mais 
favorecem a infestação pela broca, enquanto adubação nitrogenada e altitude atuam 
como fatores de proteção. Densidade de plantio e idade da lavoura apresentam efeitos 
positivos, porém de magnitude reduzida.")

# =======================================================================================
#                  Creative Commons License 4.0
#                       (CC BY-NC-SA 4.0)
#
#  This is a humam-readable summary of (and not a substitute for) the
#  license (https://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
#
#  You are free to:
#
#  Share - copy and redistribute the material in any medium or format.
#
#  The licensor cannot revoke these freedoms as long as you follow the
#  license terms.
#
#  Under the following terms:
#
#  Attribution - You must give appropriate credit, provide a link to
#  license, and indicate if changes were made. You may do so in any
#  reasonable manner, but not in any way that suggests the licensor
#  endorses you or your use.
#
#  NonCommercial - You may not use the material for commercial
#  purposes.
#
#  ShareAlike - If you remix, transform, or build upon the material,
#  you must distributive your contributions under the same license
#  as the  original.
#
#  No additional restrictions — You may not apply legal terms or
#  technological measures that legally restrict others from doing
#  anything the license permits.
# =======================================================================================