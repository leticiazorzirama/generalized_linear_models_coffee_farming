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
       datasets, DCluster, stargazer, olsrr, performance, report, faraway)

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
cat("\nEm comparação ao modelo nulo (mod00), o modelo completo SEM esforço amostral (mod01) reduziu a deviância nula em ",
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
cat("\nEm comparação ao modelo nulo (mod02), o modelo completo COM esforço amostral passado como OFFSET (mod03) reduziu 
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
    Isso sugere ajustar um modelo com o coeficiente fixado em 1. Dessa forma, modelou-se a taxa de captura de brocas mantendo, ao mesmo tempo, 
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

# Teste de Verossimilhanca entre os modelos nulo e poisson com minimo AICc
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

# Teste de Verossimilhanca entre os modelos nulo e poisson final
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
cat("\nEm comparação ao modelo nulo (mod02), o modelo completo COM esforço amostral passado como OFFSET, com o AICc mínimo e com variáveis
    de efeitos significativo selecionadas (mod06), reduziu a deviância nula em ", residuos_mod06, "%, indicando que o modelo explicou esse 
    percentual da deviância nula. O modelo de AICc mínimo antes do teste de efeito significativo das variáveis (mod05) reduziu a deviância em",
    residuos_mod05,". Apesar de uma redução maior, o modelo final com as variáveis sem efeito significativo removidas, resultou em um AICc 
    ligeiramente menor (2261.2) quando comparado aquele em que as mesmas foram mantidas (2263.1).")

# =======================================================================================
# DIAGNOSTICO E ANALISE DE RESIDUOS
# =======================================================================================
secao("DIAGNÓSTICO DOS RESÍDUOS")

# Detectar valores atipicos
# Faraway (2006) recomenda detectar valores atipicos com distribuicao meia-normal
# porque a resolucao do plot e dobrada ao ter todos os pontos em uma cauda

# Checar valores atipicos na variavel resposta
# Rstudent
halfnorm(rstudent(mod06))
# Como interpretar
# Os pontos devem se alinhar de forma proxima a uma trajetoria retilinea
# Por se tratar de um modelo de Poisson, pequenas variacoes sao comuns, mas a 
# estrutura geral deve ser uniforme
# Quaisquer pontos isolados no canto superior direito que se afastem da tendencia
# linear e apresentem um salto acentuado para cima sao valores atipicos
# Geralmente, residuos estudentizados com valores absolutos superiores a 2 ou 3 
# justificam uma inspecao cuidadosa, pois representam contagens que o modelo de Poisson 
# nao conseguiu prever com precisao.

# Em outras palavras
# O modelo tenta prever quantas brocas vão ser capturadas em cada talhão, com base nas condições dele 
# (umidade, altitude, etc.). Para cada talhão, compara-se o valor real capturado com o valor que o modelo previu.
# Se a diferença for pequena, tudo certo — o modelo acertou bem.
# Se a diferença for gigante (o modelo previu 5 brocas e apareceram 40), esse talhão é um "ponto estranho" nesse gráfico. 
# Ele está indicando: "a minha resposta (número de brocas) não combinou com o que era esperado."

# Parecer
cat("Presença de valores atípicos preditos, principalmente as predições do número de brocas para os talhões
    321 e 287.")

# Checar valores atipicos nos coeficientes das variaveis preditoras
# Leverage
halfnorm(influence(mod06)$hat) 
# Como interpretar
# A maioria dos pontos deve estar agrupada na parte inferior esquerda, 
# representando configuracoes de dados padrao 
# Pontos na parte superior direita que estejam separados do restante do conjunto 
# por um intervalo visual perceptivel são pontos de alta alavancagem 
# Uma regra pratica geral para o limite de alta alavancagem e 2p/n 
# p = numero de parametros estimados 
# n = tamanho da amostra
# Se um ponto se desviar visualmente da linha de meia-normal e exceder esse valor, 
# seus valores sao extremos
p <- 6
n <- nrow(caf)
alavanca <- 2 * p / n

# Em outras palavras
# "Esse talhão é diferente dos outros"
# Não se está olhando pro resultado (quantas brocas), e sim pras características do talhão: a altitude dele, a umidade, o pH, etc.
# Se um talhão tem uma combinação de características bem rara — tipo, é o único com altitude super alta E super seco ao mesmo tempo, 
# enquanto todos os outros são parecidos entre si — esse talhão tem alta alavancagem. 
# Ele é "diferente" porque ele mesmo, como talhão, já é fora do padrão.

cat("Para um ponto de alavancagem de", alavanca, "contata-se talhões atípicos nos, sobretudo os talhões 214 e 148.")

# Checar o quanto os valores atípicos influenciam no modelo 
# A distancia de Cook combina magnitude residual e alavancagem para medir a influencia global
# Ela quantifica o quanto as estimativas dos coeficientes do modelo mudariam se aquela observacao
# especifica fosse completamente alterada
# Distancia de cook 
halfnorm(cooks.distance(mod06)) 
# Como analisar
# Pontos individuais no canto superior direito que se destaquem muito acima de todas as outras observacoes
# Se um ponto apresentar valores elevados tanto de alavancagem quanto de residuo, sua distancia de Cook ira disparar
# Em um grafico de meia-normal, observar a distancia relativa em vez de limites absolutos fixos (como 1,0)
# Se um ou dois pontos estiverem situados muito acima do restante da curva, eles estarao exercendo uma influencia significativa 
# sobre os parametros do modelo poisson.

# Em outras palavras
# "Esse ponto realmente bagunçou o resultado final"
# Essa é a junção das duas ideias anteriores. 
# A distância de Cook pergunta: "e se eu tirar esse talhão da análise, o modelo muda muito ou fica quase igual?"
# Ela combina:
# o quão estranho foi o resultado dele (resíduo), com
# o quão diferente ele é dos outros (alavancagem).
# Se um ponto tem as duas coisas ao mesmo tempo — um resultado estranho e características raras — ele pode estar puxando o 
# modelo inteiro para um lado, tipo um aluno que, se ele saísse da turma, a média da sala mudaria bastante.
# É o ponto que realmente tem poder de influenciar a conclusão final, não só de "chamar atenção".

# Parecer
cat("O talhão 321 apresentou uma distância de Cook acima do restante da curva, 
    indicando uma influência global na capacidade preditiva do modelo.")

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
caf2$n_brocas_capt_freq_relat <- caf2$count/sum(caf2$count)

# Media de brocas capturadas
n_brocas_capt_media <- sum(caf2$n_brocas_capturadas * caf2$n_brocas_capt_freq_relat)

# Variancia de brocas capturadas
n_brocas_capt_variancia <- sum((caf2$n_brocas_capturadas - n_brocas_capt_media)^2)/(sum(caf2$count)-1)

# Verificar igualdade entre media e variancia
n_brocas_capt_media == n_brocas_capt_variancia

# Divisao da variancia pela media
brocas_var_media <- n_brocas_capt_variancia/n_brocas_capt_media

# Parecer
cat("Diante da proporção de", brocas_var_media, "entre a variância e a média do número de brocas no nível 
    marginal dos dados brutos, as contagens apresentam mais variabilidade do que uma distribuição de Poisson 
    permitiria.")

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

# Parecer
cat("As inspeções visuais também mostram uma desencaixe entre os valores observados e esperados para um 
    modelo poisson.")

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

# Testar a dispersao com pearson
mod06_dis_pear <- sum(residuals(mod06, type = "pearson")^2) / df.residual(mod06)
summary(mod06, dispersion = mod06_dis_pear)
# Como interpretar
# dp ≈ 1 = dispersao consistente com poisson
# dp > 1 = superdispersao
# dp < 1 = subdispersao
?dispersiontest

# Parecer
cat("Os testes formais de equidispersão evidenciam a superdispersão nos dados.
    O teste formal de dispersão indicou um nível de significância superior a zero
    e a variância sendo uma função quadrática. Neste caso, recomenda-se a binomial negativa.")

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

# Parecer
cat("As inspeções visuais mostram um melhor encaixe entre os valores observados e esperados para um 
    modelo binomial negativa, bem como uma diminuição dos resíduos.")

# =======================================================================================
# MODELOS LINEARES GENERALIZADOS - BINOMIAL NEGATIVA
# mod07: modelo nulo
# mod08: modelo completo
# mod09: modelo com minimo AICC
# mod10: modelo final
# =======================================================================================
secao("MODELOS LINEARES GENERALIZADOS - BINOMIAL NEGATIVA")

# --------------------------------------
# mod07: modelo nulo
# --------------------------------------
mod07 <- glm.nb(
  n_brocas_capturadas ~ 1 + offset(log(esforco_amostral)),
  data = caf1,
  na.action = "na.fail"
)

# --------------------------------------
# mod08: modelo completo
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

# Calculo manual das deviances (2 x log-lik) e do LR stat, para conferencia do teste acima
logLik_mod07 <- as.numeric(logLik(mod07))
logLik_mod08 <- as.numeric(logLik(mod08))

dev_mod07 <- 2 * logLik_mod07   # equivalente a "2 x log-lik." do mod07 no anova()
dev_mod08 <- 2 * logLik_mod08   # equivalente a "2 x log-lik." do mod08 no anova()

LR_stat_mod07_mod08 <- dev_mod08 - dev_mod07
df_LR_mod07_mod08 <- length(coef(mod08)) - length(coef(mod07))
p_valor_mod07_mod08 <- pchisq(LR_stat_mod07_mod08, df = df_LR_mod07_mod08, lower.tail = FALSE)

cat("\n2 x log-lik. mod07 (nulo):   ", round(dev_mod07, 3),
    "\n2 x log-lik. mod08 (completo):", round(dev_mod08, 3),
    "\nLR stat. (diferenca):         ", round(LR_stat_mod07_mod08, 4),
    "\ngraus de liberdade:           ", df_LR_mod07_mod08,
    "\np-valor:                      ", format.pval(p_valor_mod07_mod08, digits = 4, eps = .Machine$double.eps))

# Parecer
cat("\nO teste de razão de verossimilhanças entre o modelo nulo (mod07) e o modelo completo (mod08),
    ambos COM esforço amostral passado como OFFSET, resultou em 2 x log-lik. =", dev_mod07,
    "para o modelo nulo e", dev_mod08, "para o modelo completo, cuja diferença fornece LR stat. =", 
    LR_stat_mod07_mod08 ,"(9 g.l.), com p-valor <", p_valor_mod07_mod08, "(essencialmente 0, confirmado pelo cálculo manual acima). 
    Portanto, rejeita-se a hipótese nula de que os coeficientes das covariáveis agronômicas e ambientais são simultaneamente
    iguais a zero, ou seja, o modelo completo apresenta redução altamente significativa da deviância em
    relação ao modelo nulo. O conjunto de covariáveis incluídas contribui de forma significativa para
    explicar a variação no número de brocas capturadas, e o modelo completo (mod08) deve ser preferido
    em relação ao nulo (mod07).")

# Analise da deviancia das covariaveis do modelo completo COM esforco amostral passado como OFFSET
secao("Analise da deviância das covariáveis do modelo completo COM medida de esforço amostral")
anova(mod08, test = "Chisq")

# Parecer
cat("\nPara o modelo completo COM esforço amostral passado como OFFSET, em termos de redução de resíduos, 
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
# mod09: modelo binomial negativa 
# com minimo AICc
# --------------------------------------
secao("Modelo binomial negativa com minimo AICc")
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

# Teste de Verossimilhanca entre os modelos binomial negativa nulo e modelo com minimo AICc
secao("Teste de verossimilhança entre os modelos binomial negativa nulo e com minimo AICc")
anova(mod07, mod09, test = "Chisq")

# Calculo manual das deviances (2 x log-lik) e do LR stat, para conferencia do teste acima
logLik_mod09 <- as.numeric(logLik(mod09))

dev_mod09 <- 2 * logLik_mod09   # equivalente a "2 x log-lik." do mod09 no anova()

LR_stat_mod07_mod09 <- dev_mod09 - dev_mod07
df_LR_mod07_mod09 <- length(coef(mod09)) - length(coef(mod07))
p_valor_mod07_mod09 <- pchisq(LR_stat_mod07_mod09, df = df_LR_mod07_mod09, lower.tail = FALSE)

cat("\n2 x log-lik. mod07 (nulo):        ", round(dev_mod07, 3),
    "\n2 x log-lik. mod09 (minimo AICc):  ", round(dev_mod09, 3),
    "\nLR stat. (diferenca):              ", round(LR_stat_mod07_mod09, 4),
    "\ngraus de liberdade:                ", df_LR_mod07_mod09,
    "\np-valor:                           ", format.pval(p_valor_mod07_mod09, digits = 4, eps = .Machine$double.eps))

# Parecer
cat("\nO teste de razão de verossimilhanças entre o modelo nulo (mod07) e o modelo com mínimo AICc (mod09),
    ambos COM esforço amostral passado como OFFSET, resultou em 2 x log-lik. =", dev_mod07,
    "para o modelo nulo e", dev_mod09, "para o modelo mod09, cuja diferença fornece LR stat. =", 
    LR_stat_mod07_mod09 ,"(9 g.l.), com p-valor <", p_valor_mod07_mod09, "(essencialmente 0, confirmado pelo cálculo manual acima). 
    Portanto, rejeita-se a hipótese nula de que os coeficientes das covariáveis selecionadas por mínimo
    AICc são simultaneamente iguais a zero, ou seja, o modelo mod09 apresenta redução altamente
    significativa da deviância em relação ao modelo nulo, devendo ser preferido a este. Note-se que os
    valores de deviância e de LR stat. de mod09 coincidem com os de mod08 (mesmo conjunto de 9
    covariáveis), diferindo apenas na estimativa inicial de theta fornecida.")

# Analise da deviancia das covariaveis do modelo binomial negativa com minimo AICc 
secao("Analise da deviância das covariáveis do modelo binomial negativa com minimo AICc")
anova(mod09, test = "Chisq")

# Parecer 
cat("\nPara o modelo binomial negativa com minimo AICc, em termos de redução de resíduos, as variáveis com redução significativa foram:
    - adubação, 
    - altitude,
    - densidade de plantio,
    - matéria orgânica,
    - precipitação, e
    - umidade."
  )

# Faixa de coeficientes possiveis para os parametros do modelo binomial negativa com minimo AICc
secao("Faixa de coeficientes possíveis para os parâmetros do modelo binomial negativa com minimo AICc")
confint(mod09)

# Tabela de estimacao dos parametros
secao("Estimação dos parâmetros do modelo binomial negativa com minimo AICc")
summary(mod09)

# Visualizar efeitos
plot(effects::allEffects(mod09))

# Parecer 
cat("\nPara o modelo binomial negativa com minimo AICc, em termos de efeito, as variáveis com efeito significativo foram:
  - adubação, 
  - altitude, 
  - densidade de plantio, 
  - idade da lavoura, 
  - matéria orgânica, e 
  - umidade.
  Um modelo final será ajustado apenas com essas variáveis a fim de se observar se há um ajuste ainda melhor.")

# --------------------------------------
# mod10: modelo binomial negativa final 
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

# Teste de Verossimilhanca entre os modelos binomial negativa nulo e modelo final
secao("Teste de verossimilhança entre os modelos binomial negativa nulo e final")
anova(mod07, mod10, test = "Chisq")

# Calculo manual das deviances (2 x log-lik) e do LR stat, para conferencia do teste acima
logLik_mod10 <- as.numeric(logLik(mod10))
dev_mod10 <- 2 * logLik_mod10   # equivalente a "2 x log-lik." do mod10 no anova()

LR_stat_mod07_mod10 <- dev_mod10 - dev_mod07
df_LR_mod07_mod10 <- length(coef(mod10)) - length(coef(mod07))
p_valor_mod07_mod10 <- pchisq(LR_stat_mod07_mod10, df = df_LR_mod07_mod10, lower.tail = FALSE)

cat("\n2 x log-lik. mod07 (nulo): ", round(dev_mod07, 3),
    "\n2 x log-lik. mod10 (final): ", round(dev_mod10, 3),
    "\nLR stat. (diferenca):       ", round(LR_stat_mod07_mod10, 4),
    "\ngraus de liberdade:         ", df_LR_mod07_mod10,
    "\np-valor:                    ", format.pval(p_valor_mod07_mod10, digits = 4, eps = .Machine$double.eps))

# Parecer
cat("\nO teste de razão de verossimilhanças entre o modelo nulo (mod07) e o modelo final (mod10),
    ambos COM esforço amostral passado como OFFSET, resultou em 2 x log-lik. =", dev_mod07,
    "para o modelo nulo e", dev_mod10, "para o modelo mod10, cuja diferença fornece LR stat. =", 
    LR_stat_mod07_mod10 ,"(6 g.l.), com p-valor <", p_valor_mod07_mod10, "(essencialmente 0, confirmado pelo cálculo manual acima). 
    Portanto, rejeita-se a hipótese nula de que os coeficientes das seis covariáveis do modelo final são
    simultaneamente iguais a zero, ou seja, o modelo mod10 apresenta redução altamente significativa da
    deviância em relação ao modelo nulo, devendo ser preferido a este. Vale notar que, com apenas 6
    covariáveis (3 a menos que mod09), mod10 obteve um LR stat. praticamente idêntico ao de mod09,
    indicando que a remoção de declividade, ph do solo e precipitação não comprometeu de forma relevante
    a capacidade explicativa do modelo, o que é consistente com a escolha de mod10 como modelo mais
    parcimonioso.")

# Analise da deviancia das covariaveis do modelo binomial negativa final 
secao("Analise da deviância das covariáveis do modelo binomial negativa final")
anova(mod10, test = "Chisq")

# Parecer 
cat("\nPara o modelo binomial negativa final, com o AICc mínimo e com variaveis
    de efeitos significativo selecionadas, as variáveis com redução significativa foram:
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

# Parecer
cat("\nEm comparação ao modelo completo (mod09), o modelo com o AICc mínimo e com variáveis de efeitos significativo selecionadas (mod10), 
    apresentou uma ligeira diferença na verossimilhança. O mod10 resultou em um AICc ligeiramente menor (1950) quando comparado com o 
    modelo 09 cujas variáveis sem efeito significativo foram mantidas (1956).")

# =======================================================================================
# DIAGNOSTICO E ANALISE DE RESIDUOS
# =======================================================================================
secao("DIAGNÓSTICO DOS RESÍDUOS")

# Checar valores atipicos na variavel resposta
# Rstudent
halfnorm(rstudent(mod10))

# Parecer
cat("Presença de valores atípicos preditos, principalmente as predições do número de brocas para os talhões
    213 e 287.")

# Checar valores atipicos nos coeficientes das variaveis preditoras
# Leverage
halfnorm(influence(mod10)$hat) 
p <- 6
n <- nrow(caf)
alavanca <- 2 * p / n

cat("Para um ponto de alavancagem de", alavanca, "contata-se talhões atípicos nos, sobretudo os talhões 85 e 148.")

# Distancia de cook 
halfnorm(cooks.distance(mod10)) 

# Parecer
cat("Os talhões 84 e 272 apresentaram uma distância de Cook acima do restante da curva, 
    indicando uma influência global na capacidade preditiva do modelo.")

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