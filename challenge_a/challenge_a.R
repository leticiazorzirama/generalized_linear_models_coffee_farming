# ==========================================================================
# UNIVERSIDADE DO VALE DO ITAJAI - UNIVALI
# ESCOLA POLITECNICA
# PROGRAMA DE POS-GRADUACAO EM COMPUTACAO APLICADA - PPGCA
# MESTRADO EM COMPUTACAO APLICADA
# Disciplina: Modelagem estatistica
# Prof. Dr.: Rodrigo Sant'Ana
# Discentes: Andre Lucas Ribeiro, Leticia Zorzi Rama, Matheus Neis
# Itajai, Santa Catarina, Brasil

# ==========================================================================
# DESAFIO A
# ==========================================================================
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
# - Distribuicao Poisson: 
#   - Modelos lineares generalizados nulos e completos sem e com medidas de esforco amostral, com e sem offset
#   - Ajuste de modelos lineares generalizados exaurindo as combinacoes entre as variaveis preditoras com esforco amostral como offset
#   - Diagnostico e analise dos residuos 
#   - Verificacao da equidispersao
# - Distribuinao Binomial negativa
#   - Modelos lineares generalizados nulos e completo
#   - Ajuste de modelos lineares generalizados exaurindo as combinacoes entre as variaveis preditoras
#   - Diagnostico e analise dos residuos 
# - Parecer final (resposta as perguntas do desafio)
#
# # Base de dados: cafeicultura.csv
# ==========================================================================

# ==========================================================================
# CONFIGURAR AMBIENTE DE TRABALHO
# ==========================================================================

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
p_load(ggplot2, dplyr, MuMIn, ggcorrplot, patchwork, readxl, AER,
       DataExplorer, ggpubr, scatterplot3d, effects, car, hnp, statmod,
       datasets, DCluster, stargazer, olsrr, performance, report)

# Opcoes gerais
options(scipen = 10)

# Funcao auxiliar de impressao de texto
secao <- function(texto) {
barra <- paste(rep("=", 150), collapse = "")
  cat("\n", barra, "\n", texto, "\n", barra, "\n\n", sep = "")
}

# ==========================================================================
# BASE DE DADOS
# ==========================================================================
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

# ==========================================================================
# ANALISE EXPLORATORIA MULTIVARIADA 
# Variavel resposta com cada variavel candidata
# ==========================================================================
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

# ==========================================================================
# ANALISE EXPLORATORIA GLOBAL 
# Variavel resposta com todas as variaveis candidatas
# ==========================================================================
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

# ==========================================================================
# MODELOS LINEARES GENERALIZADOS - POISSON
# Modelos:
# - nulo
# - completo SEM medida de esforco amostral
# - completo COM medida de esforco amostral como OFFSET 
# - completo COM medida de esforco amostral como preditor livre
# ==========================================================================
secao("MODELOS LINEARES GENERALIZADOS - POISSON")

# n_brocas_capturadas e uma variavel de contagem o que significa dizer que e uma variavel discreta
# Poisson para modelar variaveis discretas

# ============
# Modelo nulo
# ============
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

# ====================
# Modelo completo
# sem esforco amostral
# ====================
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
anova(mod00, mod01, test = "Chisq")
residuos_mod01 <- ((deviance(mod00) - deviance(mod01)) / deviance(mod00)) * 100
#?anova

# Parecer
cat("\nEm comparação ao modelo nulo, o modelo completo SEM esforço amostral reduziu os resíduos em ",
    residuos_mod01, "%.")

# Analise da deviancia das covariaveis do modelo completo SEM esforco amostral
secao("Analise da deviância das covariáveis do modelo completo SEM medida de esforço amostral")
anova(mod01, test = "Chisq")

# Parecer 
cat("\nPara o modelo completo SEM esforço amostral, em termos de redução de resíduos, as variáveis com 
    redução significativa foram todas, exceto pH solo.")

# Tabela de estimacao dos parametros do modelo completo SEM esforco amostral
secao("Estimação dos parâmetros do modelo completo SEM medida de esforço amostral")
summary(mod01)

# Parecer 
cat("\nPara o modelo completo SEM esforço amostral, em termos de efeito, as variáveis com efeito 
    significativo foram quase todas, exceto declividade, pH solo e precipitação.")

# ====================
# Construcao
# do esforco amostral
# ====================

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

# ====================
# Modelo completo
# com esforco amostral
# passado como offset
# ====================
secao("Modelo completo COM medida de esforço amostral passada como OFFSET")
mod02 <- glm(
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
secao("Teste de verossimilhança entre os modelos nulo e completo COM medida de esforço amostral passado como OFFSET")
anova(mod00, mod02, test = "Chisq")
residuos_mod02 <- ((deviance(mod00) - deviance(mod02)) / deviance(mod00)) * 100

# Parecer
cat("\nEm comparação ao modelo nulo, o modelo completo COM esforço amostral passado como OFFSET reduziu 
    os resíduos em ", residuos_mod02, "%. A redução de resíduos no modelo completo sem medida de esforço 
    amostral foi de", residuos_mod01, "%. Portanto, a incorporação do offset no mecanismo do modelo
    linear generalizado resultou em uma diferença de", residuos_mod02 - residuos_mod01, "%.")

# Analise da deviancia das covariaveis do modelo completo COM esforco amostral passado como OFFSET
secao("Analise da deviância das covariáveis do modelo completo COM medida de esforço amostral")
anova(mod02, test = "Chisq")

# Parecer 
cat("\nPara o modelo completo COM esforço amostral passado como OFFSET, em termos de redução de resíduos, 
    as variáveis com redução significativa foram todas, exceto declividade e ph solo.")

# Tabela de estimacao dos parametros do modelo completo COM esforco amostral passado como OFFSET
secao("Estimação dos parâmetros do modelo completo COM medida de esforço amostral passada como OFFSET")
summary(mod02)

# Parecer 
cat("\nPara o modelo completo COM esforço amostral passado como OFFSET, em termos de efeito, as variáveis com 
    efeito significativo foram todas, exceto declividade, ph solo e precipitação.")

# ==================================
# Modelo completo
# com esforco amostral
# passado como preditor livre
# ==================================
secao("Modelo completo COM medida de esforço amostral passada como PREDITOR LIVRE")
mod03 <- glm(
  n_brocas_capturadas ~ adubacao_n_kg_ha
    + altitude_m
    + declividade_pct
    + densidade_plantio
    + idade_lavoura_anos
    + materia_organica_pct
    + ph_solo
    + precipitacao_safra_mm
    + umidade_relativa_pct
    + esforco_amostral, 
  family = poisson(link = "log"),
  data = caf1,
  na.action = "na.fail"
)

# Teste de verossimilhanca entre os modelos nulo e completo COM esforco amostral passado como PREDITOR LIVRE
secao("Teste de verossimilhança entre os modelos nulo e completo COM medida de esforço amostral passado como PREDITOR LIVRE")
anova(mod00, mod03, test = "Chisq")
residuos_mod03 <- ((deviance(mod00) - deviance(mod03)) / deviance(mod00)) * 100

# # Teste de verossimilhanca entre os modelos completo COM esforco amostral passado como PREDITOR LIVRE e como OFFSET
# secao("Teste de verossimilhança entre os modelos completos COM esforco amostral passado como PREDITOR LIVRE e como OFFSET")
# anova(mod02, mod03, test = "Chisq")
# residuos2_mod03 <- ((deviance(mod02) - deviance(mod03)) / deviance(mod02)) * 100

# Parecer
cat("\nO modelo completo COM esforço amostral passado como PREDITOR LIVRE aumentou os resíduos em ", residuos_mod02-residuos_mod03, "%.")

# Tabela de estimacao dos parametros do modelo completo COM esforco amostral passado como PREDITOR LIVRE
secao("Estimação dos parâmetros do modelo completo COM medida de esforço amostral passada como PREDITOR LIVRE")
summary(mod03)

# TO DO Parecer (comparacao do offset e do preditor livre)
cat("")

# ==========================================================================
# MODELOS LINEARES GENERALIZADOS - POISSON
# Exaurindo as combinacoes entre as variaveis preditoras 
# ==========================================================================
secao("MODELOS LINEARES GENERALIZADOS - POISSON - Exaurindo as combinações entre as variáveis preditoras")

# Ajustar diferentes modelos de forma iterativa
tab01 <- dredge(mod02, extra = "R^2")

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

# ==================================
# Modelo poisson final
# ==================================
secao("Modelo poisson final")
mod04 <- glm(
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

# OBS.: modelo poisson final (mod04) identifica-se com o modelo completo com medida de esforco amostral como offset (mod02)

# Teste de Verossimilhança entre os modelos nulo e modelo poisson final
secao("Teste de verossimilhança entre os modelos poisson nulo e final")
anova(mod00, mod04, test = "Chisq")
residuos_mod04 <- ((deviance(mod00) - deviance(mod04)) / deviance(mod00)) * 100

# Parecer
cat("\nEm comparação ao modelo nulo, o modelo poisson final reduziu os resíduos em ", residuos_mod04, "%.")

# Analise da deviancia das covariaveis do modelo poisson final 
secao("Analise da deviância das covariáveis do modelo poisson final")
anova(mod04, test = "Chisq")

# Parecer 
cat("\nPara o modelo poisson final, em termos de redução de resíduos, as variáveis com redução significativa foram:
    - adubação, 
    - altitude,
    - densidade de plantio,
    - idade da lavoura,
    - matéria orgânica,
    - precipitação, e
    - umidade."
  )

# Faixa de coeficientes possiveis para os parametros do modelo poisson final
secao("Faixa de coeficientes possíveis para os parâmetros do modelo poisson final")
confint(mod04)

# Tabela de estimacao dos parametros
secao("Estimação dos parâmetros do modelo poisson final")
summary(mod04)

# Visualizar efeitos
plot(effects::allEffects(mod04))

# Parecer 
cat("\nPara o modelo poisson final, em termos de efeito, as variáveis com efeito significativo foram:
  - adubação, 
  - altitude, 
  - densidade de plantio, 
  - idade da lavoura, 
  - matéria orgânica, e 
  - umidade.")

# ==========================================================================
# DIAGNOSTICO E ANALISE DE RESIDUOS
# ==========================================================================
secao("DIAGNÓSTICO DOS RESÍDUOS")

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

# ==========================================================================
# VERIFICACAO DA EQUIDISPERSAO DO MODELO POISSON
# ==========================================================================
secao("VERIFICAÇÃO DA EQUIDISPERSÃO")

#?dispersiontest

# Conforme a documentacao, poission assume que a esperanca condicional E[y] = μ 
# e a variancia VAR[y] = μ sao iguais

# ==================================
# Calcular a esperanca e a variancia 
# e observar se sao iguais
# ==================================

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

# ==================================
# Verificar a frequencia esperada
# e observada com um modelo poisson
# ==================================

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

# ==================================
# Teste formal de equidispersao
# ==================================

# Razao do desvio residual e dos graus de liberdade dos residuos
deviance(mod04)
df.residual(mod04)
deviance(mod04)/df.residual(mod04)

# Testar sobre a hipotese alternativa da variancia ser uma funcao linear
dispersiontest(mod04, trafo = 1)

# Testar sobre a hipotese alternativa da variancia ser uma funcao quadratica
dispersiontest(mod04, trafo = 2)

# TO DO Parecer
cat("")

# ==================================
# Verificar a frequencia esperada
# e observada com um modelo binomial
# negativo
# ==================================

# Parametro de superdispersao k
# calculado pelo metodo dos momentos
# um parametro que entra na binomial negativa
k <- (n_brocas_capt_media ^ 2) / (n_brocas_capt_variancia - n_brocas_capt_media)

# Frequencia esperada de acordo com um modelo binomial negativo
n_brocas_capt_freq_binom <- dnbinom(x = caf2$n_brocas_capturadas, size = k, mu = n_brocas_capt_media)

# Visualizar diferenca entre frequencia relativa esperada e observada no modelo binomial negativo
barplot(t(cbind(caf2$n_brocas_capt_freq_relat, n_brocas_capt_freq_binom)),
        names.arg = caf2$n_brocas_capt_freq_relat,
        ylim = c(0, 1), xlab = "Número de brocas capturadas",
        ylab = "Frequência relativa", beside = TRUE,
        legend.text = c("Observado", "Esperado"),
        main = "Binomial Negativa")

# Qui-quadrado para verificar se o modelo binomial negativo se ajusta aos dados
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

# ==========================================================================
# MODELOS LINEARES GENERALIZADOS - BINOMIAL NEGATIVA
# ==========================================================================
secao("MODELOS LINEARES GENERALIZADOS - BINOMIAL NEGATIVA")

# ============
# Modelo nulo
# ============
mod05 <- glm.nb(
  n_brocas_capturadas ~ 1,
  data = caf1,
  link = "log")

# ============
# Modelo completo
# ============
mod06 <- glm.nb(
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
secao("Teste de verossimilhança entre os modelos nulo e completo COM medida de esforço amostral passado como OFFSET")
anova(mod05, mod06, test = "Chisq")
residuos_mod06 <- ((deviance(mod05) - deviance(mod06)) / deviance(mod05)) * 100

# Parecer
cat("\nEm comparação ao modelo nulo, o modelo completo COM esforço amostral passado como OFFSET reduziu 
    os resíduos em ", residuos_mod06, "%.")

# Analise da deviancia das covariaveis do modelo completo COM esforco amostral passado como OFFSET
secao("Analise da deviância das covariáveis do modelo completo COM medida de esforço amostral")
anova(mod06, test = "Chisq")

# Parecer 
cat("\nPara o modelo completo COM esforço amostral passado como OFFSET, em termos de redução de resíduos, 
    as variáveis com redução significativa foram todas, exceto declividade, idade da lavoura e ph solo.")

# Tabela de estimacao dos parametros do modelo completo COM esforco amostral passado como OFFSET
secao("Estimação dos parâmetros do modelo completo COM medida de esforço amostral passada como OFFSET")
summary(mod06)

# Parecer 
cat("\nPara o modelo completo COM esforço amostral passado como OFFSET, em termos de efeito, as variáveis com 
    efeito significativo foram todas, exceto declividade, ph solo e precipitação.")

# Verificando a superdispersao ou sobredispersao
# Razao do desvio residual e dos graus de liberdade dos residuos
deviance(mod06)
df.residual(mod06)
deviance(mod06)/df.residual(mod06)

# ==========================================================================
# MODELOS LINEARES GENERALIZADOS - BINOMIAL NEGATIVA
# Exaurindo as combinacoes entre as variaveis preditoras 
# ==========================================================================
secao("MODELOS LINEARES GENERALIZADOS - BINOMIAL NEGATIVA - Exaurindo as combinações entre as variáveis preditoras")

# Ajustar diferentes modelos de forma iterativa
tab02 <- dredge(mod06, extra = "R^2")

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

# ==================================
# Modelo binomial negativa final
# ==================================
secao("Modelo binomial negativa final")
mod07 <- glm.nb(
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
  init.theta = 4.377708848, 
  link = log
)

# Teste de Verossimilhança entre os modelos binomial negativo nulo e modelo final
secao("Teste de verossimilhança entre os modelos binomial negativo nulo e final")
anova(mod05, mod07, test = "Chisq")
residuos_mod07 <- ((deviance(mod05) - deviance(mod07)) / deviance(mod05)) * 100

# Parecer
cat("\nEm comparação ao modelo nulo, o modelo binomial negativo reduziu os resíduos em ", residuos_mod07, "%.")

# Analise da deviancia das covariaveis do modelo binomial negativo final 
secao("Analise da deviância das covariáveis do modelo binomial negativo final")
anova(mod07, test = "Chisq")

# Parecer 
cat("\nPara o modelo binomial negativo final, em termos de redução de resíduos, as variáveis com redução significativa foram:
    - adubação, 
    - altitude,
    - densidade de plantio,
    - matéria orgânica,
    - precipitação, e
    - umidade."
  )

# Faixa de coeficientes possiveis para os parametros do modelo binomial negativo final
secao("Faixa de coeficientes possíveis para os parâmetros do modelo binomial negativo final")
confint(mod07)

# Tabela de estimacao dos parametros
secao("Estimação dos parâmetros do modelo binomial negativo final")
summary(mod07)

# Visualizar efeitos
plot(effects::allEffects(mod07))

# Parecer 
cat("\nPara o modelo binomial negativo final, em termos de efeito, as variáveis com efeito significativo foram:
  - adubação, 
  - altitude, 
  - densidade de plantio, 
  - idade da lavoura, 
  - matéria orgânica, e 
  - umidade.")

# ==========================================================================
# DIAGNOSTICO E ANALISE DE RESIDUOS
# ==========================================================================
secao("DIAGNÓSTICO DOS RESÍDUOS")

plot(mod07) 

# Envelope simulado dos residuos do modelo poisson final
res03 <- hnp(mod07, plot.sim = FALSE)

# Transformar em dataframe
res03 <- data.frame(x = res03$x, median = res03$median,
                    lower = res03$lower, upper = res03$upper,
                    residuals = res03$residuals)

# Grafico do envelope simulado dos residuos
p03 <- ggplot(data = res03) +
    geom_ribbon(aes(x = x, ymin = lower, ymax = upper), alpha = 0.8) +
    geom_line(aes(x = x, y = median), colour = "white") +
    geom_point(aes(x = x, y = residuals), pch = 21, fill = "white",
               colour = "black", size = 5, alpha = 0.5) +
    labs(x = "Quantis teóricos", y = "Resíduos") +
    theme_gray(base_size = 18)
p03

# Exportar os residuos do modelo
res04 <- fortify(mod07)
res04$ID <- 1:nrow(res04)

# Grafico de dispersao dos residuos
p04 <- ggplot(data = res04, aes(x = ID, y = .stdresid)) +
    geom_point(pch = 21, fill = "white", colour = "black", size = 5,
               alpha = 0.8) +
    geom_hline(yintercept = 0, colour = "red") +
    labs(x = "Índice da amostra", y = "Resíduos padronizados") +
    theme_gray(base_size = 18)
p04

# Histograma dos residuos
p05 <- ggplot(data = res04, aes(x = .stdresid)) +
    geom_histogram(binwidth = 1, boundary = 1, closed = "right",
                   fill = "white", colour = "black") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 25)) +
    labs(x = "Resíduos padronizados", y = "Frequência") +
    theme_gray(base_size = 18)
p05 # TO DO Verificar 

# Visualizar os graficos lado a lado
(p05 | p04 | p03)

# ==========================================================================
# PARECER FINAL
# ==========================================================================

# TO DO 
# Pergunta do corpo técnico:
# Quais condições de talhão favorecem a infestação pela broca, 
# e quanto se ganha em pressão de praga ao alterar cada uma delas?

# ==========================================================================
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
# ==========================================================================