# =====================================================================
# UNIVERSIDADE DO VALE DO ITAJAI - UNIVALI
# ESCOLA POLITECNICA
# PROGRAMA DE POS-GRADUACAO EM COMPUTACAO APLICADA - PPGCA
# MESTRADO EM COMPUTACAO APLICADA
# Disciplina: Modelagem estatistica
# Prof. Dr.: Rodrigo Sant'Ana
# Discentes: Andre Lucas Ribeiro, Leticia Zorzi Rama, Matheus Neis
# Itajai, Santa Catarina, Brasil

# =====================================================================
# DESAFIO A
# =====================================================================
# Descricao:
# Modelo Poisson (ou extensão para contagens)
# Quais condições de talhão favorecem a infestação pela broca, 
# e quanto se ganha em pressão de praga ao alterar cada uma delas?
#
# Variavel resposta:
# `n_brocas_capturadas`
# 
# Variaveis candidatas:
# `altitude_m`, `declividade_pct`, `idade_lavoura_anos`, `densidade_plantio`, `adubacao_n_kg_ha`
# `precipitacao_safra_mm`, `umidade_relativa_pct`, `ph_solo`, `materia_organica_pct`
#
# O problema oculto do esforço amostral: 
# o número de brocas capturadas depende diretamente de quantas armadilhas foram 
# instaladas e por quantos dias ficaram expostas. Um talhão com 12 armadilhas por 
# 30 dias tem seis vezes mais esforço de captura que um com 4 armadilhas por 15 dias. 
# Comparar contagens brutas entre talhões é, portanto, comparar coisas diferentes.
#
# 1. Construam a medida de esforço amostral adequada a partir de `n_armadilhas` e
#    `dias_exposicao` e incorporem-na ao modelo pelo mecanismo correto do MLG.
#    Justifiquem por que esse termo entra com coeficiente fixado em 1 e não
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
# - Ajuste de modelos lineares generalizados nulos e completos sem e com medidas de esforco amostral
# - Ajuste de modelos lineares generalizados exaurindo as combinacoes entre as variaveis preditoras com esforco amostral
# - Analise dos residuos 
# - Parecer final (resposta às perguntas do desafio)
# 
# Base de dados: cafeicultura.csv
# =====================================================================

# =====================================================================
# CONFIGURAR AMBIENTE DE TRABALHO
# =====================================================================

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
barra <- paste(rep("=", 80), collapse = "")
  cat("\n", barra, "\n", texto, "\n", barra, "\n\n", sep = "")
}

# =====================================================================
# BASE DE DADOS
# =====================================================================

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

# =====================================================================
# ANALISE EXPLORATORIA MULTIVARIADA 
# Variavel resposta com cada variavel candidata
# =====================================================================

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
cat("\nOs resultados das correlações indicam uma associação moderada positiva do número de 
    brocas capturadas com a umidade relativa (%). As demais variáveis possuem associação fraca 
    com a variável resposta.")

# =====================================================================
# ANALISE EXPLORATORIA GLOBAL 
# Variavel resposta com todas as variaveis candidatas
# =====================================================================

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
cat("\nA associação moderada e positiva do número de brocas capturadas com a umidade relativa demonstrada
    nas análises individuais se demonstrou significativa a partir do teste de hipóteses. Das demais associações 
    avaliadas como fracas, aquelas significativas são precipitação, matéria orgânica e densidade de plantio com 
    direção positiva e adubação e altitude com direção negativa. A matriz de correlação também revelou uma 
    associação moderada entre umidade e precipitação.")

# =====================================================================
# MODELOS LINEARES GENERALIZADOS
# Modelos nulo e completos sem e com medida esforco amostral
# a partir de 'n_armadilhas' e 'dias_exposicao'
# =====================================================================

# n_brocas_capturadas e uma variavel de contagem o que significa dizer que e uma variavel discreta
# Poisson para modelar variaveis discretas

# Relembrando o problema
# 1. Construam a medida de esforço amostral adequada a partir de `n_armadilhas` e
#    `dias_exposicao` e incorporem-na ao modelo pelo mecanismo correto do MLG.
#    Justifiquem por que esse termo entra com coeficiente fixado em 1 e não
#    como preditor livre. (Verifiquem essa afirmação empiricamente: ajustem
#    também o modelo com o log do esforço como preditor livre e comparem o
#    coeficiente estimado com 1.)

# Criar medida de esforco amostral
caf <- caf %>%
  mutate(
    esforco_amostral = n_armadilhas * dias_exposicao
  )

# Adicionar a medida de esforco amostral na base de dados que esta sendo manipulada
caf1$esforco_amostral <- caf$esforco_amostral

# Calcular a frequencia esperada para um modelo poisson
# Frequencia relativa
caf1$n_brocas_capturadas_freq_relativa <- caf1$n_brocas_capturadas/sum(caf1$n_brocas_capturadas)
# Media de brocas capturadas
n_brocas_capturadas_media <- sum(caf1$n_brocas_capturadas * caf1$n_brocas_capturadas_freq_relativa)
# Frequencia esperada para um modelo poisson
n_brocas_capturadas_freq_esp_pois <- dpois(x = caf1$n_brocas_capturadas, lambda = n_brocas_capturadas_media)
print(n_brocas_capturadas_freq_esp_pois)

# MODELO nulo
# Para este primeiro caso e para fins de aprendizagem, o significado
# de cada argumento sera comentado conforme constante em ?glm
?glm

mod00 <- glm(
  formula = n_brocas_capturadas ~ 1, # descricao simbolica do modelo, a funcao a ser modelada
  family = poisson(link = "log"), # descricao da distribuicao do erro e da funcao de ligacao 
  data = caf1, # base de dados com as variaveis de interesse
  weights = NULL, # pesos para as variaveis preditoras
  # subset, # subset de observacoes a ser utilizadas no ajuste do modelo
  na.action = "na.fail", # indica o que deve acontecer se a base de dados conter dados ausentes
  # start, # valores de partida para os parametros do preditor linear
  # etastart, # valores de partida para o preditor (qual a diferenca deste com o start?)
  offset = NULL, # especifica um componente a priori para ser incluido no preditor linear durante o ajuste
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

# TO DO continuar lendo a documentacao de ?glm

# MODELO completo SEM esforco amostral 
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

# Teste de Verossimilhança entre os modelos nulo e completo SEM esforco amostral
#?anova
anova(mod00, mod01, test = "Chisq")

# Parecer
cat("\nEm comparação ao modelo nulo, o modelo completo SEM esforço amostral reduziu os resíduos em ",
    ((2108.2 - 1479.1) / 2108.2) * 100, "%.")

# Analise da deviancia das covariáveis do modelo completo SEM esforco amostral
anova(mod01, test = "Chisq")

# Parecer 
cat("\nPara o modelo completo SEM esforço amostral, em termos de redução de resíduos, as variáveis com 
    redução significativa foram todas, exceto pH_solo.")

# Tabela de estimacao dos parametros
summary(mod01)

# Parecer 
cat("\nPara o modelo completo SEM esforço amostral, em termos de efeito, as variáveis com efeito 
    significativo foram quase todas, exceto declividade_pct, pH_solo e precipitacao_safra_mm.")

# MODELO completo COM esforco amostral passado como OFFSET 
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

# Teste de Verossimilhança entre os modelos nulo e completo COM esforco amostral passado como OFFSET
#?anova
anova(mod00, mod02, test = "Chisq")

# Parecer
cat("\nEm comparação ao modelo nulo, o modelo completo COM esforço amostral passado como OFFSET reduziu 
    os resíduos em ", ((2108.2 - 1008.5) / 2108.2) * 100, "%.")

# Analise da deviancia das covariáveis do modelo completo COM esforco amostral passado como OFFSET
anova(mod02, test = "Chisq")

# Parecer 
cat("\nPara o modelo completo COM esforço amostral passado como OFFSET, em termos de redução de resíduos, 
    as variáveis com redução significativa foram todas, exceto declividade_pct e ph_solo.")

# Tabela de estimacao dos parametros
summary(mod02)

# Parecer 
cat("\nPara o modelo completo COM esforço amostral passado como OFFSET, em termos de efeito, as variáveis com 
    efeito significativo foram todas, exceto declividade_pct, ph_solo e precipitacao_safra_mm.")

# TO DO modelo COM esforco amostral SEM offset

# TO DO
# 2. Verifiquem a equidispersão. Reportem a razão de Pearson sobre os graus de
#    liberdade e conduzam um teste formal de superdispersão.
# 3. Caso haja superdispersão, comparem quasi-Poisson e binomial negativa,
#    explicando como cada uma modela a variância e por que isso importa aqui.
# 4. Selecionem variáveis e interpretem os efeitos como razões de taxas de captura
#    por unidade de esforço.

# =====================================================================
# MODELOS LINEARES GENERALIZADOS
# Exaurindo as combinacoes entre as variaveis preditoras 
# =====================================================================

# Ajustar diferentes modelos de forma iterativa
tab01 <- dredge(mod02, extra = "R^2")

# Quantidade de modelos ajustados
dim(tab01)

# Identificar o melhor modelo com o criterio de informacao de akaike (AICc)
filter(tab01, AICc == min(AICc))

# Parecer
cat("\nConforme o AICc, o modelo que melhor se ajusta aos dados tem adubacao_n_kg_ha, altitude_m,  
    declividade_pct, densidade_plantio, idade_lavoura_anos, materia_organica_pct, ph_solo, 
    precipitacao_safra_mm e umidade_relativa_pct como variáveis preditoras.")

# MODELO final
mod <- glm(formula = n_brocas_capturadas ~ adubacao_n_kg_ha 
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

# Teste de Verossimilhança entre os modelos nulo e modelo final
#?anova
anova(mod00, mod, test = "Chisq")

# Parecer
cat("\nEm comparação ao modelo nulo, o modelo final reduziu os resíduos em ", 
    ((2108.2 - 1008.5) / 2108.2) * 100, "%.")

# Analise da deviancia das covariáveis do modelo final 
anova(mod, test = "Chisq")

# Parecer 
cat("\nPara o modelo final, em termos de redução de resíduos, as variáveis com redução 
    significativa foram adubacao_n_kg_ha, altitude_m, densidade_plantio, idade_lavoura_anos, 
    materia_organica_pct, precipitacao_safra_mm e umidade_relativa_pct.")

# Tabela de estimacao dos parametros
summary(mod)

# Parecer 
cat("\nPara o modelo completo com esforço amostral, em termos de efeito, as variáveis com efeito 
    significativo foram adubacao_n_kg_ha, altitude_m, densidade_plantio, idade_lavoura_anos, 
    materia_organica_pct e umidade_relativa_pct.")

# =====================================================================
# ANALISE DE RESIDUOS
# =====================================================================

# Envelope simulado dos residuos do modelo final
res01 <- hnp(mod, plot.sim = FALSE)

#####@> Transformando em data.frame...
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
res02 <- fortify(mod)
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
p02

# Visualizar os graficos lado a lado
(p01 | p02 | p00)

# =====================================================================
# PARECER FINAL
# =====================================================================

# Pergunta do corpo técnico:
# Quais condições de talhão favorecem a infestação pela broca, 
# e quanto se ganha em pressão de praga ao alterar cada uma delas?

# =====================================================================
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
# =====================================================================