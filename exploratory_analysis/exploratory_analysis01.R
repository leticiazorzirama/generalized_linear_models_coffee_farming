# ===============================
# CONFIGURAR AMBIENTE DE TRABALHO
# ===============================

# Conferir caminhos
getwd()
setwd("inserir caminho")
getwd()

# Pacotes
library(spdep) ## I de Moran, LISA e matriz de pesos
library(energy) ## Correlacao de distancia (dCor)
library(boot) ## Bootstrap dos intervalos de confianca
library(DescTools) ## V de Cramer com correcao de vies
library(psych) ## Correlacao parcial
library(car) ## Fator de inflacao de variancia (VIF)
library(moments) ## Assimetria e curtose
library(dplyr) ## Manipulacao de dados
library(tidyr) ## Reorganizacao de dados
library(ggplot2) ## Graficos
library(GGally) ## Matriz de dispersao
library(corrplot) ## Mapa de calor de matrizes de correlacao
library(patchwork) ## Composicao de graficos
library(stringr) ## Manipulacao de caracteres

# Opcoes gerais
options(width = 80, scipen = 6, digits = 4)

# Semente de reprodutibilidade
set.seed(20260818)

# Tema para ggplot
rgb01 <- "#303033"
seta <- grid::arrow(length = grid::unit(0.2, "cm"), type = "open")
my_theme <- function(base_size = 14, base_family = "Helvetica") {
theme_bw(base_size = base_size, base_family = base_family) %+replace%
theme(axis.ticks = element_blank(),
axis.line = element_line(arrow = seta, colour = rgb01),
legend.background = element_blank(),
legend.key = element_blank(),
panel.background =
element_rect(fill = ggplot2::alpha(rgb01, 0.05),
colour = "white"),
panel.border = element_blank(),
panel.grid = element_line(linetype = "solid",
linewidth = 0.2,
colour = "white"),
strip.text = element_text(colour = "white",
margin = margin(0.3, 0.3, 0.3, 0.3,
"cm"),
face = "bold"),
strip.background = element_rect(fill = rgb01,
colour = rgb01),
plot.background = element_blank(),
complete = TRUE)
}

# Definir o tema padrao
theme_set(my_theme())

# Funcao auxiliar
secao <- function(texto) {
barra <- paste(rep("=", 70), collapse = "")
cat("\n", barra, "\n", texto, "\n", barra, "\n\n", sep = "")
}

# Teste de Fischer
ic_fisher <- function(r, n, conf = 0.95) {
z <- atanh(r)
ep <- 1 / sqrt(n - 3)
q <- qnorm(1 - (1 - conf) / 2)
data.frame(r = r,
n = n,
inferior = tanh(z - q * ep),
superior = tanh(z + q * ep))
}

# ==============
# BASE DE DADOS
# ==============

# Carregar a base de dados
caf <- read.csv("data/cafeicultura.csv", sep=",")

# Inspecao da base de dados
secao("Inspeção da base de dados")

# Estrutura geral
cat("\nNúmero de observações (linhas):", nrow(caf))
cat("\nNúmero de variáveis (colunas):", ncol(caf), "\n")
cat("\nNome das variáveis:")
print(names(caf))
cat("\nEstrutura geral:\n")
print(str(caf))

# Cabecalho, primeiras 10 linhas
cat("\nCabeçalho (primeiras 10 linhas):\n")
print(head(caf, 10))

# Resumo estatistico preliminar
cat("\nResumo estatístico preliminar:\n")
print(summary(caf))
# Possui de dados faltantes? Nao
# Suspeita de valores extremos? A umidade relativa em 92.2% pode ser um valor extremo
# Suspeita de valores impossiveis? Nao

# Diagnostico global de dados faltantes
total_celulas <- prod(dim(caf))
total_na <- sum(is.na(caf))
pct_geral_na <- round(total_na / total_celulas * 100, 1)

cat("\nDiagnóstico global de dados faltantes\n")
cat("Total de células na base de dados:", total_celulas, "\n")
cat("Total de dados faltantes (NA):", total_na, "\n")
cat("Porcentagem geral de dados faltantes:", pct_geral_na, "%\n\n")

# =====================
# ANALISES UNIVARIADAS
# =====================

# Variavel 'id_talhao'
secao("Análise univariada - variável: id_talhao")

# Verificar se cada observacao corresponde a um talhao
length(unique(caf$id_talhao))
# Conclusao: cada observacao corresponde a um talhao
cat("\nA base de dados contém", length(unique(caf$id_talhao)), "observações,",
"cada uma correspondente a um talhão.")

# Variavel 'regiao_produtora'
secao("Análise univariada - variável: regiao_produtora")
table(caf$regiao_produtora)
names(table(caf$regiao_produtora))

# Moda
regiao_produtora_moda <- names(which.max(table(caf$regiao_produtora)))

# Ranking
regiao_produtora_rank <- caf %>%
  group_by(regiao_produtora) %>%
  summarise(frequencia = n()) %>%
  arrange(desc(frequencia))

# Grafico de barras
# TO DO

# Exibir resultados preliminares
# Moda
cat("\nDas", dim(table(caf$regiao_produtora)), "regiões produtoras,",
"a moda é:", str_to_title(regiao_produtora_moda))

# Ranking
cat("\nO ranking das regiões produtoras é:")
print(regiao_produtora_rank)

# Variavel 'cultivar'
secao("Análise univariada - variável: cultivar")
table(caf$cultivar)
names(table(caf$cultivar))

# Moda
cultivar_moda <- names(which.max(table(caf$cultivar)))

# Ranking
cultivar_rank <- caf %>%
  group_by(cultivar) %>%
  summarise(frequencia = n()) %>%
  arrange(desc(frequencia))

# Grafico de barras
# TO DO

# Exibir resultados preliminares
# Moda
cat("\nDos", dim(table(caf$cultivar)), "cultivares,",
"a moda é:", str_to_title(cultivar_moda))

# Ranking
cat("\nO ranking dos cultivares é:")
print(cultivar_rank)

# Variavel 'manejo'
secao("Análise univariada - variável: manejo")
table(caf$manejo)
names(table(caf$manejo))

# Moda
manejo_moda <- names(which.max(table(caf$manejo)))

# Ranking
manejo_rank <- caf %>%
  group_by(manejo) %>%
  summarise(frequencia = n()) %>%
  arrange(desc(frequencia))

# Grafico de barras
# TO DO

# Exibir resultados preliminares
# Moda
cat("\nDos", dim(table(caf$manejo)), "manejos,",
"a moda é:", str_to_title(manejo_moda))

# Ranking
cat("\nO ranking dos manejos é:")
print(manejo_rank)

# Variavel 'irrigacao'
secao("Análise univariada - variável: irrigação")
table(caf$irrigacao)
names(table(caf$irrigacao))

# Moda
irrigacao_moda <- names(which.max(table(caf$irrigacao)))

# Ranking
irrigacao_rank <- caf %>%
  group_by(irrigacao) %>%
  summarise(frequencia = n()) %>%
  arrange(desc(frequencia))

# Grafico de barras
# TO DO

# Exibir resultados preliminares
# Moda
cat("\nA irrigação predomina?", str_to_title(irrigacao_moda))

# Ranking
cat("\nO ranking da irrigação é:")
print(irrigacao_rank)

# Variavel 'altitude'
secao("Análise univariada - variável: altitude")

# Altitude do parque cafeeiro de Minas gerais, segundo Bernardes et al. (2012)
# O parque cafeeiro esta distribuido em altitudes variando entre 500 e 1.200 m.  

# Resumo estatistico
summary(caf$altitude_m)

# Histograma
altitude_m_plot <- ggplot(caf, aes(altitude_m)) +
  geom_histogram(binwidth = 50) +
  theme_minimal()

print(altitude_m_plot)

# Resultados
# Extrair estatisticas descritivas
alt_min <- min(caf$altitude_m)
alt_max <- max(caf$altitude_m)
alt_mediana <- median(caf$altitude_m)
alt_media <- mean(caf$altitude_m)
alt_sd <- sd(caf$altitude_m)
alt_q1 <- quantile(caf$altitude_m, 0.25)
alt_q3 <- quantile(caf$altitude_m, 0.75)

# Exibir resultados preliminares
cat("\nAltitude mínima:", alt_min, "metros\n")
cat("nAltitude máxima:", alt_max, "metros\n")
cat("1º quartil:", alt_q1, "metros\n")
cat("Mediana:", alt_mediana, "metros\n")
cat("Média:", round(alt_media, 2), "metros\n")
cat("3º quartil:", alt_q3, "metros\n")
cat("Desvio padrão:", round(alt_sd, 2), "metros\n")

# Parecer
cat("\nA altitude tem uma distribuição normal,
sem a presença de valores extremos ou impossíveis. Segundo Bernard et al. (2012),
a altitude no parque cafeeiro em Minas gerais varia entre 500 a 1.200 metros. 
Os dados indicam uma frequência maior de altitudes entre 830 e 1046 metros, com altitudes mínimas
de 496 e máximas de 1343 metros.")

# Variavel 'declividade_pct'
secao("Análise univariada - variável: declividade")

# Classes de declividade de Minas Gerais, segundo INPE (2023)
# https://idesisema.meioambiente.mg.gov.br/geonetwork/srv/api/records/46b18109-06ab-4529-9401-5a662f1e39b2
# Plano: ate 3%
# Suave-ondulado: 3-8%
# Ondulado: 8-20%
# Forte-ondulado: 20-45%
# Montanhoso: 45-75%
# Escarpado: >75%

# Declividade do parque cafeeiro de Minas Gerais, segundo Bernard et al. (2012)
# Sao encontradas lavouras em praticamente todas as faixas de declividade,
# porem ha um predominio de lavouras em declividades entre 5 e 15%.

# Resumo estatistico
summary(caf$declividade_pct)

# Histograma
declividade_pct_plot <- ggplot(caf, aes(declividade_pct)) +
  geom_histogram(binwidth = 3) +
  theme_minimal()

print(declividade_pct_plot)

# Resultados
# Extrair estatisticas descritivas
decl_min <- min(caf$declividade_pct)
decl_max <- max(caf$declividade_pct)
decl_mediana <- median(caf$declividade_pct)
decl_media <- mean(caf$declividade_pct)
decl_sd <- sd(caf$declividade_pct)
decl_q1 <- quantile(caf$declividade_pct, 0.25)
decl_q3 <- quantile(caf$declividade_pct, 0.75)

# Exibir resultados preliminares
cat("\nDeclividade mínima:", decl_min)
cat("Declividade máxima:", decl_max)
cat("1º quartil:", decl_q1)
cat("Mediana:", decl_mediana)
cat("Média:", round(decl_media, 2))
cat("3º quartil:", decl_q3)
cat("Desvio padrão:", round(decl_sd, 2))

# Parecer
cat("\nA declividade tem uma distribuição normal,
com a predominância de terrenos forte-ondulados, segundo
a classificação do INPE (2023). Há, entretanto, uma
heterogeneidade de terrenos visto que a declividade varia de
0% a 51.4%, caracterizando a ocorrência de terrenos planos à
montanhosos, respectivamente. Os resultados diferem da caracterização de 
Bernard et al. (2012) cuja declividae do parque cafeeiro em Minas Gerais
predomina entre 5 e 15%.")

# Variavel 'idade_lavoura_anos'
secao("Análise univariada - variável: idade_lavoura_anos")

# Resumo estatistico
summary(caf$idade_lavoura_anos)

# Histograma
idade_lavoura_anos_plot <- ggplot(caf, aes(idade_lavoura_anos)) +
  geom_histogram(binwidth = 2) +
  theme_minimal()

print(idade_lavoura_anos_plot)

# Resultados
# Extrair estatisticas descritivas
idade_min <- min(caf$idade_lavoura_anos)
idade_max <- max(caf$idade_lavoura_anos)
idade_mediana <- median(caf$idade_lavoura_anos)
idade_media <- mean(caf$idade_lavoura_anos)
idade_sd <- sd(caf$idade_lavoura_anos)
idade_q1 <- quantile(caf$idade_lavoura_anos, 0.25)
idade_q3 <- quantile(caf$idade_lavoura_anos, 0.75)

# Exibir resultados preliminares
cat("\nIdade mínima:", idade_min)
cat("Idade máxima:", idade_max)
cat("1º quartil:", idade_q1)
cat("Mediana:", idade_mediana)
cat("Média:", round(idade_media, 2))
cat("3º quartil:", idade_q3)
cat("Desvio padrão:", round(idade_sd, 2))

# Parecer
# TO DO 

# ======================
# ANALISES MULTIVARIADAS
# ======================

# TO DO

# ================
# ANALISES GLOBAIS
# ================

# Variaveis numericas
num_caf <- caf %>%
  select(where(is.numeric))
num_caf <- names(num_caf)
 
secao("Forma das distibuições - Cafeicultura")

diagnostico_forma <- data.frame(
variavel = num_caf,
media = sapply(caf[num_caf], mean, na.rm = TRUE),
mediana = sapply(caf[num_caf], median, na.rm = TRUE),
assimetria = sapply(caf[num_caf], skewness, na.rm = TRUE),
curtose = sapply(caf[num_caf], kurtosis, na.rm = TRUE),
shapiro_p = sapply(caf[num_caf],
function(x) shapiro.test(x)$p.value))
print(diagnostico_forma, row.names = FALSE)

# Matriz de dispersao

p_pares <- ggpairs(
caf[, num_caf],
lower = list(
continuous = wrap("points", alpha = 0.5, size = 1)
),
upper = list(
continuous = wrap("cor", method = "pearson", size = 4)
),
diag = list(
continuous = wrap("densityDiag", alpha = 0.5)
),
title = "Matriz de dispersão - Pearson - Cafeicultura"
)

p_pares

ggsave(
"matriz_dispersao_pearson_cafeicultura.png",
p_pares,
width = 22,
height = 22,
units = "in",
dpi = 300
)