# =====================================================================
# ANALISE ESTATISTICA DESCRITIVA DA CAFEICULTURA
# =====================================================================
#
# UNIVERSIDADE DO VALE DO ITAJAI - UNIVALI
# ESCOLA POLITECNICA
# PROGRAMA DE POS-GRADUACAO EM COMPUTACAO APLICADA - PPGCA
# MESTRADO EM COMPUTACAO APLICADA
# Disciplina: Modelagem estatistica
# Prof. Dr.: Rodrigo Sant'Ana
# Discentes: Andre Lucas Ribeiro, Leticia Zorzi Rama, Matheus Neis
# Itajai, Santa Catarina, Brasil
#
# Descrição:
# Este script realiza a inspecao e a analise estatistica descritiva 
# de uma base de dados referente a cafeicultura.
# A analise estatistica dos dados serve como etapa previa para posteriores
# analises inferenciais e modelagens dos dados conforme perguntas especificas 
# pre-formuladas e propostas ao grupo de trabalho para responde-las.
#
# O script esta estruturado em:
# - Configurar ambiente de trabalho
# - Carregamento e inspecao inicial da base de dados
# - Analise univariada por ordem das colunas da base de dados
# - Analises multivariadas
# - Analises globais
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

# =====================================================================
# BASE DE DADOS
# =====================================================================

# Carregar a base de dados
caf <- read.csv("cafeicultura.csv", sep=",")

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

# =====================================================================
# ANALISES UNIVARIADAS
# =====================================================================

# VARIAVEL 'id_talhao'
secao("Análise univariada - variável: id_talhao")

# O que e um talhao?
# (...) o talhao refere-se a determinada área da lavoura cafeeira que pode ser considerada homogênea, 
# por ter a mesma variedade de café plantada na mesma data, em um solo de características físicas, 
# químicas e topográficas semelhantes, além de receber o mesmo manejo agronômico e o mesmo tratamento 
# administrativo (Santos et al., 2009).

# Verificar se cada observacao corresponde a um talhao
length(unique(caf$id_talhao))
# Conclusao: cada observacao corresponde a um talhao
cat("\nA base de dados contém", length(unique(caf$id_talhao)), "observações,",
"cada uma correspondente a um talhão.")

# VARIAVEL 'regiao_produtora'
secao("Análise univariada - variável: regiao_produtora")
table(caf$regiao_produtora)
names(table(caf$regiao_produtora))

# Moda
regiao_produtora_moda <- names(which.max(table(caf$regiao_produtora)))

# Ranking absoluto
regiao_produtora_rank <- caf %>%
  group_by(regiao_produtora) %>%
  summarise(frequencia = n()) %>%
  arrange(desc(frequencia))

# Grafico de barras
regiao_produtora_plot <- ggplot(
  regiao_produtora_rank,
  aes(x = regiao_produtora, y = frequencia)
) +
  geom_col() +
  labs(
    x = "Região produtora",
    y = "Frequência absoluta"
  ) +
  coord_flip() +
  theme_minimal()

print(regiao_produtora_plot)

# Exibir resultados
# Moda
cat("\nDas", dim(table(caf$regiao_produtora)), "regiões produtoras,",
"a moda é:", str_to_title(regiao_produtora_moda))

# Ranking
cat("\nO ranking das regiões produtoras é:")
print(regiao_produtora_rank)

# VARIAVEL 'cultivar'
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
cultivar_plot <- ggplot(
  cultivar_rank,
  aes(
    x = reorder(cultivar, frequencia),
    y = frequencia
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Cultivar",
    y = "Frequência absoluta"
  ) +
  theme_minimal()

print(cultivar_plot)

# Exibir resultados
# Moda
cat("\nDos", dim(table(caf$cultivar)), "cultivares,",
"a moda é:", str_to_title(cultivar_moda))

# Ranking
cat("\nO ranking dos cultivares é:")
print(cultivar_rank)

# VARIAVEL 'manejo'
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
manejo_plot <- ggplot(
  manejo_rank, 
  aes(
    x = reorder(manejo, frequencia),
    y = frequencia
    )
) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Manejo",
    y = "Frequência absoluta"
  ) +
  theme_minimal()

print(manejo_plot)

# Exibir resultados
# Moda
cat("\nDos", dim(table(caf$manejo)), "manejos,",
"a moda é:", str_to_title(manejo_moda))

# Ranking
cat("\nO ranking dos manejos é:")
print(manejo_rank)

# VARIAVEL 'irrigacao'
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
irrigacao_plot <- ggplot(
  irrigacao_rank, 
  aes(
    x = reorder(irrigacao, frequencia),
    y = frequencia
    )
) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Irrigação",
    y = "Frequência absoluta"
  ) +
  theme_minimal()

print(irrigacao_plot)

# Exibir resultados
# Moda
cat("\nA irrigação predomina?", str_to_title(irrigacao_moda))

# Ranking
cat("\nO ranking da irrigação é:")
print(irrigacao_rank)

# VARIAVEL 'altitude'
secao("Análise univariada - variável: altitude")

# Altitude do parque cafeeiro de Minas gerais, segundo Bernardes et al. (2012)
# O parque cafeeiro esta distribuido em altitudes variando entre 500 e 1.200 m.  

# Resumo estatistico
summary(caf$altitude_m)

# Extrair estatisticas descritivas
alt_min <- min(caf$altitude_m)
alt_max <- max(caf$altitude_m)
alt_mediana <- median(caf$altitude_m)
alt_media <- mean(caf$altitude_m)
alt_sd <- sd(caf$altitude_m)
alt_q1 <- quantile(caf$altitude_m, 0.25)
alt_q3 <- quantile(caf$altitude_m, 0.75)

# Histograma
altitude_m_plot <- ggplot(caf, aes(x = altitude_m)) +
  geom_histogram(binwidth = 50, colour = "white") +
  geom_vline(
  xintercept = alt_media,
  linetype = "dashed"
) +
annotate(
  "text",
  x = alt_media,
  y = Inf,
  label = paste0("Média = ", round(alt_media, 2), " m"),
  vjust = 1.5,
  hjust = -0.05
) +
  labs(
    x = "Altitude (m)",
    y = "Frequência absoluta"
  ) +
  theme_minimal()

print(altitude_m_plot)

# Resultados
# Exibir resultados
cat("\nAltitude mínima:", alt_min, "metros\n")
cat("nAltitude máxima:", alt_max, "metros\n")
cat("1º quartil:", alt_q1, "metros\n")
cat("Mediana:", alt_mediana, "metros\n")
cat("Média:", round(alt_media, 2), "metros\n")
cat("3º quartil:", alt_q3, "metros\n")
cat("Desvio padrão:", round(alt_sd, 2), "metros\n")

# Parecer
cat("\nA altitude tem uma distribuição normal, sem a presença de valores extremos ou impossíveis. 
Segundo Bernard et al. (2012), a altitude no parque cafeeiro em Minas gerais varia entre 500 a 1.200 
metros. Os dados indicam uma frequência maior de altitudes entre 830 e 1046 metros, com altitudes 
mínimas de 496 e máximas de 1343 metros.")

# VARIAVEL 'declividade_pct'
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

# Extrair estatisticas descritivas
decl_min <- min(caf$declividade_pct)
decl_max <- max(caf$declividade_pct)
decl_mediana <- median(caf$declividade_pct)
decl_media <- mean(caf$declividade_pct)
decl_sd <- sd(caf$declividade_pct)
decl_q1 <- quantile(caf$declividade_pct, 0.25)
decl_q3 <- quantile(caf$declividade_pct, 0.75)

# Histograma
declividade_pct_plot <- ggplot(caf, aes(x = declividade_pct)) +
  geom_histogram(binwidth = 3, colour = "white") +
  geom_vline(
  xintercept = c(3, 8, 20, 45, 75),
  linetype = "dashed"
) +
annotate(
  "text",
  x = c(3, 8, 20, 45, 75),
  y = Inf,
  label = c("3%", "8%", "20%", "45%", "75%"),
  vjust = 1.5,
  hjust = -0.1
) +
  labs(
    x = "Declividade (%)",
    y = "Frequência absoluta"
  ) +
  theme_minimal()

print(declividade_pct_plot)

# Resultados
# Exibir resultados
cat("\nDeclividade mínima:", decl_min)
cat("Declividade máxima:", decl_max)
cat("1º quartil:", decl_q1)
cat("Mediana:", decl_mediana)
cat("Média:", round(decl_media, 2))
cat("3º quartil:", decl_q3)
cat("Desvio padrão:", round(decl_sd, 2))

# Parecer
cat("\nA declividade tem uma distribuição normal, levemente assimetrica à direita,
com a predominância de terrenos forte-ondulados, segundo a classificação do INPE (2023). 
Há, entretanto, uma heterogeneidade de terrenos visto que a declividade varia de
0% a 51.4%, caracterizando a ocorrência de terrenos planos à montanhosos, respectivamente. 
Os resultados diferem da caracterização de  Bernard et al. (2012) cuja declividae do parque 
cafeeiro em Minas Gerais predomina entre 5 e 15%.")

# VARIAVEL 'idade_lavoura_anos'
secao("Análise univariada - variável: idade_lavoura_anos")

# Idade da lavoura
# Até  os  dois  anos, chamada de lavoura jovem e após os dois anos, 
# conhecida como manejo em café adulto (Nascimento et al., 2020).
# Uma vez plantada, a lavoura permanece por 10, 20, 30 anos ou mais (Silva et al., 2002).

# Resumo estatistico
summary(caf$idade_lavoura_anos)

# Extrair estatisticas descritivas
idade_min <- min(caf$idade_lavoura_anos)
idade_max <- max(caf$idade_lavoura_anos)
idade_mediana <- median(caf$idade_lavoura_anos)
idade_media <- mean(caf$idade_lavoura_anos)
idade_sd <- sd(caf$idade_lavoura_anos)
idade_q1 <- quantile(caf$idade_lavoura_anos, 0.25)
idade_q3 <- quantile(caf$idade_lavoura_anos, 0.75)

# Histograma
idade_lavoura_anos_plot <- ggplot(caf, aes(x = idade_lavoura_anos)) +
  geom_histogram(binwidth = 2, colour = "white") +
  geom_vline(
  xintercept = 2,
  linetype = "dashed"
) +
annotate(
  "text",
  x = 2,
  y = Inf,
  label = "Jovens: até 2 anos",
  vjust = 1.5,
  hjust = -0.1
) +
  labs(
    x = "Idade da lavoura (anos)",
    y = "Frequência absoluta"
  ) +
  theme_minimal()

print(idade_lavoura_anos_plot)

# Resultados
# Exibir resultados
cat("\nIdade mínima:", idade_min)
cat("Idade máxima:", idade_max)
cat("1º quartil:", idade_q1)
cat("Mediana:", idade_mediana)
cat("Média:", round(idade_media, 2))
cat("3º quartil:", idade_q3)
cat("Desvio padrão:", round(idade_sd, 2))

# Parecer
cat("\n A idade da lavoura de café apresenta uma distribuição normal com assimetria
à direita. Existem lavouras jovens, com menos de dois anos e lavoura adultas com 31 anos.
Predominam lavouras adultas, com idade média de 9.25 anos.")

# VARIAVEL densidade_plantio
secao("Análise univariada - variável: densidade_plantio")

# Densidade de plantio sao plantas por alguma unidade de espaco
# Exemplo: plantas/metro quadrado, plantas/hectare
# Segundo Malta et al. (2008), 
# um hectare tradicional tem 2.500 plantas/ha, 
# semi-adensado possui 2.500 a 5.000 plantas/ha e 
# adensado possui mais de 5.000 plantas/ha

# Resumo estatistico
summary(caf$densidade_plantio)

# Densidade de plantio seguindo a classificacao de Malta et al. (2008)
caf2 <- caf %>%
  select(densidade_plantio) %>%
  mutate(
    densidade_plantio_classes = cut(
      densidade_plantio,
      breaks = c(-Inf, 2500, 5000, Inf),
      labels = c("Tradicional", "Semi-adensado", "Adensado"),
      right = TRUE
    ),
    densidade_plantio_classes = factor(
      densidade_plantio_classes,
      levels = c("Tradicional", "Semi-adensado", "Adensado")
    )
  )

# Grafico de barras agrupando talhoes pelas classes de densidade
densidade_plantio_plot <- ggplot(caf2, aes(x = densidade_plantio_classes)) +
  geom_bar(width = 0.7) +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    vjust = -0.3
  ) +
  labs(
    x = "Sistema de densidade por plantio",
    y = "Talhões"
  ) +
  theme_minimal(base_size = 13)

print(densidade_plantio_plot)

# Moda
densidade_moda <- names(which.max(table(caf2$densidade_plantio_classes)))

# Ranking
densidade_rank <- caf2 %>%
  group_by(densidade_plantio_classes) %>%
  summarise(frequencia = n()) %>%
  arrange(desc(frequencia))

# Exibir resultados
# Moda
cat("\nDos sistemas de plantio, predomina aqueles com densidade do tipo", 
str_to_title(densidade_moda))

# Ranking
cat("\nO ranking dos sistemas de densidadade de plantio é:")
print(densidade_rank)

# VARIAVEL adubacao_n_kg_h
secao("Análise univariada - variável: adubacao_n_kg_ha")

# Resumo estatistico
summary(caf$adubacao_n_kg_ha)

# Extrair estatisticas descritivas
adubacao_min <- min(caf$adubacao_n_kg_ha)
adubacao_max <- max(caf$adubacao_n_kg_ha)
adubacao_mediana <- median(caf$adubacao_n_kg_ha)
adubacao_media <- mean(caf$adubacao_n_kg_ha)
adubacao_sd <- sd(caf$adubacao_n_kg_ha)
adubacao_q1 <- quantile(caf$adubacao_n_kg_ha, 0.25)
adubacao_q3 <- quantile(caf$adubacao_n_kg_ha, 0.75)

# Histograma
adubacao_plot <- ggplot(caf, aes(x = adubacao_n_kg_ha)) +
  geom_histogram(binwidth = 20, colour = "white") +
  geom_vline(
  xintercept = adubacao_media,
  linetype = "dashed"
) +
annotate(
  "text",
  x = adubacao_media,
  y = Inf,
  label = paste0("Média = ", round(adubacao_media, 2), " kg/ha"),
  vjust = 1.5,
  hjust = -0.05
) +
  labs(
    x = "Adubação nitrogenada (kg/ha)",
    y = "Frequência absoluta"
  ) +
  theme_minimal()

print(adubacao_plot)

# Resultados
# Exibir resultados
cat("\nAdubação mínima:", adubacao_min)
cat("Adubação máxima:", adubacao_max)
cat("1º quartil:", adubacao_q1)
cat("Mediana:", adubacao_mediana)
cat("Média:", round(adubacao_media, 2))
cat("3º quartil:", adubacao_q3)
cat("Desvio padrão:", round(adubacao_sd, 2))

# Parecer
cat("\n.")

# VARIAVEL precipitacao_safra_mm
secao("Análise univariada - variável: precipitacao_safra_mm")

# Resumo estatistico
summary(caf$precipitacao_safra_mm)

# Extrair estatisticas descritivas
precipitacao_min <- min(caf$precipitacao_safra_mm)
precipitacao_max <- max(caf$precipitacao_safra_mm)
precipitacao_mediana <- median(caf$precipitacao_safra_mm)
precipitacao_media <- mean(caf$precipitacao_safra_mm)
precipitacao_sd <- sd(caf$precipitacao_safra_mm)
precipitacao_q1 <- quantile(caf$precipitacao_safra_mm, 0.25)
precipitacao_q3 <- quantile(caf$precipitacao_safra_mm, 0.75)

# Histograma
precipitacao_plot <- ggplot(caf, aes(x = precipitacao_safra_mm)) +
  geom_histogram(binwidth = 20, colour = "white") +
  geom_vline(
  xintercept = precipitacao_media,
  linetype = "dashed"
) +
annotate(
  "text",
  x = precipitacao_media,
  y = Inf,
  label = paste0("Média = ", round(precipitacao_media, 2), " mm"),
  vjust = 1.5,
  hjust = -0.05
) +
  labs(
    x = "Precipitação na safra (mm)",
    y = "Frequência absoluta"
  ) +
  theme_minimal()

print(precipitacao_plot)

# Resultados
# Exibir resultados
cat("\nPrecipitação mínima:", precipitacao_min)
cat("Precipitação máxima:", precipitacao_max)
cat("1º quartil:", precipitacao_q1)
cat("Mediana:", precipitacao_mediana)
cat("Média:", round(precipitacao_media, 2))
cat("3º quartil:", precipitacao_q3)
cat("Desvio padrão:", round(precipitacao_sd, 2))

# Parecer
cat("\n.")

# VARIAVEL umidade_relativa_pct
secao("Análise univariada - variável: umidade_relativa_pct")

# Resumo estatistico
summary(caf$umidade_relativa_pct)

# Extrair estatisticas descritivas
umidade_min <- min(caf$umidade_relativa_pct)
umidade_max <- max(caf$umidade_relativa_pct)
umidade_mediana <- median(caf$umidade_relativa_pct)
umidade_media <- mean(caf$umidade_relativa_pct)
umidade_sd <- sd(caf$umidade_relativa_pct)
umidade_q1 <- quantile(caf$umidade_relativa_pct, 0.25)
umidade_q3 <- quantile(caf$umidade_relativa_pct, 0.75)

# Histograma
umidade_plot <- ggplot(caf, aes(x = umidade_relativa_pct)) +
  geom_histogram(binwidth = 2) +
  geom_vline(
  xintercept = umidade_media,
  linetype = "dashed"
) +
annotate(
  "text",
  x = umidade_media,
  y = Inf,
  label = paste0("Média = ", round(umidade_media, 2), "%"),
  vjust = 1.5,
  hjust = -0.05
) +
  labs(
    x = "Umidade relativa (%)",
    y = "Frequência absoluta"
  ) +
  theme_minimal()

print(umidade_plot)

# Resultados
# Exibir resultados
cat("\nUmidade mínima:", umidade_min)
cat("Umidade máxima:", umidade_max)
cat("1º quartil:", umidade_q1)
cat("Mediana:", umidade_mediana)
cat("Média:", round(umidade_media, 2))
cat("3º quartil:", umidade_q3)
cat("Desvio padrão:", round(umidade_sd, 2))

# Parecer
cat("\n.")

# VARIAVEL ph_solo
secao("Análise univariada - variável: ph_solo")

# Resumo estatistico
summary(caf$ph_solo)

# Extrair estatisticas descritivas
ph_solo_min <- min(caf$ph_solo)
ph_solo_max <- max(caf$ph_solo)
ph_solo_mediana <- median(caf$ph_solo)
ph_solo_media <- mean(caf$ph_solo)
ph_solo_sd <- sd(caf$ph_solo)
ph_solo_q1 <- quantile(caf$ph_solo, 0.25)
ph_solo_q3 <- quantile(caf$ph_solo, 0.75)

# Histograma
ph_solo_plot <- ggplot(caf, aes(x = ph_solo)) +
  geom_histogram(binwidth = 0.25, colour = "white") +
  labs(
    x = "pH do solo",
    y = "Frequência absoluta"
  ) +
  theme_minimal()

print(ph_solo_plot)

# Resultados
# Exibir resultados
cat("\npH mínimo:", ph_solo_min)
cat("pH máximo:", ph_solo_max)
cat("1º quartil:", ph_solo_q1)
cat("Mediana:", ph_solo_mediana)
cat("Média:", round(ph_solo_media, 2))
cat("3º quartil:", ph_solo_q3)
cat("Desvio padrão:", round(ph_solo_sd, 2))

# Parecer
cat("\n.")

# VARIAVEL materia_organica_pct
secao("Análise univariada - variável: materia_organica_pct")

# Resumo estatistico
summary(caf$materia_organica_pct)

# Extrair estatisticas descritivas
materia_organica_min <- min(caf$materia_organica_pct)
materia_organica_max <- max(caf$materia_organica_pct)
materia_organica_mediana <- median(caf$materia_organica_pct)
materia_organica_media <- mean(caf$materia_organica_pct)
materia_organica_sd <- sd(caf$materia_organica_pct)
materia_organica_q1 <- quantile(caf$materia_organica_pct, 0.25)
materia_organica_q3 <- quantile(caf$materia_organica_pct, 0.75)

# Histograma
materia_organica_plot <- ggplot(caf, aes(x = materia_organica_pct)) +
  geom_histogram(binwidth = 0.5, colour = "white") +
  geom_vline(
  xintercept = materia_organica_media,
  linetype = "dashed"
) +
annotate(
  "text",
  x = materia_organica_media,
  y = Inf,
  label = paste0("Média = ", round(materia_organica_media, 2), "%"),
  vjust = 1.5,
  hjust = -0.05
) +
  labs(
    x = "Matéria orgânica (%)",
    y = "Frequência absoluta"
  ) +
  theme_minimal()

print(materia_organica_plot)

# Resultados
# Exibir resultados
cat("\nMatéria orgânica mínima:", materia_organica_min)
cat("Matéria orgânica máxima:", materia_organica_max)
cat("1º quartil:", materia_organica_q1)
cat("Mediana:", materia_organica_mediana)
cat("Média:", round(materia_organica_media, 2))
cat("3º quartil:", materia_organica_q3)
cat("Desvio padrão:", round(materia_organica_sd, 2))

# Parecer
cat("\n.")

# VARIAVEL 'n_armadilhas'
secao("Análise univariada - variável: n_armadilhas")
table(caf$n_armadilhas)

# Media
armadilhas_media <- mean(caf$n_armadilhas)

# Moda
armadilhas_moda <- names(which.max(table(caf$n_armadilhas)))

# Mediana
armadilhas_mediana <- median(caf$n_armadilhas)

# Ranking
armadilhas_rank <- caf %>%
  group_by(n_armadilhas) %>%
  summarise(frequencia = n()) %>%
  arrange(desc(frequencia))

# Grafico de barras
armadilhas_plot <- ggplot(caf, aes(x = n_armadilhas)) +
  geom_bar() +
  scale_x_continuous(breaks = seq(
    min(caf$n_armadilhas),
    max(caf$n_armadilhas),
    by = 1
  )) +
  theme_minimal() +
  labs(
    x = "Quantidade de armadilhas",
    y = "Frequência"
  )

print(armadilhas_plot)

# Exibir resultados
# Media, moda e mediana
cat("\nQuantidade média de armadilhas:", round(armadilhas_media, 2))
cat("\nQuantidade de armadilhas mais frequente (moda):", armadilhas_moda)
cat("\nMediana da quantidade de armadilhas:", armadilhas_mediana)

# Ranking
cat("\nO ranking da quantidade de armadilhas é:")
print(armadilhas_rank)

# Parecer
cat("\n.")

# VARIAVEL 'dias_exposicao'
secao("Análise univariada - variável: dias_exposicao")
table(caf$dias_exposicao)

# Media
exposicao_media <- mean(caf$dias_exposicao)

# Moda
exposicao_moda <- names(which.max(table(caf$dias_exposicao)))

# Mediana
exposicao_mediana <- median(caf$dias_exposicao)

# Ranking
dias_exposicao_rank <- caf %>%
  group_by(dias_exposicao) %>%
  summarise(frequencia = n()) %>%
  arrange(desc(frequencia))

# Grafico de barras
dias_exposicao_plot <- ggplot(
  dias_exposicao_rank,
  aes(x = dias_exposicao, y = frequencia)
) +
  geom_col() +
  labs(
    x = "Dias de exposição",
    y = "Frequência absoluta"
  ) +
  coord_flip() +
  theme_minimal()

print(dias_exposicao_plot)

# Exibir resultados
# Media, moda e mediana
cat("\nDias de exposição em média:", round(exposicao_media, 2))
cat("\nDias de exposição mais frequente (moda):", exposicao_moda)
cat("\nMediana de dias de exposição:", exposicao_mediana)

# Ranking
cat("\nO ranking dos dias de exposição é:")
print(dias_exposicao_rank)

# Parecer
cat("\n.")

# VARIAVEL 'n_brocas_capturadas'
secao("Análise univariada - variável: n_brocas_capturadas")
table(caf$n_brocas_capturadas)

# Moda
brocas_moda <- names(which.max(table(caf$n_brocas_capturadas)))

# Mediana
brocas_mediana <- median(caf$n_brocas_capturadas)

# Ranking
n_brocas_capturadas_rank <- caf %>%
  group_by(n_brocas_capturadas) %>%
  summarise(frequencia = n()) %>%
  arrange(desc(frequencia))

# Grafico de barras
n_brocas_capturadas_plot <- ggplot(caf, aes(x = n_brocas_capturadas)) +
  geom_bar() +
  scale_x_continuous(breaks = seq(
    min(caf$n_brocas_capturadas),
    max(caf$n_brocas_capturadas),
    by = 5
  )) +
  theme_minimal() +
  labs(
    x = "Número de brocas capturadas",
    y = "Frequência absoluta"
  )

print(n_brocas_capturadas_plot)

# Exibir resultados
# Moda e mediana
cat("\nQuantidade mais frequente de brocas capturadas (moda):", brocas_moda)
cat("\nMediana do número de brocas capturadas:", brocas_mediana)

# Ranking
cat("\nO ranking do número de brocas capturadas é:")
print(n_brocas_capturadas_rank)

# =====================================================================
# ANALISES MULTIVARIADAS
# =====================================================================

# TO DO

# =====================================================================
# ANALISES GLOBAIS
# =====================================================================

# TO DO

# Variaveis numericas
# num_caf <- caf %>%
#   select(where(is.numeric))
# num_caf <- names(num_caf)
 
# secao("Forma das distibuições - Cafeicultura")

# diagnostico_forma <- data.frame(
# variavel = num_caf,
# media = sapply(caf[num_caf], mean),
# mediana = sapply(caf[num_caf], median),
# assimetria = sapply(caf[num_caf], skewness),
# curtose = sapply(caf[num_caf], kurtosis),
# shapiro_p = sapply(caf[num_caf],
# function(x) shapiro.test(x)$p.value))
# print(diagnostico_forma, row.names = FALSE)

# # Matriz de dispersao

# p_pares <- ggpairs(
# caf[, num_caf],
# lower = list(
# continuous = wrap("points", alpha = 0.5, size = 1)
# ),
# upper = list(
# continuous = wrap("cor", method = "pearson", size = 4)
# ),
# diag = list(
# continuous = wrap("densityDiag", alpha = 0.5)
# ),
# title = "Matriz de dispersão - Pearson - Cafeicultura"
# )

# p_pares

# ggsave(
# "matriz_dispersao_pearson_cafeicultura.png",
# p_pares,
# width = 22,
# height = 22,
# units = "in",
# dpi = 300
# )
