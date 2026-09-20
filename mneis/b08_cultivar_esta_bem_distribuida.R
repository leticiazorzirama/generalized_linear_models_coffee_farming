# =====================================================================
# DESAFIO B - PASSO 3: A CULTIVAR ESTA BEM DISTRIBUIDA?
# Se o Icatu estiver concentrado em lugares favoraveis, parte da vantagem
# dele vem de ONDE foi plantado, nao do que ele e.
# =====================================================================
caf <- read.csv(file.path(
  if (file.exists("data/cafeicultura.csv")) "." else "..", "data/cafeicultura.csv"))
secao <- function(t) cat("\n", strrep("=", 68), "\n", t, "\n", strrep("=", 68), "\n", sep="")

secao("1. A CULTIVAR CRUZADA COM OS OUTROS FATORES")
for (f in c("regiao_produtora", "manejo", "irrigacao")) {
  cat("\ncultivar x", f, "- em % de cada linha:\n")
  tb <- table(caf$cultivar, caf[[f]])
  print(round(100 * prop.table(tb, 1), 1))
  cs <- suppressWarnings(chisq.test(tb))
  cat("qui-quadrado de Pearson: X2 =", round(cs$statistic, 2),
      ", gl =", cs$parameter, ", p =", format.pval(cs$p.value, digits = 3), "\n")
  cat(if (cs$p.value < 0.05) "  -> DESBALANCEADO\n" else "  -> equilibrado\n")
}

secao("2. A CULTIVAR CRUZADA COM AS NUMERICAS")
cat("Media de cada covariavel dentro de cada cultivar:\n\n")
num <- c("altitude_m","umidade_relativa_pct","precipitacao_safra_mm",
         "densidade_plantio","idade_lavoura_anos","adubacao_n_kg_ha",
         "ph_solo","materia_organica_pct","declividade_pct")
tab <- t(sapply(num, function(v) tapply(caf[[v]], caf$cultivar, mean)))
print(round(tab, 2))

cat("\nHa diferenca entre as cultivares? (Kruskal-Wallis por covariavel)\n\n")
res <- data.frame(
  covariavel = num,
  p = sapply(num, function(v) kruskal.test(caf[[v]] ~ factor(caf$cultivar))$p.value),
  amplitude_pct = sapply(num, function(v) {
    m <- tapply(caf[[v]], caf$cultivar, mean); round(100*(max(m)/min(m) - 1), 1)}))
res <- res[order(res$p), ]
res$situacao <- ifelse(res$p < 0.05, "<<< DESBALANCEADO", "equilibrado")
res$p <- format.pval(res$p, digits = 2)
print(res, row.names = FALSE)

secao("3. O TESTE QUE IMPORTA: o Icatu esta em lugar melhor?")
cat("A umidade foi a covariavel mais associada a severidade (Spearman 0,389).\n")
cat("Se o Icatu estiver em talhoes mais secos, parte da vantagem e do lugar.\n\n")
u <- tapply(caf$umidade_relativa_pct, caf$cultivar, mean)
s <- tapply(caf$severidade_ferrugem,  caf$cultivar, mean)
print(data.frame(umidade_media = round(u,2), severidade_media = round(s,4)))
cat("\ncorrelacao entre as duas colunas acima (4 pontos):",
    round(cor(u, s), 3), "\n")
cat("Kruskal-Wallis umidade ~ cultivar: p =",
    format.pval(kruskal.test(caf$umidade_relativa_pct ~ factor(caf$cultivar))$p.value,
                digits = 3), "\n")

secao("PARECER DO PASSO 3")
cat("O desenho e proximo de equilibrado: as cultivares se distribuem de forma\n")
cat("parecida entre regioes, manejos e irrigacao, e as covariaveis numericas\n")
cat("nao diferem sistematicamente entre elas.\n\n")
cat("Isso e boa noticia para a interpretacao: a vantagem do Icatu NAO parece\n")
cat("ser artefato de ele ter sido plantado em lugares melhores.\n\n")
cat("MAS isso nao dispensa o modelo. Equilibrio marginal nao garante ausencia\n")
cat("de confundimento conjunto, e so o ajuste multivariado estima o efeito da\n")
cat("cultivar mantendo as demais constantes. O que este passo faz e reduzir a\n")
cat("suspeita, nao eliminar a necessidade.\n")


# =====================================================================
# GRAFICOS
# =====================================================================
if (!require(pacman)) { install.packages("pacman"); library(pacman) }
p_load(ggplot2, dplyr, tidyr, patchwork)

RAIZ <- if (file.exists("data/cafeicultura.csv")) "." else ".."
FIG  <- file.path(RAIZ, "mneis", "figuras")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
salvar <- function(g, nome, w = 9, h = 6) {
  ggsave(file.path(FIG, nome), g, width = w, height = h, dpi = 150, bg = "white")
  cat("  figura salva:", nome, "\n")
}
VERDE <- "#18675A"; VERM <- "#A3352A"

secao("GRAFICO - a cultivar esta bem distribuida?")

paineis <- lapply(c("regiao_produtora", "manejo", "irrigacao"), function(f) {
  d <- as.data.frame(prop.table(table(caf$cultivar, caf[[f]]), 1))
  names(d) <- c("cultivar", "nivel", "prop")
  cs <- suppressWarnings(chisq.test(table(caf$cultivar, caf[[f]])))
  ggplot(d, aes(x = cultivar, y = prop, fill = nivel)) +
    geom_col(width = .7, colour = "white", linewidth = .3) +
    scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
    scale_fill_brewer(palette = "Greens", direction = -1) +
    labs(x = NULL, y = NULL, fill = NULL, title = f,
         subtitle = paste0("qui-quadrado p = ", format.pval(cs$p.value, digits = 2),
                           if (cs$p.value >= .05) "  (equilibrado)" else "  (DESBALANCEADO)")) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", legend.text = element_text(size = 8),
          axis.text.x = element_text(angle = 20, hjust = 1))
})
g <- wrap_plots(paineis, ncol = 3) +
  plot_annotation(
    title = "Composicao de cada cultivar",
    subtitle = "se as barras tem a mesma reparticao, o desenho esta equilibrado",
    theme = theme(plot.title = element_text(size = 15, face = "bold")))
print(g); salvar(g, "03_composicao_das_cultivares.png", 12, 6)

cat("\n# Parecer\n")
cat("As quatro barras de cada painel tem reparticao parecida. Nao ha cultivar\n")
cat("concentrada numa regiao, num manejo ou na irrigacao.\n")

secao("GRAFICO - o teste que mais importa: umidade contra severidade")
u <- tapply(caf$umidade_relativa_pct, caf$cultivar, mean)
s <- tapply(caf$severidade_ferrugem,  caf$cultivar, mean)
d <- data.frame(cultivar = names(u), umidade = as.vector(u), severidade = as.vector(s))

g2 <- ggplot(d, aes(x = umidade, y = severidade)) +
  geom_point(size = 5, pch = 21, fill = VERDE, colour = "black") +
  geom_text(aes(label = cultivar), vjust = -1.3, size = 4) +
  scale_x_continuous(limits = c(73.0, 74.8)) +
  scale_y_continuous(limits = c(0.04, 0.24)) +
  labs(x = "umidade relativa media do talhao (%)", y = "severidade media",
       title = "O Icatu esta num ambiente melhor?",
       subtitle = "a umidade varia 1,4% entre cultivares; a severidade varia 3,2x",
       caption = "ATENCAO: sao 4 pontos. A correlacao de 0,893 entre eles nao e evidencia de nada.") +
  theme_bw(base_size = 12)
print(g2); salvar(g2, "03_umidade_vs_severidade.png", 8, 6)

cat("\n# Parecer\n")
cat("O eixo x cobre menos de 2 pontos percentuais de umidade; o eixo y cobre\n")
cat("uma variacao de 3,2 vezes. As escalas denunciam: a diferenca de ambiente\n")
cat("e pequena demais para explicar a diferenca de doenca.\n")
