# =====================================================================
# DESAFIO B - PASSO 2: A SEVERIDADE MUDA CONFORME OS FATORES?
# Ainda sem modelo. So descricao e um teste nao parametrico.
# =====================================================================

caf <- read.csv(file.path(
  if (file.exists("data/cafeicultura.csv")) "." else "..", "data/cafeicultura.csv"))
y <- caf$severidade_ferrugem

secao <- function(t) cat("\n", strrep("=", 66), "\n", t, "\n", strrep("=", 66), "\n", sep="")

# Funcao auxiliar: descreve y dentro de cada nivel de um fator
por_fator <- function(f, nome) {
  s <- split(y, caf[[f]])
  tab <- do.call(rbind, lapply(s, function(v) data.frame(
    n = length(v), media = round(mean(v), 4), mediana = round(median(v), 4),
    dp = round(sd(v), 4), min = round(min(v), 4), max = round(max(v), 4))))
  tab <- tab[order(tab$media), ]
  kw <- kruskal.test(y ~ factor(caf[[f]]))
  secao(paste("SEVERIDADE POR", toupper(nome)))
  print(tab)
  cat("\nKruskal-Wallis: chi2 =", round(kw$statistic, 2),
      ", gl =", kw$parameter,
      ", p =", format.pval(kw$p.value, digits = 3), "\n")
  cat(if (kw$p.value < 0.05)
        "-> ha diferenca entre os niveis.\n"
      else "-> sem evidencia de diferenca entre os niveis.\n")
  invisible(tab)
}

t_cult <- por_fator("cultivar",         "cultivar")
t_man  <- por_fator("manejo",           "manejo")
t_irr  <- por_fator("irrigacao",        "irrigacao")
t_reg  <- por_fator("regiao_produtora", "regiao produtora")

secao("COMPARANDO A FORCA DOS QUATRO FATORES")
amp <- function(t) round(max(t$media) / min(t$media), 2)
cmp <- data.frame(
  fator = c("cultivar", "manejo", "irrigacao", "regiao_produtora"),
  niveis = c(nrow(t_cult), nrow(t_man), nrow(t_irr), nrow(t_reg)),
  menor_media = c(min(t_cult$media), min(t_man$media), min(t_irr$media), min(t_reg$media)),
  maior_media = c(max(t_cult$media), max(t_man$media), max(t_irr$media), max(t_reg$media)),
  razao = c(amp(t_cult), amp(t_man), amp(t_irr), amp(t_reg)),
  p_kw = c(kruskal.test(y ~ factor(caf$cultivar))$p.value,
           kruskal.test(y ~ factor(caf$manejo))$p.value,
           kruskal.test(y ~ factor(caf$irrigacao))$p.value,
           kruskal.test(y ~ factor(caf$regiao_produtora))$p.value))
cmp$p_holm <- format.pval(p.adjust(cmp$p_kw, method = "holm"), digits = 2)
cmp$p_bonf <- format.pval(p.adjust(cmp$p_kw, method = "bonferroni"), digits = 2)
cmp$p_kw <- format.pval(cmp$p_kw, digits = 2)
cmp <- cmp[order(-cmp$razao), ]
print(cmp, row.names = FALSE)

secao("PARECER DO PASSO 2")
cat("A cultivar e, de longe, o fator que mais move a severidade: o pior nivel\n")
cat("tem media", amp(t_cult), "vezes a do melhor. Nenhum outro fator chega perto.\n\n")
cat("O Icatu tem a menor severidade media (", min(t_cult$media), ") e tambem o\n", sep="")
cat("menor desvio-padrao (", t_cult$dp[1], "). Ele nao so adoece menos: adoece\n", sep="")
cat("de forma mais previsivel. Guardar isso - vai importar na Tarefa 3.\n\n")
cat("ATENCAO: isto e descricao, nao conclusao. As cultivares podem estar\n")
cat("distribuidas de forma desigual entre regioes, altitudes ou manejos, e\n")
cat("parte da diferenca pode vir dai. So o modelo separa o efeito da cultivar\n")
cat("do efeito das companhias dela.\n")


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

secao("GRAFICO - severidade por cada fator categorico")

fatores <- c("cultivar", "regiao_produtora", "irrigacao", "manejo")
paineis <- lapply(fatores, function(f) {
  ordem <- names(sort(tapply(y, caf[[f]], median)))
  d <- caf; d$g <- factor(d[[f]], levels = ordem)
  kw <- kruskal.test(y ~ factor(caf[[f]]))
  ggplot(d, aes(x = g, y = severidade_ferrugem)) +
    geom_boxplot(outlier.shape = NA, fill = "grey95", width = .55) +
    geom_jitter(width = .15, pch = 21, fill = "white", colour = "black",
                size = 1.4, alpha = .5) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 2.6,
                 fill = VERDE, colour = VERDE) +
    labs(x = NULL, y = NULL, title = f,
         subtitle = paste0("Kruskal-Wallis p = ", format.pval(kw$p.value, digits = 2))) +
    theme_bw(base_size = 11) +
    theme(plot.subtitle = element_text(
      colour = if (kw$p.value < .05) VERDE else VERM))
})
g <- wrap_plots(paineis, ncol = 2) +
  plot_annotation(
    title = "Severidade por fator categorico",
    subtitle = "niveis ordenados pela mediana | losango = media | cada ponto e um talhao",
    theme = theme(plot.title = element_text(size = 15, face = "bold")))
print(g); salvar(g, "02_severidade_por_fator.png", 11, 8)

cat("\n# Parecer\n")
cat("O contraste entre os quatro paineis e o achado: no de cultivar as caixas\n")
cat("estao claramente escalonadas; nos outros tres elas praticamente se sobrepoem.\n")
cat("No manejo, os tres niveis sao quase indistinguiveis - e o p = 0,889 confirma.\n")

secao("GRAFICO - a forca de cada fator, lado a lado")
amp <- data.frame(
  fator = fatores,
  razao = sapply(fatores, function(f) {
    m <- tapply(y, caf[[f]], mean); max(m)/min(m) }),
  p = sapply(fatores, function(f) kruskal.test(y ~ factor(caf[[f]]))$p.value))
amp$fator <- factor(amp$fator, levels = amp$fator[order(amp$razao)])
amp$rot <- sprintf("%.2fx   p = %s", amp$razao, format.pval(amp$p, digits = 2))

g2 <- ggplot(amp, aes(x = razao, y = fator)) +
  geom_vline(xintercept = 1, colour = "grey60", linetype = "dashed") +
  geom_segment(aes(x = 1, xend = razao, yend = fator), colour = VERDE, linewidth = 1.4) +
  geom_point(size = 4, colour = VERDE) +
  geom_text(aes(label = rot), hjust = -0.15, size = 3.5) +
  scale_x_continuous(limits = c(0.95, 4.2), breaks = 1:4) +
  labs(x = "razao entre a maior e a menor media do fator", y = NULL,
       title = "Quanto cada fator move a severidade",
       subtitle = "1,00x = nenhum efeito | a cultivar esta sozinha no topo") +
  theme_bw(base_size = 12)
print(g2); salvar(g2, "02_forca_dos_fatores.png", 9, 5)

cat("\n# Parecer\n")
cat("A cultivar sozinha move 3,19x. Os outros tres fatores ficam todos abaixo\n")
cat("de 1,5x - e o manejo, em 1,11x, esta praticamente colado no 1,00.\n")
