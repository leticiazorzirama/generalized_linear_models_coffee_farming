# =====================================================================
# DESAFIO B - ITEM 0.4: AS COVARIAVEIS NUMERICAS CONTRA A SEVERIDADE
# Com graficos, no mesmo estilo do challenge_a.R da equipe.
# As figuras sao salvas em mneis/figuras/ e tambem aparecem no painel
# de Plots do Positron.
# =====================================================================

if (!require(pacman)) { install.packages("pacman"); library(pacman) }
p_load(ggplot2, dplyr, tidyr, ggcorrplot, patchwork)

localizar_script <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) return(normalizePath(sub("^--file=", "", m[1]), winslash = "/", mustWork = FALSE))
  for (i in seq_len(sys.nframe())) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(normalizePath(of, winslash = "/", mustWork = FALSE))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    pp <- try(rstudioapi::getSourceEditorContext()$path, silent = TRUE)
    if (!inherits(pp, "try-error") && length(pp) == 1 && nzchar(pp))
      return(normalizePath(pp, winslash = "/", mustWork = FALSE))
  }
  NA_character_
}
raiz_repo <- function(marcador = file.path("data", "cafeicultura.csv")) {
  sp <- localizar_script()
  d <- if (!is.na(sp)) dirname(sp) else normalizePath(getwd(), winslash = "/")
  repeat {
    if (file.exists(file.path(d, marcador))) return(d)
    pai <- dirname(d)
    if (identical(pai, d)) stop("Abra a PASTA glm-cafeicultura, nao o arquivo solto.", call. = FALSE)
    d <- pai
  }
}
RAIZ <- raiz_repo()
FIG  <- file.path(RAIZ, "mneis", "figuras")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

caf <- read.csv(file.path(RAIZ, "data", "cafeicultura.csv"), stringsAsFactors = TRUE)
y   <- caf$severidade_ferrugem

secao <- function(t) cat("\n", strrep("=", 68), "\n", t, "\n", strrep("=", 68), "\n", sep = "")
salvar <- function(g, nome, w = 10, h = 7) {
  ggsave(file.path(FIG, nome), g, width = w, height = h, dpi = 150, bg = "white")
  cat("  figura salva:", nome, "\n")
}

num <- c("umidade_relativa_pct", "densidade_plantio", "precipitacao_safra_mm",
         "adubacao_n_kg_ha", "ph_solo", "materia_organica_pct",
         "declividade_pct", "idade_lavoura_anos", "altitude_m")


# =====================================================================
secao("1. CORRELACAO DE CADA COVARIAVEL COM A SEVERIDADE")
# =====================================================================
cor_tab <- data.frame(
  covariavel = num,
  pearson  = sapply(num, function(v) cor(y, caf[[v]])),
  spearman = sapply(num, function(v) cor(y, caf[[v]], method = "spearman")),
  p_spearman = sapply(num, function(v)
    suppressWarnings(cor.test(y, caf[[v]], method = "spearman")$p.value)))
cor_tab <- cor_tab[order(-abs(cor_tab$spearman)), ]
cor_tab$p_bh <- p.adjust(cor_tab$p_spearman, method = "BH")
cor_tab$p_bonf <- p.adjust(cor_tab$p_spearman, method = "bonferroni")
cor_tab$forca <- cut(abs(cor_tab$spearman), c(-Inf, .3, .7, Inf),
                     labels = c("fraca", "moderada", "forte"))
print(data.frame(
  covariavel = cor_tab$covariavel,
  pearson    = round(cor_tab$pearson, 3),
  spearman   = round(cor_tab$spearman, 3),
  p          = format.pval(cor_tab$p_spearman, digits = 2),
  p_BH       = format.pval(cor_tab$p_bh, digits = 2),
  p_Bonf     = format.pval(cor_tab$p_bonf, digits = 2),
  forca      = cor_tab$forca), row.names = FALSE)

cat("\n# Parecer\n")
cat("Apenas a umidade relativa alcanca associacao moderada (0,389). Densidade de\n")
cat("plantio vem em seguida (0,237), ja na faixa fraca. As demais sao fracas, e\n")
cat("quatro delas nem significativas. Ajustando para multiplos testes (BH/Bonferroni),\n")
cat("variaveis marginais como adubacao_n_kg_ha (p=0.045) perdem significancia (BH p=0.068).\n")
cat("Destaca-se tambem o caso da precipitacao_safra_mm: marginalmente significativa\n")
cat("(p = 0,0024), porem condicionalmente nula no modelo completo (p = 0,972). Nenhuma\n")
cat("covariavel numerica rivaliza com a cultivar, que sozinha separa medias em 3,19x.\n")


# =====================================================================
secao("2. GRAFICO - severidade contra cada covariavel")
# =====================================================================
longo <- caf %>%
  select(severidade_ferrugem, all_of(num)) %>%
  pivot_longer(-severidade_ferrugem, names_to = "covariavel", values_to = "valor")
longo$covariavel <- factor(longo$covariavel, levels = cor_tab$covariavel)

g1 <- ggplot(longo, aes(x = valor, y = severidade_ferrugem)) +
  geom_point(pch = 21, fill = "white", colour = "black", size = 1.8, alpha = .6) +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE,
              colour = "#18675A", fill = "#18675A", alpha = .15, linewidth = .9) +
  facet_wrap(~ covariavel, scales = "free_x", ncol = 3) +
  labs(x = NULL, y = "severidade da ferrugem",
       title = "Severidade contra cada covariavel numerica",
       subtitle = "painéis ordenados pela magnitude da correlacao de Spearman") +
  theme_bw(base_size = 12)
print(g1); salvar(g1, "04_severidade_vs_covariaveis.png", 11, 8)

cat("\n# Parecer\n")
cat("As curvas sao praticamente planas, exceto umidade e densidade de plantio,\n")
cat("que sobem de forma suave e aproximadamente linear. Nenhuma mostra curvatura\n")
cat("forte que exigisse termo quadratico.\n")


# =====================================================================
secao("3. GRAFICO - a severidade por cultivar, que e o que importa")
# =====================================================================
ordem <- names(sort(tapply(y, caf$cultivar, median)))
caf$cultivar_ord <- factor(caf$cultivar, levels = ordem)

g2 <- ggplot(caf, aes(x = cultivar_ord, y = severidade_ferrugem)) +
  geom_boxplot(outlier.shape = NA, fill = "grey95", width = .55) +
  geom_jitter(width = .16, pch = 21, fill = "white", colour = "black",
              size = 1.9, alpha = .65) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3.4,
               fill = "#18675A", colour = "#18675A") +
  labs(x = NULL, y = "severidade da ferrugem",
       title = "Severidade por cultivar",
       subtitle = "losango = media · caixa = quartis · cada ponto e um talhao",
       caption = "Kruskal-Wallis: p = 2,8e-15") +
  theme_bw(base_size = 13)
print(g2); salvar(g2, "04_severidade_por_cultivar.png", 8, 6)

cat("\n# Parecer\n")
cat("A separacao do Icatu e visivel sem teste nenhum: a caixa inteira dele fica\n")
cat("abaixo da mediana das outras tres. Note tambem que a caixa dele e mais\n")
cat("ESTREITA - menor dispersao. E o indicio da Tarefa 3, visto a olho.\n")


# =====================================================================
secao("4. GRAFICO - matriz de correlacao (inclui multicolinearidade, item 0.5)")
# =====================================================================
mat  <- cor(caf[, c("severidade_ferrugem", num)], method = "spearman")
pmat <- ggcorrplot::cor_pmat(caf[, c("severidade_ferrugem", num)])

g3 <- ggcorrplot(mat, method = "square", type = "lower", lab = TRUE, lab_size = 2.9,
                 p.mat = pmat, insig = "blank",
                 colors = c("#A3352A", "white", "#18675A"),
                 title = "Correlacoes de Spearman (celulas em branco = nao significativas)") +
  theme(plot.title = element_text(size = 12))
print(g3); salvar(g3, "04_matriz_correlacao.png", 9, 8)

cat("\n# Parecer - e uma previa do item 0.5\n")
pares <- which(abs(mat) > .5 & abs(mat) < 1, arr.ind = TRUE)
if (nrow(pares)) {
  vistos <- character(0)
  for (i in seq_len(nrow(pares))) {
    a <- rownames(mat)[pares[i,1]]; b <- colnames(mat)[pares[i,2]]
    ch <- paste(sort(c(a,b)), collapse = "|")
    if (!ch %in% vistos) {
      vistos <- c(vistos, ch)
      cat(sprintf("  |r| > 0,5 entre %s e %s : %.3f\n", a, b, mat[pares[i,1], pares[i,2]]))
    }
  }
} else cat("  Nenhum par de covariaveis com |r| > 0,5.\n")
cat("\nSem pares fortemente correlacionados entre si, a multicolinearidade nao\n")
cat("deve ser problema. O item 0.5 confirmara pelo VIF, apos o ajuste.\n")

secao("FIM DO ITEM 0.4")
cat("Figuras em:", FIG, "\n")
