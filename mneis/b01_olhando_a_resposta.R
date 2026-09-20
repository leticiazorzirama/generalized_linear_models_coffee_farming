# =====================================================================
# DESAFIO B - PASSO 1: OLHANDO A VARIAVEL RESPOSTA
# severidade_ferrugem - so ela. Nenhum modelo ainda.
# =====================================================================

caf <- read.csv(file.path(
  if (file.exists("data/cafeicultura.csv")) "." else "..", "data/cafeicultura.csv"))
y <- caf$severidade_ferrugem

secao <- function(t) cat("\n", strrep("=", 62), "\n", t, "\n", strrep("=", 62), "\n", sep="")

secao("O QUE E")
cat("severidade_ferrugem = fracao da area foliar lesionada pelo fungo\n")
cat("medida por analise de imagem, talhao a talhao\n\n")
print(head(caf[, c("id_talhao", "cultivar", "severidade_ferrugem")], 5))

secao("ONDE ELA VIVE")
print(summary(y))
cat("\nminimo :", min(y), "\n")
cat("maximo :", max(y), "\n")
cat("exatamente 0:", sum(y == 0), " | exatamente 1:", sum(y == 1), "\n")

secao("AS TRES COISAS QUE ISSO DECIDE")
cat("1. E CONTINUA - nao e contagem.\n")
cat("   valores distintos:", length(unique(y)), "em", length(y), "talhoes\n")
cat("   -> Poisson esta fora. Ela so descreve 0, 1, 2, 3...\n\n")

cat("2. E LIMITADA nos dois lados - vive em (0, 1).\n")
cat("   -> qualquer modelo que possa predizer -0,05 ou 1,3 esta errado.\n\n")

cat("3. NAO TEM DENOMINADOR - ninguem contou folhas doentes de um total.\n")
cat("   -> Binomial esta fora: ela precisa do n de tentativas.\n")

secao("SOBRA A BETA")
cat("E a unica distribuicao das vistas em aula que e:\n")
cat("  contínua  +  definida exatamente em (0,1)  +  sem precisar de n\n")


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

secao("GRAFICO - onde a severidade vive")

g1 <- ggplot(caf, aes(x = severidade_ferrugem)) +
  geom_histogram(bins = 40, fill = "grey92", colour = "black", linewidth = .25) +
  geom_vline(xintercept = c(0, 1), colour = VERM, linetype = "dashed", linewidth = .7) +
  annotate("text", x = 0.015, y = Inf, label = "0", vjust = 1.6, hjust = 0,
           colour = VERM, size = 3.4) +
  annotate("text", x = 0.985, y = Inf, label = "1", vjust = 1.6, hjust = 1,
           colour = VERM, size = 3.4) +
  geom_vline(xintercept = mean(caf$severidade_ferrugem), colour = VERDE, linewidth = .8) +
  annotate("text", x = mean(caf$severidade_ferrugem) + .015, y = Inf,
           label = paste0("media ", round(mean(caf$severidade_ferrugem), 4)),
           vjust = 3.2, hjust = 0, colour = VERDE, size = 3.4) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
  labs(x = "severidade da ferrugem (fracao da area foliar lesionada)",
       y = "talhoes",
       title = "A resposta vive no intervalo aberto (0, 1)",
       subtitle = "assimetrica a direita, e longe de ocupar todo o intervalo disponivel",
       caption = "340 talhoes | faixa observada 0,0010 a 0,6098 | nenhum valor em 0 ou 1") +
  theme_bw(base_size = 12)
print(g1); salvar(g1, "01_distribuicao_da_resposta.png")

cat("\n# Parecer\n")
cat("A massa se concentra abaixo de 0,25 e a cauda direita e longa - assimetria\n")
cat("classica de proporcao baixa. E repare: os dados nem chegam perto de 1.\n")
cat("Mesmo assim o limite superior existe e o modelo tem de respeita-lo.\n")

secao("GRAFICO - por que nao e contagem")
cmp <- rbind(
  data.frame(v = caf$severidade_ferrugem, tipo = "severidade_ferrugem (continua)"),
  data.frame(v = caf$n_brocas_capturadas / max(caf$n_brocas_capturadas),
             tipo = "n_brocas_capturadas (contagem, reescalada)"))
g2 <- ggplot(cmp, aes(x = v)) +
  geom_dotplot(binwidth = .008, dotsize = .55, fill = "white", colour = "black",
               stackratio = .9) +
  facet_wrap(~ tipo, ncol = 1, scales = "free_y") +
  labs(x = "valor (reescalado para 0-1)", y = NULL,
       title = "Continua contra discreta",
       subtitle = "a contagem empilha em poucos valores; a severidade nao repete quase nunca") +
  theme_bw(base_size = 12) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
print(g2); salvar(g2, "01_continua_vs_discreta.png", 9, 6)

cat("\n# Parecer\n")
cat("A diferenca fica obvia: a contagem forma colunas sobre valores inteiros,\n")
cat("a severidade espalha sem repetir. 320 valores distintos em 340 talhoes.\n")
cat("Nenhuma distribuicao discreta descreve isso.\n")
