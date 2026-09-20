# =====================================================================
# O MODELO NULO
# O modelo que nao explica nada - e por que ele e indispensavel.
# Rode e confira. Tudo aqui fecha na mao.
# =====================================================================

localizar_script <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) return(normalizePath(sub("^--file=", "", m[1]), winslash = "/", mustWork = FALSE))
  for (i in seq_len(sys.nframe())) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(normalizePath(of, winslash = "/", mustWork = FALSE))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- try(rstudioapi::getSourceEditorContext()$path, silent = TRUE)
    if (!inherits(p, "try-error") && length(p) == 1 && nzchar(p))
      return(normalizePath(p, winslash = "/", mustWork = FALSE))
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
caf  <- read.csv(file.path(RAIZ, "data", "cafeicultura.csv"))
caf$esforco <- caf$n_armadilhas * caf$dias_exposicao
y <- caf$n_brocas_capturadas

sec <- function(t) cat("\n", strrep("=", 68), "\n", t, "\n", strrep("=", 68), "\n", sep = "")


# =====================================================================
sec("1. O QUE E: o modelo sem nenhuma covariavel")
# =====================================================================
# A formula e  y ~ 1.  O "1" e o intercepto, e nada mais entra.
#
#     log(mu_i) = beta0          para TODO talhao i
#
# Repare que nao ha indice i do lado direito: o modelo afirma que todos
# os 340 talhoes tem a MESMA media esperada. Nada explica nada.

m_nulo <- glm(y ~ 1, family = poisson(link = "log"))

cat("coeficiente estimado (beta0) :", round(coef(m_nulo)[[1]], 6), "\n")
cat("exp(beta0)                   :", round(exp(coef(m_nulo)[[1]]), 6), "\n")
cat("media dos dados  mean(y)     :", round(mean(y), 6), "  <- o MESMO numero\n")
cat("\nvalores ajustados - todos identicos?", length(unique(round(fitted(m_nulo), 10))) == 1, "\n")
cat("  os 5 primeiros:", round(head(fitted(m_nulo), 5), 4), "\n")

cat("\n>> O modelo nulo com link log e apenas a MEDIA GERAL, escrita em\n")
cat("   linguagem de MLG. exp(beta0) = mean(y), sempre.\n")


# =====================================================================
sec("2. COM OFFSET, ele vira a TAXA GERAL")
# =====================================================================
m_nulo_off <- glm(y ~ 1, family = poisson, offset = log(caf$esforco))

cat("exp(beta0) do modelo nulo com offset :", round(exp(coef(m_nulo_off)[[1]]), 8), "\n")
cat("sum(y) / sum(esforco)                :", round(sum(y) / sum(caf$esforco), 8), "\n")
cat("  <- a taxa media global, brocas por armadilha-dia\n")

cat("\nOu seja:\n")
cat("  y ~ 1                 -> a media geral      (", round(mean(y), 3), "brocas )\n")
cat("  y ~ 1 + offset(log E) -> a taxa media geral (", round(sum(y)/sum(caf$esforco), 5),
    "por armadilha-dia )\n")


# =====================================================================
sec("3. PARA QUE SERVE: e a REGUA de comparacao")
# =====================================================================
# Sozinho ele nao interessa. Ele existe para responder:
#   "meu modelo explica alguma coisa, ou eu poderia ter chutado a media?"

f <- y ~ regiao_produtora + cultivar + manejo + irrigacao + altitude_m +
  declividade_pct + idade_lavoura_anos + densidade_plantio + adubacao_n_kg_ha +
  precipitacao_safra_mm + umidade_relativa_pct + ph_solo + materia_organica_pct

m_completo <- glm(f, family = poisson, data = caf, offset = log(esforco))

cat("Deviance = o quanto o modelo ERRA. Menor e melhor.\n\n")
cat("  deviance do modelo NULO     :", round(m_completo$null.deviance, 2), "\n")
cat("  deviance do modelo COMPLETO :", round(m_completo$deviance, 2), "\n")
cat("  reducao                     :",
    round(100 * (1 - m_completo$deviance / m_completo$null.deviance), 1), "%\n")

cat("\n>> Esses dois numeros ja vem de graca em QUALQUER glm - veja summary():\n")
cat("   'Null deviance' e 'Residual deviance'. O modelo nulo ja foi ajustado\n")
cat("   por baixo dos panos, mesmo que voce nunca o tenha escrito.\n")

cat("\nE o teste formal (a analise de deviance da aula de 08/09):\n")
m_nulo_mesmo_dado <- glm(y ~ 1, family = poisson, data = caf, offset = log(esforco))
print(anova(m_nulo_mesmo_dado, m_completo, test = "Chisq"))


# =====================================================================
sec("4. ONDE ELE JA APARECEU SEM VOCE PERCEBER")
# =====================================================================
razao <- function(m) sum(residuals(m, type = "pearson")^2) / df.residual(m)
cat("Na tabela de superdispersao, as duas primeiras linhas eram modelos nulos:\n\n")
cat("  y ~ 1                 -> dispersao", round(razao(m_nulo), 3), "\n")
cat("  y ~ 1 + offset(log E) -> dispersao", round(razao(m_nulo_off), 3), "\n")
cat("  completo              -> dispersao", round(razao(m_completo), 3), "\n")

cat("\nE no dredge() da Tarefa 4, o modelo nulo e uma das combinacoes testadas:\n")
cat("  ele e o caso 'nenhuma variavel'. Se ele vencer no AIC, nenhuma\n")
cat("  covariavel valeu o parametro que gastou.\n")

cat("\nAIC do nulo    :", round(AIC(m_nulo_off), 1), "\n")
cat("AIC do completo:", round(AIC(m_completo), 1), "\n")

sec("FIM")
