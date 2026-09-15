# =====================================================================
# DESAFIO A - PONTO DE PARTIDA
# Objetivo: ver com os proprios olhos os 4 fatos que governam o Desafio A.
# Nao produz o modelo final. Produz ENTENDIMENTO.
# Rode inteiro (Ctrl+Shift+Enter no Positron) e leia a saida de cima a baixo.
# =====================================================================

library(MASS)   # glm.nb
library(AER)    # dispersiontest

sec <- function(txt) cat("\n", strrep("=", 72), "\n", txt, "\n",
                         strrep("=", 72), "\n", sep = "")

# --- Localizacao dos dados -------------------------------------------------
# O script acha o CSV sozinho, nao importa como voce executa: pasta aberta como
# projeto, arquivo solto com Source, Ctrl+Enter linha a linha, ou Rscript.

localizar_script <- function() {
  # a) Rscript / R CMD BATCH
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) return(normalizePath(sub("^--file=", "", m[1]), winslash = "/", mustWork = FALSE))
  # b) source() - e o que o botao Source do Positron/RStudio usa
  for (i in seq_len(sys.nframe())) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(normalizePath(of, winslash = "/", mustWork = FALSE))
  }
  # c) editor ativo da IDE - cobre o Ctrl+Enter linha a linha
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- try(rstudioapi::getSourceEditorContext()$path, silent = TRUE)
    if (!inherits(p, "try-error") && length(p) == 1 && nzchar(p))
      return(normalizePath(p, winslash = "/", mustWork = FALSE))
  }
  NA_character_   # colado direto no console: nao ha arquivo
}

raiz_repo <- function(marcador = file.path("data", "cafeicultura.csv")) {
  # Parte do diretorio do proprio script; se nao der para saber, do cwd.
  spath <- localizar_script()
  d <- if (!is.na(spath)) dirname(spath) else normalizePath(getwd(), winslash = "/")
  repeat {
    if (file.exists(file.path(d, marcador))) return(d)
    pai <- dirname(d)
    if (identical(pai, d))
      stop("Nao encontrei '", marcador, "'.
",
           "  Abra a PASTA glm-cafeicultura no Positron (nao o arquivo solto),
",
           "  ou use o botao Source em vez de colar o codigo no console.",
           call. = FALSE)
    d <- pai
  }
}

RAIZ <- raiz_repo()
caf  <- read.csv(file.path(RAIZ, "data", "cafeicultura.csv"))
caf$esforco <- caf$n_armadilhas * caf$dias_exposicao

# Formula com TODAS as candidatas numericas + as categoricas do delineamento.
f <- n_brocas_capturadas ~ regiao_produtora + cultivar + manejo + irrigacao +
  altitude_m + declividade_pct + idade_lavoura_anos + densidade_plantio +
  adubacao_n_kg_ha + precipitacao_safra_mm + umidade_relativa_pct +
  ph_solo + materia_organica_pct


# ---------------------------------------------------------------------
sec("FATO 1 - o esforco amostral nao e constante entre talhoes")
# ---------------------------------------------------------------------
cat("armadilhas variam de", min(caf$n_armadilhas), "a", max(caf$n_armadilhas), "\n")
cat("dias variam de      ", min(caf$dias_exposicao), "a", max(caf$dias_exposicao), "\n")
cat("esforco = armadilhas x dias, de", min(caf$esforco), "a", max(caf$esforco),
    "->", max(caf$esforco) / min(caf$esforco), "vezes\n\n")
cat("Dois talhoes reais para comparar:\n")
print(caf[c(which.min(caf$esforco)[1], which.max(caf$esforco)[1]),
          c("id_talhao", "n_armadilhas", "dias_exposicao", "esforco",
            "n_brocas_capturadas")], row.names = FALSE)
cat("\n-> comparar as contagens brutas desses dois e comparar coisas diferentes.\n")


# ---------------------------------------------------------------------
sec("FATO 2 - ignorar o esforco contamina ate a analise exploratoria")
# ---------------------------------------------------------------------
caf$taxa <- caf$n_brocas_capturadas / caf$esforco

# Spearman, e nao Pearson, por um motivo especifico: ele so olha a ORDEM dos
# valores (postos), entao nao muda sob qualquer transformacao monotonica.
# Como vamos modelar na escala log, isso importa - veja a comparacao abaixo.

cat("Ranking de associacao com a CONTAGEM BRUTA (Spearman):

")
vs <- c("esforco", "dias_exposicao", "n_armadilhas", "umidade_relativa_pct",
        "densidade_plantio", "altitude_m", "materia_organica_pct",
        "precipitacao_safra_mm", "adubacao_n_kg_ha", "declividade_pct",
        "idade_lavoura_anos", "ph_solo")
r <- sapply(vs, function(v) cor(caf$n_brocas_capturadas, caf[[v]], method = "spearman"))
o <- order(-abs(r))
print(data.frame(
  variavel = vs[o],
  spearman = round(r[o], 3),
  natureza = ifelse(vs[o] %in% c("esforco", "n_armadilhas", "dias_exposicao"),
                    "<<< PROTOCOLO", "    agronomica")), row.names = FALSE)

cat("
-> AQUI esta a conclusao: as tres primeiras posicoes em magnitude sao
")
cat("   esforco (0.488), dias_exposicao (0.350) e n_armadilhas (0.343). Nenhuma
")
cat("   delas e propriedade do talhao - sao decisoes de quem foi a campo montar
")
cat("   o monitoramento. A melhor variavel agronomica (umidade, 0.343) apenas
")
cat("   empata com n_armadilhas.
")
cat("   Quem ordenar variaveis por associacao com a contagem crua vai priorizar
")
cat("   o protocolo de coleta achando que esta priorizando agronomia.

")

cat("Por que Spearman e nao Pearson - o teste:
")
cmp <- data.frame(
  escala   = c("bruto (esforco, contagem)", "log   (log esforco, log contagem+1)"),
  pearson  = round(c(cor(caf$esforco, caf$n_brocas_capturadas),
                     cor(log(caf$esforco), log(caf$n_brocas_capturadas + 1))), 3),
  spearman = round(c(cor(caf$esforco, caf$n_brocas_capturadas, method = "spearman"),
                     cor(log(caf$esforco), log(caf$n_brocas_capturadas + 1),
                         method = "spearman")), 3))
print(cmp, row.names = FALSE)
cat("-> Pearson muda conforme a escala; Spearman e identico. Como o modelo vive
")
cat("   na escala log, Spearman mede na mesma lingua do modelo.

")

vs2 <- c("densidade_plantio", "umidade_relativa_pct", "ph_solo", "altitude_m")
cmp2 <- data.frame(
  variavel   = vs2,
  r_contagem = round(sapply(vs2, function(v) cor(caf$n_brocas_capturadas, caf[[v]], method = "spearman")), 3),
  r_taxa     = round(sapply(vs2, function(v) cor(caf$taxa, caf[[v]], method = "spearman")), 3))
print(cmp2, row.names = FALSE)
cat("
-> olhando a TAXA, o protocolo sai de cena e a hierarquia agronomica muda.

")


# ---------------------------------------------------------------------
sec("FATO 3 - offset sem log falha EM SILENCIO")
# ---------------------------------------------------------------------
mod_errado <- suppressWarnings(
  glm(f, family = poisson(link = "log"), data = caf, offset = esforco))

cat("mod_errado <- glm(..., offset = esforco)   # sem log\n\n")
cat("  converged :", mod_errado$converged, "  <- FALSO\n")
cat("  iteracoes :", mod_errado$iter, "(bateu no teto de", mod_errado$control$maxit, ")\n")
cat("  AIC       :", format(AIC(mod_errado), digits = 3), "\n")
cat("  class     :", paste(class(mod_errado), collapse = ", "), "\n\n")
cat("summary() roda normalmente e imprime coeficientes e p-valores.\n")
cat("Nada grita. Por isso a linha abaixo entra no topo do script final:\n")
cat("    options(warn = 2)   # transforma warning em erro\n")

# Por que precisa do log: o preditor linear esta na escala log.
#   log(mu_i) = log(E_i) + x_i'beta
# Sem o log, o R soma 60, 150, 360... direto em eta -> exp(360) estoura.


# ---------------------------------------------------------------------
sec("FATO 4 - com o offset correto, ha superdispersao severa")
# ---------------------------------------------------------------------
mod_pois <- glm(f, family = poisson(link = "log"), data = caf,
                offset = log(esforco))
stopifnot(mod_pois$converged)   # agora sim

razao <- sum(residuals(mod_pois, type = "pearson")^2) / df.residual(mod_pois)
cat("Tarefa 2 do enunciado:\n")
cat("  razao de Pearson / gl =", round(razao, 3), " (equidispersao seria ~1)\n\n")
print(dispersiontest(mod_pois))

mod_nb <- glm.nb(update(f, . ~ . + offset(log(esforco))), data = caf)
cat("\nTarefa 3 - binomial negativa:\n")
cat("  theta estimado :", round(mod_nb$theta, 3), "\n")
cat("  AIC Poisson    :", round(AIC(mod_pois), 1), "\n")
cat("  AIC Bin.Neg.   :", round(AIC(mod_nb), 1), "\n")


# ---------------------------------------------------------------------
sec("FATO 5 - a resposta da Tarefa 1 depende da especificacao do modelo")
# ---------------------------------------------------------------------
# O enunciado pede: soltar o log(esforco) como preditor livre e comparar com 1.
# Fazemos isso sob DUAS especificacoes da media e DUAS familias, porque o
# resultado nao e o mesmo nas quatro combinacoes - e esse e o achado.

f_num <- n_brocas_capturadas ~ altitude_m + declividade_pct +
  idade_lavoura_anos + densidade_plantio + adubacao_n_kg_ha +
  precipitacao_safra_mm + umidade_relativa_pct + ph_solo + materia_organica_pct

teste_gamma <- function(modelo, nome) {
  g  <- coef(modelo)[["log(esforco)"]]
  ic <- confint.default(modelo)["log(esforco)", ]
  cat(sprintf("    %-16s gamma = %.3f  IC95%% [%.3f ; %.3f]  -> %s
",
              nome, g, ic[1], ic[2],
              ifelse(ic[1] <= 1 & ic[2] >= 1, "contem 1", "REJEITA 1")))
}

disp <- function(f) {
  m <- glm(f, family = poisson, data = caf, offset = log(esforco))
  sum(residuals(m, type = "pearson")^2) / df.residual(m)
}

cat("(A) so covariaveis numericas  [ = o que challenge_a.R faz hoje ]
")
cat(sprintf("    dispersao = %.2f
", disp(f_num)))
teste_gamma(glm(update(f_num, . ~ . + log(esforco)), family = poisson, data = caf), "Poisson")
teste_gamma(glm.nb(update(f_num, . ~ . + log(esforco)), data = caf), "Bin.Negativa")

cat("
(B) com as categoricas do delineamento
")
cat(sprintf("    dispersao = %.2f
", disp(f)))
teste_gamma(glm(update(f, . ~ . + log(esforco)), family = poisson, data = caf), "Poisson")
teste_gamma(glm.nb(update(f, . ~ . + log(esforco)), data = caf), "Bin.Negativa")

cat("
Quanto as categoricas explicam:
")
print(anova(glm(f_num, family = poisson, data = caf, offset = log(esforco)),
            glm(f,     family = poisson, data = caf, offset = log(esforco)),
            test = "Chisq"))

cat("
LEITURA:
")
cat("  A unica celula que rejeita gamma=1 e Poisson + so numericas. E um falso
")
cat("  positivo com DUAS causas somadas: media mal especificada (as categoricas
")
cat("  explicam 188 de deviance em 9 gl) e variancia mal especificada (Poisson
")
cat("  forca var=media). Cada erro sozinho ja inflaria a dispersao; juntos,
")
cat("  estreitam o IC a ponto de inverter a conclusao.

")
cat("  Corrigir QUALQUER UM dos dois ja basta para gamma=1 deixar de ser
")
cat("  rejeitado. A binomial negativa e robusta as duas especificacoes - mais
")
cat("  um argumento para ela, alem do AIC.

")
cat("  => O offset se justifica. Mas so da para afirmar isso depois de acertar
")
cat("     media e variancia. O enunciado pede a Tarefa 1 primeiro; responde-la
")
cat("     nessa ordem, ingenuamente, leva a conclusao errada.
")


sec("FIM - proximo passo: selecao de variaveis (dredge) sobre mod_nb")
