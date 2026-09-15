# =====================================================================
# POR QUE 19 PARAMETROS E NAO 13
# Rode e olhe. O R mostra sozinho.
# =====================================================================
caf <- read.csv(file.path(
  if (file.exists("data/cafeicultura.csv")) "." else "..", "data/cafeicultura.csv"))

# ---------------------------------------------------------------------
# O R nao sabe fazer conta com a palavra "Catuai". Ele precisa de numero.
# Entao ele TRADUZ cada fator em colunas de 0 e 1.
# ---------------------------------------------------------------------

cat("O que voce escreveu na formula - uma coluna de texto:\n\n")
print(head(caf[, c("id_talhao", "cultivar")], 6))

cat("\nO que o R construiu para poder calcular:\n\n")
X <- model.matrix(~ cultivar, data = caf)
print(head(X, 6))

cat("\nRepare:\n")
cat("  - 4 cultivares viraram 3 colunas, nao 4\n")
cat("  - Bourbon NAO tem coluna. Ele e o nivel de REFERENCIA.\n")
cat("  - Bourbon e a linha com 0 em todas as tres colunas (TAL003).\n")
cat("  - o intercepto (coluna de 1s) representa o Bourbon.\n")
cat("  - os outros coeficientes medem a DIFERENCA em relacao ao Bourbon.\n")

cat("\nPor que nao 4 colunas? Porque a quarta seria redundante:\n")
cat("sabendo que nao e Catuai, nem Icatu, nem Mundo Novo, so sobra Bourbon.\n")
cat("A quarta coluna nao acrescenta informacao - e por isso que cada fator\n")
cat("gasta (niveis - 1) parametros.\n")

cat("\n--- A conta dos 19 ---\n\n")
f <- n_brocas_capturadas ~ regiao_produtora + cultivar + manejo + irrigacao +
  altitude_m + declividade_pct + idade_lavoura_anos + densidade_plantio +
  adubacao_n_kg_ha + precipitacao_safra_mm + umidade_relativa_pct +
  ph_solo + materia_organica_pct

Xc <- model.matrix(f, data = caf)
cat("colunas que o R realmente criou:", ncol(Xc), "\n\n")
print(colnames(Xc))

cat("\nSomando:\n")
cat("  intercepto                    :  1\n")
cat("  regiao_produtora  (4 niveis)  :  3\n")
cat("  cultivar          (4 niveis)  :  3\n")
cat("  manejo            (3 niveis)  :  2\n")
cat("  irrigacao         (2 niveis)  :  1\n")
cat("  numericas (uma coluna cada)   :  9\n")
cat("                                  --\n")
cat("                                  19\n")
