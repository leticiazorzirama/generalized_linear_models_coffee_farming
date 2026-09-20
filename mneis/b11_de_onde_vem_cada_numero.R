# =====================================================================
# DESAFIO B - DE ONDE VEM CADA NUMERO
# Do coeficiente ate a severidade em pontos percentuais, sem pular passo.
# Tudo conferivel: cada linha mostra o valor usado, o coeficiente e o
# produto, e a soma bate com o que o predict() devolve.
# =====================================================================
if (!require(pacman)) { install.packages("pacman"); library(pacman) }
p_load(betareg)

RAIZ <- if (file.exists("data/cafeicultura.csv")) "." else ".."
caf  <- read.csv(file.path(RAIZ, "data", "cafeicultura.csv"), stringsAsFactors = TRUE)

covs <- ~ cultivar + manejo + irrigacao + regiao_produtora + umidade_relativa_pct +
  densidade_plantio + precipitacao_safra_mm + altitude_m + idade_lavoura_anos +
  adubacao_n_kg_ha + ph_solo + materia_organica_pct + declividade_pct
m <- betareg(update(covs, severidade_ferrugem ~ .), data = caf, link = "logit")
b <- coef(m)

titulo <- function(t) cat("\n", strrep("=", 74), "\n", t, "\n", strrep("=", 74), "\n", sep = "")

# ---------------------------------------------------------------------
titulo("1. TODOS OS COEFICIENTES DO MODELO")
print(round(b, 6))

# ---------------------------------------------------------------------
titulo("2. O TALHAO DE REFERENCIA - todos os valores")
nums <- c("umidade_relativa_pct", "densidade_plantio", "precipitacao_safra_mm",
          "altitude_m", "idade_lavoura_anos", "adubacao_n_kg_ha", "ph_solo",
          "materia_organica_pct", "declividade_pct")
ref <- caf[1, ]
for (v in nums) ref[[v]] <- median(caf[[v]])
fx <- function(col, val) factor(val, levels = levels(caf[[col]]))
ref$cultivar         <- fx("cultivar", "Bourbon")
ref$manejo           <- fx("manejo", "Convencional")
ref$irrigacao        <- fx("irrigacao", "Nao")
ref$regiao_produtora <- fx("regiao_produtora", "Sul de Minas")

cat("categoricas - escolhi o nivel de referencia de cada fator:\n")
for (v in c("cultivar", "manejo", "irrigacao", "regiao_produtora"))
  cat(sprintf("  %-22s = %s\n", v, as.character(ref[[v]])))
cat("\nnumericas - usei a MEDIANA de cada uma:\n")
for (v in nums) cat(sprintf("  %-22s = %10.2f\n", v, ref[[v]]))

# ---------------------------------------------------------------------
titulo("3. O PREDITOR LINEAR eta, TERMO A TERMO (cultivar = Bourbon)")
X <- model.matrix(delete.response(terms(m)), data = ref)
cat(sprintf("%-32s %12s %13s %14s\n", "termo", "valor de x", "coef beta", "x * beta"))
cat(strrep("-", 74), "\n")
eta <- 0
for (nm in colnames(X)) {
  contrib <- X[1, nm] * b[[nm]]
  eta <- eta + contrib
  cat(sprintf("%-32s %12.4f %13.6f %14.6f\n", nm, X[1, nm], b[[nm]], contrib))
}
cat(strrep("-", 74), "\n")
cat(sprintf("%-32s %12s %13s %14.6f\n", "SOMA = eta", "", "", eta))
cat("\nconferindo com predict(type='link') :",
    round(as.numeric(predict(m, ref, type = "link")), 6), "\n")

# ---------------------------------------------------------------------
titulo("4. DE eta PARA A SEVERIDADE - o inverso do logit")
cat("   mu = 1 / (1 + exp(-eta))\n\n")
cat(sprintf("   eta            = %.6f\n", eta))
cat(sprintf("   -eta           = %.6f\n", -eta))
cat(sprintf("   exp(-eta)      = %.6f\n", exp(-eta)))
cat(sprintf("   1 + exp(-eta)  = %.6f\n", 1 + exp(-eta)))
cat(sprintf("   mu = 1 / %.6f = %.6f\n", 1 + exp(-eta), 1 / (1 + exp(-eta))))
cat(sprintf("   em pontos percentuais: %.4f%%\n", 100 / (1 + exp(-eta))))
cat("\nconferindo com predict(type='response'):",
    round(as.numeric(predict(m, ref, type = "response")), 6), "\n")

# ---------------------------------------------------------------------
titulo("5. AS QUATRO CULTIVARES - so o termo da cultivar muda")
cat(sprintf("%-12s %13s %14s %13s %13s\n", "cultivar", "coef", "eta", "mu", "pontos %"))
cat(strrep("-", 70), "\n")
mus <- c()
for (cv in levels(caf$cultivar)) {
  r <- ref; r$cultivar <- fx("cultivar", cv)
  e  <- as.numeric(predict(m, r, type = "link"))
  mu <- 1 / (1 + exp(-e)); mus[cv] <- mu
  co <- if (cv == "Bourbon") 0 else b[[paste0("cultivar", cv)]]
  cat(sprintf("%-12s %13.6f %14.6f %13.6f %13.4f\n", cv, co, e, mu, 100 * mu))
}

# ---------------------------------------------------------------------
titulo("6. DE ONDE SAEM OS PONTOS PERCENTUAIS")
mb <- mus[["Bourbon"]]; mi <- mus[["Icatu"]]
cat(sprintf("   severidade do Bourbon  = %.6f   ->  %.4f%%\n", mb, 100 * mb))
cat(sprintf("   severidade do Icatu    = %.6f   ->  %.4f%%\n", mi, 100 * mi))
cat("  ", strrep("-", 56), "\n")
cat(sprintf("   diferenca              = %.6f   ->  %.4f pontos percentuais\n",
            mb - mi, 100 * (mb - mi)))
cat(sprintf("   razao Icatu / Bourbon  = %.6f   ->  Icatu tem %.2f%% da severidade\n",
            mi / mb, 100 * mi / mb))
cat(sprintf("   razao Bourbon / Icatu  = %.6f   ->  Bourbon tem %.2fx a do Icatu\n",
            mb / mi, mb / mi))
cat("\n   ATENCAO: 'ponto percentual' e a diferenca absoluta entre duas\n")
cat("   porcentagens. De 16,97%% para 5,27%% sao 11,70 PONTOS percentuais,\n")
cat("   e nao 11,70%% - que seria outra conta.\n")

# ---------------------------------------------------------------------
titulo("7. OS TRES MANEJOS - so o termo do manejo muda")
ref2 <- ref; ref2$cultivar <- fx("cultivar", "Catuai")
cat("(cultivar fixada em Catuai, para nao usar um extremo)\n\n")
cat(sprintf("%-14s %13s %14s %13s %13s\n", "manejo", "coef", "eta", "mu", "pontos %"))
cat(strrep("-", 72), "\n")
mus2 <- c()
for (mj in levels(caf$manejo)) {
  r <- ref2; r$manejo <- fx("manejo", mj)
  e  <- as.numeric(predict(m, r, type = "link"))
  mu <- 1 / (1 + exp(-e)); mus2[mj] <- mu
  co <- if (mj == "Convencional") 0 else b[[paste0("manejo", mj)]]
  cat(sprintf("%-14s %13.6f %14.6f %13.6f %13.4f\n", mj, co, e, mu, 100 * mu))
}
cat(strrep("-", 72), "\n")
cat(sprintf("   maior - menor = %.6f  ->  %.4f pontos percentuais\n",
            max(mus2) - min(mus2), 100 * (max(mus2) - min(mus2))))

# ---------------------------------------------------------------------
titulo("8. A COMPARACAO")
cat(sprintf("   cultivar (Bourbon -> Icatu)  : %.4f pontos percentuais\n", 100 * (mb - mi)))
cat(sprintf("   manejo   (pior -> melhor)    : %.4f pontos percentuais\n",
            100 * (max(mus2) - min(mus2))))
cat(sprintf("   %.6f / %.6f = %.2f\n", mb - mi, max(mus2) - min(mus2),
            (mb - mi) / (max(mus2) - min(mus2))))
cat("\n   A cultivar move", round((mb - mi) / (max(mus2) - min(mus2))),
    "vezes mais que o manejo.\n")

# ---------------------------------------------------------------------
titulo("9. RESSALVA")
cat("   Estes numeros valem PARA ESTE TALHAO DE REFERENCIA. Como a ligacao\n")
cat("   logit e nao linear, a diferenca em pontos percentuais muda conforme\n")
cat("   o nivel das demais variaveis. Exemplo, variando so a umidade:\n\n")
cat(sprintf("   %-12s %14s %14s %16s\n", "umidade", "Bourbon %", "Icatu %", "diferenca (pp)"))
cat("  ", strrep("-", 58), "\n")
for (u in c(60, 70, 74.1, 80, 90)) {
  r1 <- ref; r1$umidade_relativa_pct <- u; r1$cultivar <- fx("cultivar", "Bourbon")
  r2 <- ref; r2$umidade_relativa_pct <- u; r2$cultivar <- fx("cultivar", "Icatu")
  p1 <- as.numeric(predict(m, r1, type = "response"))
  p2 <- as.numeric(predict(m, r2, type = "response"))
  cat(sprintf("   %-12.1f %14.4f %14.4f %16.4f\n", u, 100 * p1, 100 * p2, 100 * (p1 - p2)))
}
cat("\n   E por isso que a Tarefa 4 pede efeitos marginais em CENARIOS\n")
cat("   agronomicos concretos, e nao um numero unico.\n")
