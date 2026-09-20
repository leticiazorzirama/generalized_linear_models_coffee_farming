# =======================================================================================
# UNIVERSIDADE DO VALE DO ITAJAI - UNIVALI
# PROGRAMA DE POS-GRADUACAO EM COMPUTACAO APLICADA - PPGCA
# Disciplina: Modelagem Estatistica | Prof. Dr. Rodrigo Sant'Ana
# Discentes: Matheus Piotroski Neis, Leticia Zorzi Rama, Andre Ribeiro
# =======================================================================================
#
# EXTRATOR DE NUMEROS PARA O RELATORIO
#
# O que este script faz:
#   Le a base, reajusta os modelos FINAIS dos tres desafios e escreve
#   'relatorio/numeros.tex' com um \newcommand por valor.
#
# Por que ele existe:
#   Para que nenhum numero do relatorio seja digitado a mao. Se um modelo
#   mudar, o texto acompanha na proxima compilacao. Numero digitado envelhece
#   em silencio; macro nao.
#
# O que ele NAO exige:
#   Nada dos scripts dos colegas. Este arquivo e autossuficiente: le o CSV e
#   reajusta os modelos por conta propria. Ninguem precisa alterar
#   challenge_a.R nem challenge_c.R.
#
# Uso:
#   Rscript relatorio/numeros.R
#   (ou automaticamente, via relatorio/build.R)
#
# ATENCAO - o unico ponto de manutencao:
#   As formulas da secao "MODELOS FINAIS" abaixo espelham o que cada desafio
#   concluiu. Se um desafio trocar de modelo final, atualize AQUI.
# =======================================================================================

options(scipen = 10, warn = 1)
set.seed(20260920)

# ---------------------------------------------------------------------------------------
# Localizacao da raiz do repositorio - roda de qualquer diretorio, sem setwd()
# ---------------------------------------------------------------------------------------
localizar_script <- function() {
    a <- commandArgs(trailingOnly = FALSE)
    m <- grep("^--file=", a, value = TRUE)
    if (length(m))
        return(normalizePath(sub("^--file=", "", m[1]), winslash = "/", mustWork = FALSE))
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
        if (identical(pai, d))
            stop("Nao encontrei '", marcador, "'. Rode de dentro do repositorio.", call. = FALSE)
        d <- pai
    }
}

RAIZ <- raiz_repo()
SAIDA <- file.path(RAIZ, "relatorio", "numeros.tex")

if (!require(pacman)) { install.packages("pacman", dependencies = TRUE); library(pacman) }
p_load(MASS, betareg, AER, pROC)

# ---------------------------------------------------------------------------------------
# Coletor de macros
# ---------------------------------------------------------------------------------------
MACROS <- new.env(parent = emptyenv())
MACROS$lista <- list()

# Numero formatado no padrao brasileiro (virgula decimal), pronto para o LaTeX
num <- function(nome, valor, casas = 2, milhar = FALSE) {
    txt <- formatC(valor, format = "f", digits = casas,
                   decimal.mark = ",", big.mark = if (milhar) "." else "")
    MACROS$lista[[nome]] <<- trimws(txt)
    invisible(txt)
}

# Texto livre (nomes de cultivar, rotulos)
txt <- function(nome, valor) {
    MACROS$lista[[nome]] <<- as.character(valor)
    invisible(valor)
}

# p-valor com a convencao academica: abaixo do limiar, reporta "< 0,001"
pval <- function(nome, p, limiar = 0.001) {
    MACROS$lista[[nome]] <<- if (p < limiar) {
        paste0("$<$ ", formatC(limiar, format = "f", digits = 3, decimal.mark = ","))
    } else {
        formatC(p, format = "f", digits = 3, decimal.mark = ",")
    }
    invisible(p)
}

cat("\n", strrep("=", 72), "\nEXTRAINDO NUMEROS PARA O RELATORIO\n", strrep("=", 72), "\n\n", sep = "")

# ---------------------------------------------------------------------------------------
# BASE
# ---------------------------------------------------------------------------------------
caf <- read.csv(file.path(RAIZ, "data", "cafeicultura.csv"), stringsAsFactors = TRUE)

num("baseN",         nrow(caf), 0)
num("baseVars",      ncol(caf), 0)   # antes de derivar o esforco - senao conta 20

caf$esforco_amostral <- caf$n_armadilhas * caf$dias_exposicao
num("baseEsforcoMin", min(caf$esforco_amostral), 0)
num("baseEsforcoMax", max(caf$esforco_amostral), 0)
num("baseEsforcoRaz", max(caf$esforco_amostral) / min(caf$esforco_amostral), 0)
cat("  base: n =", nrow(caf), "| esforco varia",
    round(max(caf$esforco_amostral)/min(caf$esforco_amostral)), "vezes\n")

# =======================================================================================
# MODELOS FINAIS - o unico ponto de manutencao deste arquivo
# =======================================================================================

# --- DESAFIO A: binomial negativa com offset composto (espelha mod08) -------------------
FORM_A <- n_brocas_capturadas ~ adubacao_n_kg_ha + altitude_m + densidade_plantio +
    idade_lavoura_anos + materia_organica_pct + umidade_relativa_pct +
    offset(log(esforco_amostral))

# --- DESAFIO B: regressao beta, ligacao logit, phi constante ----------------------------
FORM_B <- severidade_ferrugem ~ cultivar + irrigacao + umidade_relativa_pct +
    densidade_plantio + altitude_m + idade_lavoura_anos + adubacao_n_kg_ha

# --- DESAFIO C: logistica. Ainda nao definido; ver 'Plano de Implementacao - Desafio C' --
FORM_C <- NULL

# =======================================================================================
# DESAFIO A - contagens
# =======================================================================================
cat("\n--- Desafio A ---\n")

# Superdispersao: o argumento que motiva sair da Poisson
mA_pois <- glm(FORM_A, data = caf, family = poisson)
razao_pearson <- sum(residuals(mA_pois, type = "pearson")^2) / df.residual(mA_pois)
num("aPearsonRazao", razao_pearson, 3)
dt <- AER::dispersiontest(mA_pois)
pval("aDispersaoP", dt$p.value)
num("aDispersaoEst", dt$estimate, 3)

# O offset: se o log do esforco entrar como preditor LIVRE, o coeficiente estimado
# deve ser compativel com 1 - e isso que justifica fixa-lo.
#
# CUIDADO: update() nao remove termo de offset. Se a formula for montada com
# update(), o modelo fica com o offset E com o preditor, e o coeficiente passa a
# medir o DESVIO em relacao a 1, nao o gama. A formula e montada a mao.
termos_a <- attr(terms(FORM_A), "term.labels")
termos_a <- termos_a[!grepl("^offset\\(", termos_a)]
mA_livre <- glm(as.formula(paste("n_brocas_capturadas ~",
                paste(c(termos_a, "log(esforco_amostral)"), collapse = " + "))),
                data = caf, family = poisson)
ic_g <- suppressMessages(confint(mA_livre)["log(esforco_amostral)", ])
num("aGamaEst", coef(mA_livre)[["log(esforco_amostral)"]], 3)
num("aGamaIcLo", ic_g[1], 3)
num("aGamaIcHi", ic_g[2], 3)
txt("aGamaContemUm", if (ic_g[1] <= 1 && ic_g[2] >= 1) "contem" else "nao contem")
cat(sprintf("  gama livre = %.3f  IC95%% [%.3f ; %.3f]  -> %s o valor 1\n",
    coef(mA_livre)[["log(esforco_amostral)"]], ic_g[1], ic_g[2],
    if (ic_g[1] <= 1 && ic_g[2] >= 1) "CONTEM" else "NAO CONTEM"))

# Modelo final: binomial negativa
mA <- MASS::glm.nb(FORM_A, data = caf)
stopifnot(mA$converged)
num("aTheta",    mA$theta, 3)
num("aThetaSE",  mA$SE.theta, 3)
num("aAIC",      AIC(mA), 2)
num("aAICpois",  AIC(mA_pois), 2)
num("aNpar",     length(coef(mA)), 0)
razao_nb <- sum(residuals(mA, type = "pearson")^2) / df.residual(mA)
num("aPearsonRazaoNB", razao_nb, 3)

# Efeitos como razoes de taxa (IRR)
for (v in names(coef(mA))[-1]) {
    nome <- gsub("_", "", gsub("_pct|_kg_ha|_m$|_anos", "", v))
    nome <- paste0("aIRR", toupper(substr(nome, 1, 1)), substr(nome, 2, nchar(nome)))
    num(nome, exp(coef(mA)[[v]]), 4)
}
cat("  theta =", round(mA$theta, 3),
    "| razao de Pearson: Poisson", round(razao_pearson, 2), "-> NB", round(razao_nb, 2), "\n")

# =======================================================================================
# DESAFIO B - proporcao continua
# =======================================================================================
cat("\n--- Desafio B ---\n")

mB_full <- betareg(severidade_ferrugem ~ cultivar + manejo + irrigacao + regiao_produtora +
    umidade_relativa_pct + densidade_plantio + precipitacao_safra_mm + altitude_m +
    idade_lavoura_anos + adubacao_n_kg_ha + ph_solo + materia_organica_pct +
    declividade_pct, data = caf, link = "logit")
mB <- betareg(FORM_B, data = caf, link = "logit")
stopifnot(mB$converged, mB_full$converged)

y <- caf$severidade_ferrugem
num("bYmin",      min(y), 4)
num("bYmax",      max(y), 4)
num("bYdistintos", length(unique(y)), 0)
num("bYzeros",    sum(y == 0), 0)
num("bYuns",      sum(y == 1), 0)

num("bPhi",       mB$coefficients$precision, 2)
num("bPhiFull",   mB_full$coefficients$precision, 2)
num("bPseudoR",   mB$pseudo.r.squared, 4)
num("bLogLik",    as.numeric(logLik(mB)), 2)
num("bAIC",       AIC(mB), 2)
num("bAICfull",   AIC(mB_full), 2)
num("bBIC",       BIC(mB), 2)
num("bBICfull",   BIC(mB_full), 2)
num("bNpar",      length(coef(mB, model = "mean")), 0)
num("bNparFull",  length(coef(mB_full, model = "mean")), 0)

# Cultivar em log-odds e em razao de chances, com IC
ic_b <- confint(mB)
for (lv in c("Catuai", "Icatu", "Mundo Novo")) {
    k <- paste0("cultivar", lv)
    tag <- gsub(" ", "", lv)
    num(paste0("bCoef", tag),   coef(mB)[[k]], 4)
    num(paste0("bOR", tag),     exp(coef(mB)[[k]]), 4)
    num(paste0("bORlo", tag),   exp(ic_b[k, 1]), 4)
    num(paste0("bORhi", tag),   exp(ic_b[k, 2]), 4)
    pval(paste0("bP", tag),     summary(mB)$coefficients$mean[k, 4])
}

# Efeito marginal em pontos percentuais: cenario tipico, Bourbon contra Icatu.
# Usa emmeans com a MESMA construcao do challenge_b.R (Tarefa 4). Se aqui fosse
# media sobre os dados observados, o relatorio citaria numero diferente do script.
p_load(emmeans)
grid_tip <- ref_grid(mB, at = list(
    umidade_relativa_pct = median(caf$umidade_relativa_pct),
    densidade_plantio    = median(caf$densidade_plantio)))
em_tip <- emmeans(grid_tip, ~ cultivar, type = "response")
pr <- predict(em_tip); lv <- summary(em_tip)$cultivar
pb  <- pr[lv == "Bourbon"]
pi_ <- pr[lv == "Icatu"]
num("bBourbonPct", 100 * pb, 2)
num("bIcatuPct",   100 * pi_, 2)
num("bDifPP",      100 * (pb - pi_), 2)

# Dispersao variavel: o LRT que o enunciado exige
rhs_b  <- paste(attr(terms(FORM_B), "term.labels"), collapse = " + ")
mB_phi <- betareg(as.formula(paste("severidade_ferrugem ~", rhs_b, "| cultivar")),
                  data = caf, link = "logit")
lr <- 2 * (as.numeric(logLik(mB_phi)) - as.numeric(logLik(mB)))
num("bLRphi",  lr, 3)
pval("bPphi",  pchisq(lr, 3, lower.tail = FALSE))
cat("  phi =", round(mB$coefficients$precision, 2),
    "| pseudo-R2 =", round(mB$pseudo.r.squared, 4),
    "| Bourbon-Icatu =", round(100*(pb - pi_), 2), "p.p.\n")

# --- Comparacao das ligacoes: replica a Tarefa 2 do challenge_b.R ----------------------
# O relatorio cita o dAIC da cauchit e a faixa de pontos percentuais entre ligacoes.
# Sem macro, esses valores envelheceriam em silencio se a base ou o modelo mudassem.
rhs_full <- paste(c("cultivar", "manejo", "irrigacao", "regiao_produtora",
    "umidade_relativa_pct", "densidade_plantio", "precipitacao_safra_mm", "altitude_m",
    "idade_lavoura_anos", "adubacao_n_kg_ha", "ph_solo", "materia_organica_pct",
    "declividade_pct"), collapse = " + ")
f_full <- as.formula(paste("severidade_ferrugem ~", rhs_full))

LIGACOES <- c("logit", "probit", "cloglog", "log", "loglog", "cauchit")
aics <- sapply(LIGACOES, function(L) AIC(betareg(f_full, data = caf, link = L)))
daic <- aics - min(aics)
num("bDaicCauchit", daic[["cauchit"]], 2)
num("bDaicLoglog",  daic[["loglog"]],  2)
num("bDaicLogit",   daic[["logit"]],   2)
num("bLigacoesEquiv", sum(daic < 2), 0)

# Talhao de referencia: medianas nas numericas, primeiro nivel nos fatores
ref <- caf[1, ]
for (v in c("umidade_relativa_pct", "densidade_plantio", "precipitacao_safra_mm",
            "altitude_m", "idade_lavoura_anos", "adubacao_n_kg_ha", "ph_solo",
            "materia_organica_pct", "declividade_pct")) ref[[v]] <- median(caf[[v]])
fx <- function(col, val) factor(val, levels = levels(caf[[col]]))
ref$manejo           <- fx("manejo", "Convencional")
ref$irrigacao        <- fx("irrigacao", "Nao")
ref$regiao_produtora <- fx("regiao_produtora", "Sul de Minas")

dif_pp <- sapply(c("cloglog", "log", "logit", "probit", "loglog"), function(L) {
    m  <- betareg(f_full, data = caf, link = L)
    rb <- ref; rb$cultivar <- fx("cultivar", "Bourbon")
    ri <- ref; ri$cultivar <- fx("cultivar", "Icatu")
    100 * (as.numeric(predict(m, rb)) - as.numeric(predict(m, ri)))
})
num("bPPligacaoMin", min(dif_pp), 1)
num("bPPligacaoMax", max(dif_pp), 1)

# --- Cenarios da Tarefa 4: favoravel e desfavoravel a doenca ---------------------------
cenario_pp <- function(q) {
    g <- ref_grid(mB, at = list(
        umidade_relativa_pct = as.numeric(quantile(caf$umidade_relativa_pct, q)),
        densidade_plantio    = as.numeric(quantile(caf$densidade_plantio, q))))
    e <- emmeans(g, ~ cultivar, type = "response")
    pr <- predict(e); lv <- summary(e)$cultivar
    100 * (pr[lv == "Bourbon"] - pr[lv == "Icatu"])
}
num("bPPfavoravel",    cenario_pp(0.9), 1)
num("bPPdesfavoravel", cenario_pp(0.1), 1)

# --- Estabilidade da dispersao: replica a Etapa 4 da Tarefa 3 --------------------------
r_pad <- residuals(mB, type = "sweighted2")
faixa <- cut(fitted(mB), quantile(fitted(mB), seq(0, 1, .2)), include.lowest = TRUE)
mu_f  <- tapply(fitted(mB), faixa, mean)
var_f <- tapply(y, faixa, var)
sd_f  <- tapply(r_pad, faixa, sd)
num("bRazaoMu",      max(mu_f)  / min(mu_f),  1)
num("bRazaoVar",     max(var_f) / min(var_f), 1)
num("bRazaoSdResid", max(sd_f)  / min(sd_f),  3)
num("bBartlettP",       bartlett.test(r_pad, faixa)$p.value, 3)
num("bBrownForsytheP",
    anova(lm(abs(r_pad - ave(r_pad, faixa, FUN = median)) ~ faixa))[["Pr(>F)"]][1], 3)

# --- IC das razoes de phi entre cultivares --------------------------------------------
ic_phi <- confint(mB_phi)
rn_phi <- rownames(ic_phi)[startsWith(rownames(ic_phi), "(phi)_") &
                           !grepl("Intercept", rownames(ic_phi), fixed = TRUE)]
num("bPhiIcLo", min(exp(ic_phi[rn_phi, 1])), 2)
num("bPhiIcHi", max(exp(ic_phi[rn_phi, 2])), 2)

cat("  ligacoes: ", sum(daic < 2), " dentro de dAIC<2 | cauchit ",
    round(daic[["cauchit"]], 2), " atras\n", sep = "")
cat("  cenarios: desfavoravel ", round(cenario_pp(0.1), 1),
    " pp | tipico ", round(100*(pb - pi_), 1),
    " pp | favoravel ", round(cenario_pp(0.9), 1), " pp\n", sep = "")

# =======================================================================================
# DESAFIO C - binaria
# =======================================================================================
cat("\n--- Desafio C ---\n")

num("cExporta",    sum(caf$padrao_exportacao == 1), 0)
num("cNaoExporta", sum(caf$padrao_exportacao == 0), 0)
num("cPrevalencia", 100 * mean(caf$padrao_exportacao), 2)
num("cEPV",        floor(min(table(caf$padrao_exportacao)) / 10), 0)

# Efeito total x efeito direto - o item autoral da Tarefa 2
cTot <- glm(padrao_exportacao ~ cultivar, data = caf, family = binomial)
cDir <- glm(padrao_exportacao ~ cultivar + severidade_ferrugem, data = caf, family = binomial)
for (lv in c("Catuai", "Icatu", "Mundo Novo")) {
    k <- paste0("cultivar", lv); tag <- gsub(" ", "", lv)
    num(paste0("cTotal", tag),  coef(cTot)[[k]], 4)
    num(paste0("cDireto", tag), coef(cDir)[[k]], 4)
}
num("cSevCoef", coef(cDir)[["severidade_ferrugem"]], 4)
pval("cSevP",   summary(cDir)$coefficients["severidade_ferrugem", 4])

if (is.null(FORM_C)) {
    cat("  MODELO FINAL AINDA NAO DEFINIDO - so os numeros descritivos e de mediacao\n")
    txt("cStatus", "preliminar")
} else {
    mC <- glm(FORM_C, data = caf, family = binomial)
    stopifnot(mC$converged)
    r <- pROC::roc(caf$padrao_exportacao, fitted(mC), quiet = TRUE)
    num("cAUC",    as.numeric(pROC::auc(r)), 4)
    ci <- as.numeric(pROC::ci.auc(r))
    num("cAUClo",  ci[1], 4); num("cAUChi", ci[3], 4)
    txt("cStatus", "final")
}

# =======================================================================================
# ESCRITA DO numeros.tex
# =======================================================================================
dir.create(dirname(SAIDA), showWarnings = FALSE, recursive = TRUE)
linhas <- c(
    "% =====================================================================",
    "% ARQUIVO GERADO AUTOMATICAMENTE - NAO EDITAR A MAO",
    "% Fonte: relatorio/numeros.R  |  Regerar: Rscript relatorio/numeros.R",
    "% Qualquer edicao aqui e perdida na proxima compilacao.",
    "% =====================================================================",
    "",
    sprintf("\\newcommand{\\%s}{%s}", names(MACROS$lista), unlist(MACROS$lista))
)
writeLines(linhas, SAIDA, useBytes = TRUE)

cat("\n", strrep("=", 72), "\n", sep = "")
cat(length(MACROS$lista), "macros escritas em:", SAIDA, "\n")
cat("Use no texto como \\nomeDaMacro - ex: \\bPhi, \\aTheta, \\cPrevalencia\n")
cat(strrep("=", 72), "\n\n")
