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

# --- DESAFIO C: logistica, ligacao logit (espelha modelo_final do challenge_c.R) --------
# Resultado da selecao backward por LRT a 5%, com cultivar retida por desenho.
FORM_C <- padrao_exportacao ~ cultivar + manejo + altitude_m + idade_lavoura_anos

# Candidatas do modelo amplo do Desafio C. Ficam de fora: os identificadores e o
# esforco amostral (processo de medicao, nao caracteristica do talhao), a contagem
# de brocas (desfecho do Desafio A) e a severidade da ferrugem (desfecho do
# Desafio B, tratada como possivel mediadora na Tarefa 2).
CAND_C <- c("cultivar", "manejo", "irrigacao", "regiao_produtora", "altitude_m",
            "declividade_pct", "idade_lavoura_anos", "densidade_plantio",
            "adubacao_n_kg_ha", "precipitacao_safra_mm", "umidade_relativa_pct",
            "ph_solo", "materia_organica_pct")

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

p_load(statmod, hnp)

# --- a resposta ------------------------------------------------------------------------
num("cExporta",    sum(caf$padrao_exportacao == 1), 0)
num("cNaoExporta", sum(caf$padrao_exportacao == 0), 0)
num("cPrevalencia", 100 * mean(caf$padrao_exportacao), 2)
num("cEPV",        floor(min(table(caf$padrao_exportacao)) / 10), 0)

# --- modelo amplo e selecao backward por LRT -------------------------------------------
# A selecao e refeita aqui, e nao copiada: se o criterio mudar, o numero muda junto.
f_amplo_c <- as.formula(paste("padrao_exportacao ~", paste(CAND_C, collapse = " + ")))
mC_nulo  <- glm(padrao_exportacao ~ 1, data = caf, family = binomial)
mC_amplo <- glm(f_amplo_c, data = caf, family = binomial, na.action = na.fail)
stopifnot(mC_amplo$converged)

num("cNparFull", length(coef(mC_amplo)), 0)
num("cAICfull",  AIC(mC_amplo), 2)
num("cBICfull",  BIC(mC_amplo), 2)

# O modelo amplo explica algo alem da proporcao geral?
an_glob <- anova(mC_nulo, mC_amplo, test = "Chisq")
num("cLRTglobal", an_glob$Deviance[2], 2)
num("cGLglobal",  an_glob$Df[2], 0)
pval("cPglobal",  an_glob[["Pr(>Chi)"]][2])

# Multicolinearidade: GVIF ajustado, comparavel entre fatores e numericas
vif_c <- car::vif(mC_amplo)
gvif_aj <- if (is.matrix(vif_c)) vif_c[, 3] else sqrt(vif_c)
num("cVIFmax", max(gvif_aj), 2)
txt("cVIFmaxTermo", gsub("_", "\\\\_", names(which.max(gvif_aj))))

# Backward por LRT a 5%, mantendo cultivar por desenho (e a exposicao da Tarefa 2)
ALFA_C <- 0.05
m_at <- mC_amplo; retirados <- character(0)
repeat {
    d <- drop1(m_at, test = "Chisq")[-1, ]
    d <- d[!rownames(d) %in% "cultivar", ]
    if (!nrow(d) || max(d[["Pr(>Chi)"]]) < ALFA_C) break
    r_out <- rownames(d)[which.max(d[["Pr(>Chi)"]])]
    retirados <- c(retirados, r_out)
    m_at <- update(m_at, as.formula(paste(". ~ . -", r_out)))
}
# Conferencia: o backward tem de reproduzir a FORM_C declarada em MODELOS FINAIS
if (!setequal(attr(terms(m_at), "term.labels"), attr(terms(FORM_C), "term.labels")))
    stop("O backward nao reproduz FORM_C. Atualize o bloco MODELOS FINAIS.", call. = FALSE)

mC <- glm(FORM_C, data = caf, family = binomial, na.action = na.fail)
stopifnot(mC$converged)
num("cNpar",      length(coef(mC)), 0)
num("cAIC",       AIC(mC), 2)
num("cBIC",       BIC(mC), 2)
num("cRetirados", length(retirados), 0)
txt("cPrimeiroRetirado", gsub("_", "\\\\_", retirados[1]))
txt("cUltimoRetirado",   gsub("_", "\\\\_", retirados[length(retirados)]))

# A simplificacao custou ajuste?
an_fa <- anova(mC, mC_amplo, test = "Chisq")
num("cLRTfinalAmplo", an_fa$Deviance[2], 2)
num("cGLfinalAmplo",  an_fa$Df[2], 0)
pval("cPfinalAmplo",  an_fa[["Pr(>Chi)"]][2])

# Irrelevantes por construcao: p alto no modelo amplo, e saida cedo no backward
d_amplo <- drop1(mC_amplo, test = "Chisq")[-1, ]
irrel <- rownames(d_amplo)[d_amplo[["Pr(>Chi)"]] > 0.5]
num("cIrrelevantesN", length(irrel), 0)
txt("cIrrelevantes", paste(gsub("_", "\\\\_", irrel), collapse = ", "))

# Sensibilidade: a contagem de brocas acrescentaria algo ao modelo final?
caf$log_taxa_broca <- log((caf$n_brocas_capturadas + 0.5) / caf$esforco_amostral)
pval("cBrocaP", anova(mC, update(mC, . ~ . + log_taxa_broca), test = "Chisq")[["Pr(>Chi)"]][2])

# --- ligacao: os dados distinguem logit de probit, cloglog e cauchit? ------------------
aic_lig <- sapply(c("logit", "probit", "cloglog", "cauchit"),
                  function(L) AIC(glm(FORM_C, data = caf, family = binomial(link = L))))
daic_lig <- aic_lig - min(aic_lig)
num("cLigacoesEquiv", sum(daic_lig < 2), 0)
num("cDaicLogit",     daic_lig[["logit"]], 2)
txt("cMelhorLigacao", names(which.min(aic_lig)))

# --- Tarefa 4: razoes de chances com IC de verossimilhanca perfilada -------------------
ic_c <- suppressMessages(confint(mC))
orc <- function(termo, tag, escala = 1) {
    num(paste0("cOR",   tag), exp(escala * coef(mC)[[termo]]), 3)
    num(paste0("cORlo", tag), exp(escala * ic_c[termo, 1]), 3)
    num(paste0("cORhi", tag), exp(escala * ic_c[termo, 2]), 3)
}
for (lv in c("Catuai", "Icatu", "Mundo Novo"))
    orc(paste0("cultivar", lv), gsub(" ", "", lv))
for (lv in c("Integrado", "Organico"))
    orc(paste0("manejo", lv), lv)
orc("altitude_m",         "AltitudeCem", 100)   # por 100 m
orc("idade_lavoura_anos", "IdadeCinco",    5)   # por 5 anos

# --- Tarefa 2: efeito total x efeito direto, agora no modelo final --------------------
# Antes esta comparacao usava um modelo so com a cultivar; agora vem do modelo final,
# de modo que os coeficientes sao os mesmos reportados na Tarefa 4.
cTot <- mC
cDir <- update(mC, . ~ . + severidade_ferrugem)
stopifnot(cDir$converged)
for (lv in c("Catuai", "Icatu", "Mundo Novo")) {
    k <- paste0("cultivar", lv); tag <- gsub(" ", "", lv)
    num(paste0("cTotal", tag),  coef(cTot)[[k]], 4)
    num(paste0("cDireto", tag), coef(cDir)[[k]], 4)
}
num("cSevCoef", coef(cDir)[["severidade_ferrugem"]], 4)
pval("cSevP",   summary(cDir)$coefficients["severidade_ferrugem", 4])
# Razao de chances por 0,10 de severidade: "uma unidade" nao existe nesta base
ic_sev <- suppressMessages(confint(cDir))["severidade_ferrugem", ]
num("cSevORdez",   exp(0.10 * coef(cDir)[["severidade_ferrugem"]]), 3)
num("cSevORdezLo", exp(0.10 * ic_sev[1]), 3)
num("cSevORdezHi", exp(0.10 * ic_sev[2]), 3)
sev_media <- tapply(caf$severidade_ferrugem, caf$cultivar, mean)
num("cSevMediaBourbon", sev_media[["Bourbon"]], 3)
num("cSevMediaIcatu",   sev_media[["Icatu"]], 3)

# --- Tarefa 3: capacidade preditiva e ponto de corte ----------------------------------
p_hat <- fitted(mC)
roc_c <- pROC::roc(caf$padrao_exportacao, p_hat, levels = c(0, 1), direction = "<", quiet = TRUE)
num("cAUC", as.numeric(pROC::auc(roc_c)), 3)
ci_auc <- as.numeric(pROC::ci.auc(roc_c))
num("cAUClo", ci_auc[1], 3); num("cAUChi", ci_auc[3], 3)

# Validacao cruzada em 10 partes: o acerto medido fora da amostra de ajuste
set.seed(20260920)
K <- 10; fold <- sample(rep(1:K, length.out = nrow(caf))); p_cv <- numeric(nrow(caf))
for (k in 1:K) {
    mk <- glm(FORM_C, data = caf[fold != k, ], family = binomial)
    p_cv[fold == k] <- predict(mk, newdata = caf[fold == k, ], type = "response")
}
auc_cv <- as.numeric(pROC::auc(pROC::roc(caf$padrao_exportacao, p_cv,
                                         levels = c(0, 1), direction = "<", quiet = TRUE)))
num("cAUCcv", auc_cv, 3)
num("cFolds", K, 0)

metr <- function(corte) {
    pr <- as.integer(p_hat >= corte); ob <- caf$padrao_exportacao
    VP <- sum(pr == 1 & ob == 1); FP <- sum(pr == 1 & ob == 0)
    FN <- sum(pr == 0 & ob == 1); VN <- sum(pr == 0 & ob == 0)
    c(VP = VP, FP = FP, FN = FN, VN = VN,
      sens = VP / (VP + FN), espec = VN / (VN + FP), acur = (VP + VN) / length(ob))
}
m05 <- metr(0.5)
for (k in c("VP", "FP", "FN", "VN")) num(paste0("c", k), m05[[k]], 0)
num("cSens",     100 * m05[["sens"]], 1)
num("cEspec",    100 * m05[["espec"]], 1)
num("cAcuracia", 100 * m05[["acur"]], 1)
num("cAcuraciaTrivial",
    100 * max(mean(caf$padrao_exportacao), 1 - mean(caf$padrao_exportacao)), 1)

# Corte por custo assimetrico: p > C_FP / (C_FP + C_FN). A razao adotada e uma
# suposicao declarada da equipe, nao um resultado dos dados.
CUSTO_FP <- 2; CUSTO_FN <- 1
corte_c <- CUSTO_FP / (CUSTO_FP + CUSTO_FN)
mCC <- metr(corte_c)
num("cRazaoCusto", CUSTO_FP / CUSTO_FN, 0)
num("cCorteCusto", corte_c, 3)   # 3 casas: o corte e 2/3, e 0,67 arredondaria o registro
num("cSensCusto",  100 * mCC[["sens"]], 1)
num("cEspecCusto", 100 * mCC[["espec"]], 1)
num("cCusto",      (mCC[["FP"]] * CUSTO_FP + mCC[["FN"]] * CUSTO_FN) / nrow(caf), 3)
num("cCustoNunca", sum(caf$padrao_exportacao == 1) * CUSTO_FN / nrow(caf), 3)
num("cCorteYouden", as.numeric(pROC::coords(roc_c, "best", best.method = "youden",
                                            ret = "threshold")[1, 1]), 2)

# --- requisito (c): diagnostico do modelo final ---------------------------------------
set.seed(20260920)
rq <- statmod::qresiduals(mC)
pval("cShapiroP", shapiro.test(rq)$p.value)
env <- hnp::hnp(mC, resid.type = "deviance", sim = 99, conf = 0.95,
                how.many.out = TRUE, print.on = FALSE, plot.sim = FALSE)
num("cEnvFora",  env$out, 0)
num("cEnvTotal", env$total, 0)
num("cPearsonRazao", sum(residuals(mC, type = "pearson")^2) / df.residual(mC), 2)

# Calibracao por decis (Hosmer-Lemeshow): o teste recomendado para dados binarios,
# no lugar da deviance residual
G <- 10
dec <- cut(p_hat, quantile(p_hat, seq(0, 1, length.out = G + 1)), include.lowest = TRUE)
obs <- tapply(caf$padrao_exportacao, dec, sum); nn <- as.vector(table(dec))
pm  <- tapply(p_hat, dec, mean); esp <- nn * pm
hl  <- sum((obs - esp)^2 / (esp * (1 - pm)))
num("cHLqui", hl, 2)
num("cHLgl", G - 2, 0)
pval("cHLp", pchisq(hl, G - 2, lower.tail = FALSE))

# Influencia: identificada pelo id_talhao, como pede o enunciado
cook <- cooks.distance(mC)
num("cCookMax", max(cook), 4)
num("cInfluentesN", sum(cook > 4 / nrow(caf)), 0)
txt("cInfluentesIds", paste(caf$id_talhao[order(-cook)][1:3], collapse = ", "))

# Ligacao (eta^2 como covariavel extra) e escala das continuas (variavel construida)
caf$eta_c <- predict(mC)
pval("cEtaQuadradoP", anova(mC, update(mC, . ~ . + I(eta_c^2)), test = "Chisq")[["Pr(>Chi)"]][2])
for (v in c("altitude_m", "idade_lavoura_anos")) {
    mbt <- update(mC, as.formula(sprintf(". ~ . + I(%s * log(%s))", v, v)))
    tag <- if (v == "altitude_m") "Altitude" else "Idade"
    pval(paste0("cBT", tag), anova(mC, mbt, test = "Chisq")[["Pr(>Chi)"]][2])
}

# --- requisito (d): cenarios na escala da probabilidade -------------------------------
cenario <- function(cv, mj, alt = median(caf$altitude_m), idade = median(caf$idade_lavoura_anos)) {
    nd <- data.frame(cultivar = factor(cv, levels = levels(caf$cultivar)),
                     manejo   = factor(mj, levels = levels(caf$manejo)),
                     altitude_m = alt, idade_lavoura_anos = idade)
    pr <- predict(mC, newdata = nd, type = "link", se.fit = TRUE)
    # unname(): predict() devolve vetor nomeado, e o nome contaminaria o c(p=, lo=, hi=)
    c(p  = unname(plogis(pr$fit)),
      lo = unname(plogis(pr$fit - 1.96 * pr$se.fit)),
      hi = unname(plogis(pr$fit + 1.96 * pr$se.fit)))
}
grade <- expand.grid(cv = levels(caf$cultivar), mj = levels(caf$manejo),
                     stringsAsFactors = FALSE)
grade$p <- mapply(function(cv, mj) cenario(cv, mj)[["p"]], grade$cv, grade$mj)
i_max <- which.max(grade$p); i_min <- which.min(grade$p)
num("cPmax", 100 * grade$p[i_max], 1)
num("cPmin", 100 * grade$p[i_min], 1)
# Os niveis vem do CSV sem acentuacao; o relatorio e em portugues corrente
acentuar <- function(x) {
    mapa <- c(Catuai = "Catuaí", Organico = "orgânico",
              Convencional = "convencional", Integrado = "integrado")
    unname(ifelse(x %in% names(mapa), mapa[x], x))
}
txt("cCenarioMax", sprintf("%s com manejo %s", acentuar(grade$cv[i_max]), acentuar(grade$mj[i_max])))
txt("cCenarioMin", sprintf("%s com manejo %s", acentuar(grade$cv[i_min]), acentuar(grade$mj[i_min])))
num("cAltitudeMediana", median(caf$altitude_m), 0)
num("cIdadeMediana",    median(caf$idade_lavoura_anos), 1)

# Ganho do manejo organico sobre o convencional, por cultivar, em pontos percentuais
ganho <- sapply(levels(caf$cultivar), function(cv)
    100 * (cenario(cv, "Organico")[["p"]] - cenario(cv, "Convencional")[["p"]]))
num("cGanhoManejoMin", min(ganho), 1)
num("cGanhoManejoMax", max(ganho), 1)

# Efeito dos quartis de altitude e de idade, em pontos percentuais
q_alt <- quantile(caf$altitude_m, c(.25, .75))
q_ida <- quantile(caf$idade_lavoura_anos, c(.25, .75))
num("cAltitudePP", 100 * (cenario("Bourbon", "Convencional", alt = q_alt[2])[["p"]] -
                          cenario("Bourbon", "Convencional", alt = q_alt[1])[["p"]]), 1)
num("cIdadePP",    100 * (cenario("Bourbon", "Convencional", idade = q_ida[2])[["p"]] -
                          cenario("Bourbon", "Convencional", idade = q_ida[1])[["p"]]), 1)
txt("cStatus", "final")

cat(sprintf("  selecao: %d -> %d parametros | LRT final x amplo p = %.3f\n",
            length(coef(mC_amplo)), length(coef(mC)), an_fa[["Pr(>Chi)"]][2]))
cat(sprintf("  AUC = %.3f (amostra) | %.3f (validacao cruzada)\n",
            as.numeric(pROC::auc(roc_c)), auc_cv))
cat(sprintf("  envelope: %d de %d fora | Hosmer-Lemeshow p = %.3f\n",
            env$out, env$total, pchisq(hl, G - 2, lower.tail = FALSE)))

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
