# =======================================================================================
# UNIVERSIDADE DO VALE DO ITAJAI - UNIVALI
# ESCOLA POLITECNICA
# PROGRAMA DE POS-GRADUACAO EM COMPUTACAO APLICADA - PPGCA
# MESTRADO EM COMPUTACAO APLICADA
# Disciplina: Modelagem estatistica
# Prof. Dr.: Rodrigo Sant'Ana
# Discentes: Andre Lucas Ribeiro, Leticia Zorzi Rama, Matheus Neis
# Itajai, Santa Catarina, Brasil
#
# =======================================================================================
# DESAFIO C - REGRESSAO LOGISTICA (BINOMIAL)
# =====================================================================
# Pergunta do corpo tecnico:
# O que determina a chance de um lote atingir o padrao de exportacao?
#
# Variavel resposta:
# `padrao_exportacao` - binaria: 1 se o lote atingiu o padrao de
# exportacao, 0 caso contrario. Cada talhao contribui com UMA observacao
# de Bernoulli (n = 1), nao com "k lotes de n".
#
# Variaveis candidatas:
# `cultivar`, `manejo`, `irrigacao`, `regiao_produtora`, `altitude_m`,
# `declividade_pct`, `idade_lavoura_anos`, `densidade_plantio`,
# `adubacao_n_kg_ha`, `precipitacao_safra_mm`, `umidade_relativa_pct`,
# `ph_solo`, `materia_organica_pct`
#
# Deixadas de fora do conjunto de candidatas, com justificativa:
# - `id_talhao`: identificador, nao covariavel.
# - `n_armadilhas`, `dias_exposicao`: esforco amostral da broca, sem
#   relacao causal plausivel com a qualidade do lote.
# - `n_brocas_capturadas`: e desfecho do Desafio A e esta confundido pelo
#   esforco amostral; e tratada como verificacao de sensibilidade apenas.
# - `severidade_ferrugem`: e desfecho do Desafio B e MEDIADORA do efeito
#   da cultivar. Entra apenas na Tarefa 2, de forma deliberada.
#
# Tarefas do enunciado:
# 1. Ajustem o modelo logistico e conduzam selecao de variaveis.
# 2. Questao de raciocinio causal: `severidade_ferrugem` e ela propria um
#    desfecho (Desafio B). Ajustem o modelo COM e SEM essa variavel e
#    comparem os coeficientes de `cultivar`. Qual dos dois responde melhor
#    a pergunta do produtor que precisa ESCOLHER a cultivar antes do plantio?
# 3. Avaliem a capacidade preditiva: matriz de confusao, curva ROC e AUC,
#    discutindo a escolha do ponto de corte em funcao do custo assimetrico.
# 4. Interpretem em razoes de chances com intervalos de confianca.
#
# O script esta estruturado em steps cumulativos:
# - STEP 0: configurar ambiente, carregar dados, natureza da resposta,
#           modelo nulo
# - STEP 1: exploracao condicional e logit empirico
# - STEP 2: especificacao - por que binomial e por que logit
# - STEP 3: modelo amplo, deviance, LRT, multicolinearidade
# - STEP 4: selecao de variaveis
# - STEP 5: efeito total x direto (severidade como mediadora)
# - STEP 6: razoes de chances com IC
# - STEP 7: capacidade preditiva e ponto de corte
# - STEP 8: diagnostico de residuos e influencia
# - STEP 9: cenarios e recomendacao
#
# Base de dados: data/cafeicultura.csv
# Figuras geradas: challenge_c/figuras/
# =====================================================================


# =====================================================================
# ---- STEP 0 ---- CONFIGURAR AMBIENTE DE TRABALHO
# =====================================================================

# Pacotes - pacman verifica se ja estao instalados e instala se preciso
if (!require(pacman)) {
    print(paste0("Pacote ainda nao instalado. Instalando..."))
    install.packages("pacman", dependencies = TRUE)
    library(pacman)
} else {
    print(paste0("Pacote ja instalado e carregado"))
    library(pacman)
}

# car     : vif()
# MuMIn   : dredge()
# pROC    : roc(), auc()
# hnp     : envelope simulado meio-normal
# statmod : qresiduals() - residuos quantilicos aleatorizados
# mgcv    : gam() para checar linearidade
# lmtest  : lrtest()
p_load(car, MuMIn, pROC, hnp, statmod, mgcv, lmtest,
       ggplot2, dplyr, tidyr, patchwork)

# Opcoes gerais
options(scipen = 10, warn = 1)

# Semente de reprodutibilidade.
# OBRIGATORIA sempre que houver simulacao ou aleatorizacao: envelope do
# hnp, residuos quantilicos aleatorizados, validacao cruzada. Sem ela o
# resultado muda a cada execucao, e o enunciado exige script reprodutivel.
set.seed(20260920)

# ---------------------------------------------------------------------
# Localizacao da raiz do repositorio.
# O script encontra os dados sozinho, nao importa de onde seja executado:
# pasta aberta como projeto, arquivo solto com Source, Ctrl+Enter linha a
# linha, ou Rscript. Nao depende de setwd().
# ---------------------------------------------------------------------
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
            stop("Nao encontrei '", marcador, "'.\n",
                 "  Abra a PASTA do repositorio, nao o arquivo solto.", call. = FALSE)
        d <- pai
    }
}

RAIZ <- raiz_repo()
FIG  <- file.path(RAIZ, "challenge_c", "figuras")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------
# Funcoes auxiliares
# ---------------------------------------------------------------------

# Cabecalho de secao
secao <- function(texto) {
    barra <- paste(rep("=", 72), collapse = "")
    cat("\n", barra, "\n", texto, "\n", barra, "\n\n", sep = "")
}

# Parecer em prosa. Marcado de forma uniforme para poder ser extraido de
# uma vez so na hora de montar o relatorio:
#   grep -A20 ">> PARECER" saida.txt
parecer <- function(...) {
    txt <- paste0(...)
    cat("\n>> PARECER\n")
    cat(strwrap(txt, width = 72, prefix = "   "), sep = "\n")
    cat("\n")
}

# Salvar figura em disco E mostrar no painel de Plots
salvar <- function(g, nome, w = 9, h = 6) {
    print(g)
    ggsave(file.path(FIG, nome), g, width = w, height = h, dpi = 150, bg = "white")
    cat("   figura salva: ", nome, "\n", sep = "")
}

# Verificacao de ajuste.
# Um glm() que nao converge devolve um objeto de classe glm cujo summary()
# imprime coeficientes e p-valores normalmente - a falha e silenciosa.
# Esta funcao torna a falha ruidosa.
verificar <- function(modelo, nome) {
    conv <- if (!is.null(modelo$converged)) modelo$converged else TRUE
    if (!isTRUE(conv))
        stop("O modelo '", nome, "' NAO convergiu. Nao interprete nada dele.",
             call. = FALSE)
    cat("   [ok] ", nome, " convergiu", sep = "")
    if (!is.null(modelo$iter)) cat(" em ", modelo$iter, " iteracoes", sep = "")
    cat("\n")
    invisible(modelo)
}

# Cores usadas nos graficos
VERDE <- "#18675A"; VERMELHO <- "#A3352A"


# =====================================================================
# ---- STEP 0 ---- BASE DE DADOS E NATUREZA DA RESPOSTA
# =====================================================================
secao("BASE DE DADOS")

caf <- read.csv(file.path(RAIZ, "data", "cafeicultura.csv"), stringsAsFactors = TRUE)
y   <- caf$padrao_exportacao

cat("observacoes:", nrow(caf), "| variaveis:", ncol(caf), "\n")
cat("resposta: padrao_exportacao\n")
cat("  valores distintos :", paste(sort(unique(y)), collapse = ", "), "\n")

# A resposta e 0/1, entao a MEDIA e a proporcao de sucessos: a soma conta
# os 1s, e dividir por n da a fracao. Esse e o estimador de maxima
# verossimilhanca de p numa Bernoulli - o mesmo numero que o modelo nulo
# vai devolver (na escala logit) mais abaixo.
tab_y <- table(y)
p_obs <- mean(y)
cat("  nao atingiu (0)   :", tab_y[["0"]], "\n")
cat("  atingiu (1)       :", tab_y[["1"]], "\n")
cat(sprintf("  proporcao de 1    : %.4f\n", p_obs))

# Variancia de uma Bernoulli: Var(Y) = p(1-p). E maxima em p = 0,5
# (vale 0,25) e vai a zero nos extremos. Quanto mais perto de 0,5, mais
# "informativa" e cada observacao, no sentido de que ha mais incerteza a
# ser explicada. Isso importa para o Step 2: e essa dependencia da
# variancia em relacao a media que o modelo linear ignora.
cat(sprintf("  Var(Y) = p(1-p)   : %.4f  (maximo possivel: 0,25 em p = 0,5)\n", p_obs * (1 - p_obs)))

# ---------------------------------------------------------------------
# Niveis de referencia dos fatores.
# read.csv(stringsAsFactors = TRUE) ordena os niveis ALFABETICAMENTE e o
# primeiro vira a referencia. Todo coeficiente de fator que aparecer nos
# proximos steps e uma comparacao CONTRA esse nivel. Anotar aqui evita
# ler "cultivarIcatu = -0,75" sem saber "-0,75 em relacao a quem".
#
# Referencias resultantes:
#   cultivar         -> Bourbon
#   manejo           -> Convencional
#   irrigacao        -> Nao
#   regiao_produtora -> Cerrado Mineiro
# ---------------------------------------------------------------------
cat("\nniveis dos fatores (o primeiro e a referencia):\n")
fatores <- names(caf)[sapply(caf, is.factor) & names(caf) != "id_talhao"]
for (f in fatores)
    cat(sprintf("  %-17s: %s\n", f, paste(levels(caf[[f]]), collapse = " | ")))

# Formula das covariaveis candidatas, usada em todos os steps seguintes.
# Definida aqui, uma unica vez, para que todos os modelos partam do mesmo
# conjunto e nenhuma covariavel "entre pela janela" mais adiante.
covs <- ~ cultivar + manejo + irrigacao + regiao_produtora +
    altitude_m + declividade_pct + idade_lavoura_anos + densidade_plantio +
    adubacao_n_kg_ha + precipitacao_safra_mm + umidade_relativa_pct +
    ph_solo + materia_organica_pct

cat("\ncovariaveis candidatas:", length(attr(terms(covs), "term.labels")), "\n")


# =====================================================================
# ---- STEP 0 ---- MODELO NULO
# =====================================================================
secao("MODELO NULO")

# O modelo nulo nao tem covariavel: logit(p_i) = beta_0 para todo talhao.
# Serve para tres coisas:
#   (1) demonstrar o que a ligacao logit faz com o intercepto;
#   (2) fornecer a deviance de referencia para o LRT do Step 3;
#   (3) ser a "linha de base" de qualquer medida de ajuste.
mod_nulo <- glm(padrao_exportacao ~ 1,
                family = binomial(link = "logit"),
                data   = caf)
verificar(mod_nulo, "mod_nulo")

# Se a ligacao e logit, o intercepto do nulo tem de ser log(p/(1-p)) da
# proporcao observada: e a unica maneira de a media ajustada coincidir
# com a media amostral, que e o que a maxima verossimilhanca exige para
# a ligacao CANONICA. qlogis() e exatamente log(p/(1-p)).
b0        <- unname(coef(mod_nulo))
logit_obs <- qlogis(p_obs)

cat(sprintf("   intercepto do modelo nulo : %.6f\n", b0))
cat(sprintf("   qlogis(proporcao observada): %.6f\n", logit_obs))
# all.equal() e nao ==: o ajuste e iterativo (IRLS) e para num criterio
# de tolerancia, entao os dois numeros coincidem ate a 6a-8a casa, nao
# bit a bit. Comparar com == poderia dar FALSE por 1e-10 de diferenca.
cat("   coincidem (all.equal):", isTRUE(all.equal(b0, logit_obs, tolerance = 1e-6)), "\n")
# A volta: plogis(intercepto) devolve a proporcao.
cat(sprintf("   plogis(intercepto)        : %.4f  = proporcao observada\n", plogis(b0)))

# Deviance nula: -2 * logLik do modelo so com intercepto. Os graus de
# liberdade sao n - 1 porque um parametro (beta_0) foi estimado. Este e
# o "quanto falta explicar" antes de qualquer covariavel.
dev_nulo <- deviance(mod_nulo)
gl_nulo  <- df.residual(mod_nulo)
cat(sprintf("\n   deviance nula : %.4f  com %d graus de liberdade\n", dev_nulo, gl_nulo))
cat(sprintf("   -2 * logLik   : %.4f  (deve ser identico)\n", -2 * as.numeric(logLik(mod_nulo))))

parecer("A resposta `padrao_exportacao` e binaria: cada um dos ", nrow(caf),
        " talhoes contribui com uma unica observacao de Bernoulli, e nao com ",
        "uma contagem de sucessos em n tentativas. ", tab_y[["1"]], " lotes (",
        sprintf("%.1f%%", 100 * p_obs), ") atingiram o padrao de exportacao, de modo que ",
        "as classes estao aproximadamente balanceadas e a variancia p(1-p) = ",
        sprintf("%.3f", p_obs * (1 - p_obs)), " fica proxima do maximo teorico de 0,25. ",
        "O modelo nulo com ligacao logit devolve intercepto ", sprintf("%.4f", b0),
        ", identico a log(p/(1-p)) da proporcao observada, o que ilustra a ",
        "propriedade da ligacao canonica: a media ajustada reproduz a media ",
        "amostral. A deviance nula de ", sprintf("%.2f", dev_nulo), " com ",
        gl_nulo, " graus de liberdade e a referencia contra a qual todo modelo ",
        "com covariaveis sera testado.")


# =======================================================================================
#                  Creative Commons License 4.0
#                       (CC BY-NC-SA 4.0)
#
#  This is a human-readable summary of (and not a substitute for) the
#  license (https://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
# =======================================================================================
