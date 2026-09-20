# =====================================================================
# UNIVERSIDADE DO VALE DO ITAJAI - UNIVALI
# ESCOLA POLITECNICA
# PROGRAMA DE POS-GRADUACAO EM COMPUTACAO APLICADA - PPGCA
# MESTRADO EM COMPUTACAO APLICADA
# Disciplina: Modelagem estatistica
# Prof. Dr.: Rodrigo Sant'Ana
# Discentes: Andre Lucas Ribeiro, Leticia Zorzi Rama, Matheus Neis
# Itajai, Santa Catarina, Brasil
#
# =====================================================================
# DESAFIO B - REGRESSAO BETA
# =====================================================================
# Pergunta do corpo tecnico:
# Que fatores agronomicos determinam a severidade da ferrugem, e existe
# uma cultivar comprovadamente mais resistente?
#
# Variavel resposta:
# `severidade_ferrugem` - proporcao continua em (0,1), fracao da area
# foliar lesionada, medida por analise de imagem. NAO existe denominador.
#
# Variaveis candidatas:
# `cultivar`, `manejo`, `irrigacao`, `regiao_produtora`, `altitude_m`,
# `declividade_pct`, `idade_lavoura_anos`, `densidade_plantio`,
# `adubacao_n_kg_ha`, `precipitacao_safra_mm`, `umidade_relativa_pct`,
# `ph_solo`, `materia_organica_pct`
#
# Tarefas do enunciado cobertas neste arquivo:
# 1. Explicitem por que este caso NAO admite modelo binomial nem regressao
#    linear sobre a proporcao, e por que a transformacao arco-seno-raiz
#    (historicamente usada em fitopatologia) e uma solucao inferior a
#    regressao beta. Comparem as duas abordagens empiricamente.
# 2. Ajustem a regressao beta (betareg::betareg), justificando a ligacao
#    escolhida para a media.
#
# Tarefas 3 a 5 (dispersao variavel, efeitos marginais, diagnostico) ficam
# no arquivo seguinte.
#
# O script esta estruturado em:
# - Configurar ambiente de trabalho
# - Carregar a base de dados
# - TAREFA 1.1: por que nao binomial
# - TAREFA 1.2: por que nao regressao linear
# - TAREFA 1.3: por que arco-seno-raiz e inferior
# - TAREFA 1.4: comparacao empirica das tres abordagens
# - TAREFA 2:   a regressao beta e a ligacao logit
#
# Base de dados: data/cafeicultura.csv
# Figuras geradas: challenge_b/figuras/
# =====================================================================


# =====================================================================
# CONFIGURAR AMBIENTE DE TRABALHO
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

p_load(betareg, lmtest, sandwich, car, ggplot2, dplyr, tidyr, patchwork)

# Opcoes gerais
options(scipen = 10, warn = 1)

# Semente de reprodutibilidade.
# OBRIGATORIA sempre que houver simulacao - envelope do hnp, bootstrap,
# residuos quantilicos. Sem ela o resultado muda a cada execucao, e o
# enunciado exige script reprodutivel.
set.seed(20260919)

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
                 "  Abra a PASTA glm-cafeicultura, nao o arquivo solto.", call. = FALSE)
        d <- pai
    }
}

RAIZ <- raiz_repo()
FIG  <- file.path(RAIZ, "challenge_b", "figuras")
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
# BASE DE DADOS
# =====================================================================
secao("BASE DE DADOS")

caf <- read.csv(file.path(RAIZ, "data", "cafeicultura.csv"), stringsAsFactors = TRUE)
y   <- caf$severidade_ferrugem

cat("observacoes:", nrow(caf), "| variaveis:", ncol(caf), "\n")
cat("resposta: severidade_ferrugem\n")
cat("  faixa        :", min(y), "a", max(y), "\n")
cat("  exatamente 0 :", sum(y == 0), "| exatamente 1:", sum(y == 1), "\n")
cat("  distintos    :", length(unique(y)), "em", length(y), "\n")

# Formula das covariaveis, usada por todas as abordagens comparadas
covs <- ~ cultivar + manejo + irrigacao + regiao_produtora +
    umidade_relativa_pct + densidade_plantio + precipitacao_safra_mm +
    altitude_m + idade_lavoura_anos + adubacao_n_kg_ha + ph_solo +
    materia_organica_pct + declividade_pct

parecer("A resposta nao tem observacoes nos extremos 0 ou 1, o que permite ",
        "ajustar a regressao beta diretamente, sem a transformacao de ",
        "Smithson-Verkuilen. O piso de 0,001 e o limite de deteccao do metodo, ",
        "conforme a nota do enunciado.")


# =====================================================================
# TAREFA 1.1 - POR QUE NAO BINOMIAL
# =====================================================================
secao("TAREFA 1.1 - POR QUE ESTE CASO NAO ADMITE MODELO BINOMIAL")

# A binomial modela "k sucessos em n tentativas", e o n entra no ajuste.
# Aqui mediu-se diretamente uma fracao de area foliar por analise de
# imagem: nao houve contagem de folhas doentes sobre um total. Sem n, a
# binomial nao tem como ser escrita.
#
# Isto e literalmente o caso usado pelo professor como exemplo canonico
# em 03_MLGs.pdf, secao 2.3, no box "A CONFUSAO MAIS FREQUENTE DA
# DISCIPLINA":
#
#   "Proporcao continua nao e o mesmo que proporcao binomial. Se voce
#    contou k sucessos em n tentativas, o modelo e binomial e o n entra no
#    ajuste. Se voce mediu diretamente uma fracao, a area foliar lesionada
#    por analise de imagem, um indice de rendimento, nao existe n, e o
#    modelo e beta."

# O argumento tem TRES camadas. A primeira e definicional, a segunda
# corrobora e a terceira e a decisiva.

# ---------------------------------------------------------------------
cat("CAMADA 1 - DEFINICIONAL: o modelo nao pode sequer ser escrito\n")
cat(strrep("-", 72), "\n")
# A binomial e P(K=k) = C(n,k) p^k (1-p)^(n-k). O n e parte da formula,
# nao um detalhe de implementacao. No R, glm(family = binomial) exige
# uma de duas coisas:
#   (a) matriz de duas colunas - sucessos | fracassos
#   (b) a proporcao acompanhada de weights = n
# Nenhuma das duas existe aqui: mediu-se area foliar por analise de
# imagem, sem tentativas a contar.

aviso <- tryCatch({
    glm(severidade_ferrugem ~ cultivar, family = binomial, data = caf)
    "ajustou sem reclamar"
}, warning = function(w) paste("WARNING:", conditionMessage(w)))

cat("   glm(severidade_ferrugem ~ cultivar, family = binomial)\n")
cat("   ->", aviso, "\n")
cat("   O R ajusta, mas avisa: ele nao sabe quantas tentativas houve.\n")

# ---------------------------------------------------------------------
cat("\nCAMADA 2 - EMPIRICA: os dados nao parecem contagem disfarcada\n")
cat(strrep("-", 72), "\n")
# Se y fosse k/n, os valores cairiam sobre uma GRADE de multiplos de 1/n.
cat("   Se y = k/n, os valores cairiam em multiplos de 1/n:\n\n")
for (n_teste in c(10, 20, 50, 100)) {
    grade <- round(y * n_teste) / n_teste
    cat(sprintf("     n = %4d -> y cai exatamente na grade em %5.1f%% dos talhoes\n",
                n_teste, 100 * mean(abs(y - grade) < 1e-9)))
}
cat("\n   valores distintos :", length(unique(y)), "em", length(y), "talhoes\n")
cat("   casas decimais    :", max(nchar(sub("^[0-9]*\\.", "", as.character(y)))), "\n")
cat("   oito menores      :", paste(head(sort(unique(y)), 8), collapse = ", "), "\n")
cat("\n   Isto DESCARTA n pequeno ou moderado. Nao descarta n gigante:\n")
cat("   com n = 10.000, k/n tambem pareceria continuo. Esta camada\n")
cat("   CORROBORA a primeira, nao a substitui.\n")

# ---------------------------------------------------------------------
cat("\nCAMADA 3 - DECISIVA: inventar um n torna a inferencia arbitraria\n")
cat(strrep("-", 72), "\n")
# Suponha que alguem force a binomial inventando um n. A binomial impoe
# Var = mu(1-mu)/n. Trocar o n inventado muda a variancia, logo muda
# erros-padrao, p-valores e intervalos - sem que nenhum DADO mude.
cat("   Forcando a binomial com diferentes n inventados:\n\n")
cat(sprintf("   %-14s %13s %14s %11s %12s\n",
            "n inventado", "coef Icatu", "erro-padrao", "z", "p"))
cat("   ", strrep("-", 68), "\n", sep = "")
eps <- c()
for (n_teste in c(20, 100, 500, 2000)) {
    k  <- round(y * n_teste)
    mm <- suppressWarnings(glm(cbind(k, n_teste - k) ~ cultivar + umidade_relativa_pct,
                               family = binomial, data = caf))
    s <- summary(mm)$coefficients["cultivarIcatu", ]
    eps <- c(eps, s[2])
    cat(sprintf("   %-14d %13.5f %14.5f %11.2f %12s\n",
                n_teste, s[1], s[2], s[3], format.pval(s[4], digits = 2)))
}
mb_ref <- betareg(severidade_ferrugem ~ cultivar + umidade_relativa_pct, data = caf)
sb <- summary(mb_ref)$coefficients$mean["cultivarIcatu", ]
cat("   ", strrep("-", 68), "\n", sep = "")
cat(sprintf("   %-14s %13.5f %14.5f %11.2f %12s\n",
            "beta (sem n)", sb[1], sb[2], sb[3], format.pval(sb[4], digits = 2)))
cat(sprintf("\n   O coeficiente varia pouco; o erro-padrao varia %.0f vezes\n",
            max(eps) / min(eps)))
cat("   entre o menor e o maior n inventado.\n")

parecer("O argumento contra a binomial e estrutural, nao empirico. ",
        "PRIMEIRO, a distribuicao binomial e definida sobre k sucessos em n ",
        "tentativas, e o n entra no ajuste; aqui a severidade foi medida por ",
        "analise de imagem da area foliar, sem tentativas a contar, de modo que ",
        "o modelo nao pode ser escrito sem que se invente um denominador. ",
        "SEGUNDO, os dados corroboram essa leitura: caso fossem k/n com n ",
        "pequeno ou moderado, os valores cairiam sobre uma grade de multiplos ",
        "de 1/n, o que nao se observa - sao ", length(unique(y)), " valores ",
        "distintos entre ", length(y), " talhoes, com quatro casas decimais. ",
        "TERCEIRO, e decisivo: caso se forcasse a binomial com um n arbitrario, ",
        "a variancia imposta pela familia, mu(1-mu)/n, passaria a depender dessa ",
        "escolha. Variando o n inventado de 20 a 2000, o coeficiente da cultivar ",
        "Icatu permanece proximo de -1,30, mas seu erro-padrao varia por um ",
        "fator de aproximadamente ", round(max(eps)/min(eps)), ". A incerteza ",
        "reportada seria, portanto, decisao do analista e nao evidencia dos ",
        "dados - e nao existe criterio estatistico para fixar esse n.")


# =====================================================================
# TAREFA 1.2 - POR QUE NAO REGRESSAO LINEAR
# =====================================================================
secao("TAREFA 1.2 - POR QUE NAO REGRESSAO LINEAR SOBRE A PROPORCAO")

mod_lm <- lm(update(covs, severidade_ferrugem ~ .), data = caf)
verificar(mod_lm, "mod_lm")

# A beta e ajustada aqui apenas como REFERENCIA de comparacao na camada 3.
# Sua justificativa formal esta na Tarefa 2.
mod_beta <- betareg(update(covs, severidade_ferrugem ~ .), data = caf, link = "logit")
verificar(mod_beta, "mod_beta")

pred_lm <- fitted(mod_lm)

# Tambem aqui o argumento tem camadas. As duas primeiras sao SINTOMAS, e
# ambos admitem remendo; a terceira e a causa, e nao admite.

# ---------------------------------------------------------------------
cat("\nSINTOMA 1 - predicoes fora do dominio valido\n")
cat(strrep("-", 72), "\n")
neg <- which(pred_lm < 0)
cat("   predicoes negativas:", length(neg), "de", length(pred_lm), "talhoes\n")
cat("   menor predicao     :", round(min(pred_lm), 4), "\n")
cat("   severidade OBSERVADA nesses talhoes:",
    paste(round(y[neg], 4), collapse = ", "), "\n")
cat("   -> o modelo erra para baixo exatamente onde o piso aperta.\n\n")
cat("   Da para remendar truncando em zero? Da - mas a que custo:\n")
cat(sprintf("     vies do lm sem truncar : %.8f  (exatamente zero, por construcao)\n",
            mean(pred_lm) - mean(y)))
cat(sprintf("     vies do lm truncado    : %.6f\n", mean(pmax(pred_lm, 0)) - mean(y)))
cat("   O remendo e ad hoc, nao vem de teoria nenhuma, e destroi a unica\n")
cat("   propriedade boa do modelo linear.\n")

# ---------------------------------------------------------------------
cat("\nSINTOMA 2 - variancia nao constante (heterocedasticidade)\n")
cat(strrep("-", 72), "\n")
bp <- bptest(mod_lm)
cat("   Breusch-Pagan: BP =", round(bp$statistic, 2), ", gl =", bp$parameter,
    ", p =", format.pval(bp$p.value, digits = 3), "\n\n")
gr <- cut(pred_lm, quantile(pred_lm, seq(0, 1, .25)), include.lowest = TRUE)
cat("   desvio-padrao dos residuos por faixa de predicao:\n")
print(round(tapply(residuals(mod_lm), gr, sd), 5))

cat("\n   Da para remendar com erros-padrao robustos? Da:\n")
cl <- coeftest(mod_lm)["cultivarIcatu", ]
rb <- coeftest(mod_lm, vcov = vcovHC(mod_lm, type = "HC3"))["cultivarIcatu", ]
cat(sprintf("     erro-padrao classico (cultivarIcatu): %.6f  p = %s\n",
            cl[2], format.pval(cl[4], digits = 2)))
cat(sprintf("     erro-padrao robusto  (HC3)          : %.6f  p = %s\n",
            rb[2], format.pval(rb[4], digits = 2)))
cat("   Alguem poderia parar por aqui e declarar o problema resolvido.\n")

# ---------------------------------------------------------------------
cat("\nCAUSA - a variancia de uma proporcao DEPENDE da media\n")
cat(strrep("-", 72), "\n")
# Este e o ponto que nao tem remendo. Perto de 0 e de 1 nao ha espaco
# para variar; no meio da escala ha. Nenhum modelo de variancia constante
# representa isso, por mais robusto que seja o erro-padrao.
mu_b  <- fitted(mod_beta)
phi_b <- mod_beta$coefficients$precision
faixa <- cut(mu_b, quantile(mu_b, seq(0, 1, .2)), include.lowest = TRUE)

var_tab <- do.call(rbind, lapply(split(seq_along(y), faixa), function(ix) {
    m <- mean(mu_b[ix])
    data.frame(n = length(ix), mu_medio = m,
               observada  = var(y[ix]),
               supoe_lm   = summary(mod_lm)$sigma^2,
               supoe_beta = m * (1 - m) / (1 + phi_b))
}))
print(round(var_tab, 6))
cat(sprintf("\n   A variancia observada varia %.1f vezes entre a menor e a maior faixa.\n",
            max(var_tab$observada) / min(var_tab$observada)))
cat(sprintf("   O modelo linear supoe o MESMO valor (%.6f) nas cinco.\n",
            summary(mod_lm)$sigma^2))
cat("   A beta acompanha a subida, porque a variancia dela e funcao da media.\n")

parecer("O modelo linear falha aqui por dois sintomas e uma causa. Os sintomas ",
        "sao predicoes fora do dominio - ", length(neg), " valores negativos, o ",
        "menor deles ", round(min(pred_lm), 4), " - e heterocedasticidade ",
        "(Breusch-Pagan p = ", format.pval(bp$p.value, digits = 3), "). Cabe ",
        "registrar que ambos admitem remendo: as predicoes podem ser truncadas ",
        "em zero e os erros-padrao podem ser estimados de forma robusta. Nenhum ",
        "dos dois, porem, alcanca a causa. A variancia de uma proporcao e ",
        "estruturalmente ligada a sua media, pois nao ha espaco para variar ",
        "perto dos limites do intervalo; medida em cinco faixas de media ",
        "ajustada, a variancia observada varia ",
        sprintf("%.1f", max(var_tab$observada)/min(var_tab$observada)),
        " vezes, enquanto o modelo linear supoe um unico valor para todas. O ",
        "erro-padrao robusto corrige a inferencia sobre um modelo cuja ",
        "estrutura permanece incorreta, ao custo de estimativas menos ",
        "eficientes - e e precisamente essa dependencia que a familia beta ",
        "incorpora por construcao.")

# --- Figura: as tres variancias
vplot <- rbind(
    data.frame(mu = var_tab$mu_medio, variancia = var_tab$observada,
               fonte = "observada nos dados"),
    data.frame(mu = var_tab$mu_medio, variancia = var_tab$supoe_lm,
               fonte = "suposta pelo modelo linear"),
    data.frame(mu = var_tab$mu_medio, variancia = var_tab$supoe_beta,
               fonte = "suposta pela regressao beta"))
vplot$fonte <- factor(vplot$fonte, levels = c("observada nos dados",
                                              "suposta pela regressao beta",
                                              "suposta pelo modelo linear"))
g0 <- ggplot(vplot, aes(x = mu, y = variancia, colour = fonte, linetype = fonte)) +
    geom_line(linewidth = 1) +
    geom_point(data = subset(vplot, fonte == "observada nos dados"), size = 3.2) +
    scale_colour_manual(values = c("observada nos dados" = "#191C17",
                                   "suposta pela regressao beta" = VERDE,
                                   "suposta pelo modelo linear" = VERMELHO)) +
    scale_linetype_manual(values = c("observada nos dados" = "solid",
                                     "suposta pela regressao beta" = "solid",
                                     "suposta pelo modelo linear" = "dashed")) +
    labs(x = "media ajustada do grupo", y = "variancia", colour = NULL, linetype = NULL,
         title = "O modelo linear supoe uma variancia que os dados nao tem",
         subtitle = paste0("a variancia observada varia ",
                           sprintf("%.1f", max(var_tab$observada)/min(var_tab$observada)),
                           "x entre as faixas; o modelo linear supoe uma linha reta"),
         caption = "340 talhoes em 5 faixas de media ajustada") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
salvar(g0, "B1_tres_variancias.png", 9, 6)


# =====================================================================
# TAREFA 1.3 - POR QUE ARCO-SENO-RAIZ E INFERIOR
# =====================================================================
secao("TAREFA 1.3 - POR QUE A TRANSFORMACAO ARCO-SENO-RAIZ E INFERIOR")

# A transformacao classica da fitopatologia: z = asin(sqrt(y)), seguida de
# minimos quadrados.
#
# ATENCAO AO SENTIDO DO ARGUMENTO. Ao testar empiricamente, verificou-se
# que nesta base a transformacao FUNCIONA no que promete. Sustentar que
# ela falha seria afirmar o que os dados negam. O argumento correto e
# outro, e esta construido abaixo em quatro etapas.

caf$z_arco <- asin(sqrt(y))
mod_arco   <- lm(update(covs, z_arco ~ .), data = caf)
verificar(mod_arco, "mod_arco")
pred_arco <- sin(fitted(mod_arco))^2
vies      <- mean(pred_arco) - mean(y)

# ---------------------------------------------------------------------
cat("\nETAPA 1 - ELA FUNCIONA. Comecar reconhecendo isso.\n")
cat(strrep("-", 72), "\n")
bp_y <- bptest(mod_lm); bp_z <- bptest(mod_arco)
cat("   A promessa do arco-seno-raiz e ESTABILIZAR a variancia. E cumpre:\n\n")
cat(sprintf("     Breusch-Pagan no lm sobre y             : BP = %6.2f  p = %s\n",
            bp_y$statistic, format.pval(bp_y$p.value, digits = 3)))
cat(sprintf("     Breusch-Pagan no lm sobre asin(sqrt(y)) : BP = %6.2f  p = %s\n",
            bp_z$statistic, format.pval(bp_z$p.value, digits = 3)))
g_y <- cut(fitted(mod_lm),   quantile(fitted(mod_lm),   seq(0, 1, .25)), include.lowest = TRUE)
g_z <- cut(fitted(mod_arco), quantile(fitted(mod_arco), seq(0, 1, .25)), include.lowest = TRUE)
sd_y <- tapply(residuals(mod_lm), g_y, sd); sd_z <- tapply(residuals(mod_arco), g_z, sd)
cat(sprintf("\n     razao max/min do desvio dos residuos, sobre y      : %.2f\n",
            max(sd_y) / min(sd_y)))
cat(sprintf("     razao max/min do desvio dos residuos, sobre arcsen : %.2f\n",
            max(sd_z) / min(sd_z)))

rs_y <- resettest(mod_lm, power = 2:3, type = "fitted")
rs_z <- resettest(mod_arco, power = 2:3, type = "fitted")
cat("\n   E ainda melhora a linearidade (teste RESET):\n\n")
cat(sprintf("     lm sobre y      : F = %6.3f  p = %s\n",
            rs_y$statistic, format.pval(rs_y$p.value, digits = 3)))
cat(sprintf("     lm sobre arcsen : F = %6.3f  p = %s\n",
            rs_z$statistic, format.pval(rs_z$p.value, digits = 3)))
cat("\n   Portanto o argumento NAO pode ser que ela falha. Ela nao falha aqui.\n")

# ---------------------------------------------------------------------
cat("\nETAPA 2 - MAS FUNCIONA POR SORTE, NAO POR CONSTRUCAO.\n")
cat(strrep("-", 72), "\n")
# A transformacao foi derivada para proporcoes BINOMIAIS: se p ~ Bin(n,p)/n,
# entao Var[asin(sqrt(p))] ~ 1/(4n), constante. Nossa resposta nao e
# binomial, entao a garantia teorica simplesmente nao se aplica.
n_sim <- 50
cat("   Origem teorica: para p ~ Bin(n,p)/n, Var[asin(sqrt(p))] ~ 1/(4n).\n")
cat(sprintf("   Simulando com n = %d, esperado 1/(4n) = %.6f:\n\n", n_sim, 1 / (4 * n_sim)))
cat(sprintf("     %-24s %14s %20s\n", "origem", "var(y)", "var(asin(sqrt(y)))"))
for (p0 in c(.05, .15, .30, .50)) {
    v <- rbinom(20000, n_sim, p0) / n_sim
    cat(sprintf("     binomial p = %.2f        %14.6f %20.6f\n", p0, var(v), var(asin(sqrt(v)))))
}
phi_sim <- mod_beta$coefficients$precision
cat(sprintf("\n   Agora de uma BETA com phi = %.1f, que e o caso real:\n\n", phi_sim))
for (m0 in c(.05, .15, .30, .50)) {
    v <- rbeta(20000, m0 * phi_sim, (1 - m0) * phi_sim)
    cat(sprintf("     beta mu = %.2f           %14.6f %20.6f\n", m0, var(v), var(asin(sqrt(v)))))
}
cat("\n   A estabilizacao e parcial e sem garantia teorica. Funcionou nesta\n")
cat("   base; nada assegura que funcione em outra.\n")

# ---------------------------------------------------------------------
cat("\nETAPA 3 - O QUE ELA COBRA MESMO FUNCIONANDO.\n")
cat(strrep("-", 72), "\n")
cat("   CUSTO 1 - a volta para a escala original e enviesada (Jensen):\n")
cat(sprintf("     media observada de y        : %.6f\n", mean(y)))
cat(sprintf("     media das predicoes (arco)  : %.6f\n", mean(pred_arco)))
cat(sprintf("     vies                        : %.6f  (%.1f%% da media)\n",
            vies, 100 * abs(vies) / mean(y)))
cat(sprintf("     media das predicoes (lm)    : %.6f  <- vies zero por construcao\n",
            mean(pred_lm)))

cat("\n   CUSTO 2 - os coeficientes nao tem leitura agronomica:\n")
ca <- coef(mod_arco)[["cultivarIcatu"]]
cb <- coef(mod_beta)[["cultivarIcatu"]]
cat(sprintf("     beta do Icatu na escala arco-seno : %+.5f\n", ca))
cat("       -> nao e ponto percentual, nao e razao de chances, nao e nada\n")
cat(sprintf("     beta do Icatu na escala logit     : %+.5f\n", cb))
cat(sprintf("       -> exp(beta) = %.4f, razao de chances, converte em pontos pct\n", exp(cb)))

fora_z <- sum(fitted(mod_arco) < 0 | fitted(mod_arco) > pi / 2)
cat("\n   NAO-CUSTO - registrar o que nao se confirmou:\n")
cat(sprintf("     predicoes de z fora de [0, pi/2]: %d\n", fora_z))
cat("     O argumento de que a transformacao tambem prediz fora do intervalo\n")
cat("     NAO se sustenta nesta base. Nao deve ser usado no relatorio.\n")

# ---------------------------------------------------------------------
cat("\nETAPA 4 - O DECISIVO: ela fecha uma pergunta do enunciado.\n")
cat(strrep("-", 72), "\n")
cat("   A Tarefa 3 pede, literalmente:\n")
cat("     'investiguem se o parametro de precisao phi e constante'\n\n")
cat("   Nao ha como perguntar se a dispersao varia dentro de um metodo cujo\n")
cat("   OBJETIVO e fazer a dispersao desaparecer. A transformacao apaga\n")
cat("   exatamente a quantidade que o enunciado manda investigar.\n\n")
cat("   A regressao beta trata a dispersao como parametro e permite modela-la:\n")
cat("     betareg(y ~ x1 + x2 | z1)\n")

parecer("A avaliacao empirica desta base contraria a expectativa usual e merece ",
        "registro explicito: a transformacao arco-seno-raiz CUMPRE o que promete. ",
        "O teste de Breusch-Pagan deixa de rejeitar a homocedasticidade (p = ",
        format.pval(bp_z$p.value, digits = 3), " contra ",
        format.pval(bp_y$p.value, digits = 3), " sobre a resposta original), a razao ",
        "entre os desvios extremos dos residuos cai de ",
        sprintf("%.2f para %.2f", max(sd_y) / min(sd_y), max(sd_z) / min(sd_z)),
        ", e o teste RESET deixa de acusar nao-linearidade. O argumento contra ela, ",
        "portanto, nao pode ser o de que falha. E outro, em tres pontos. PRIMEIRO, a ",
        "estabilizacao carece de fundamento neste caso: a transformacao foi derivada ",
        "para proporcoes binomiais, para as quais a variancia transformada converge a ",
        "1/(4n); a resposta aqui e uma proporcao continua sem denominador, de modo que ",
        "a propriedade vale por coincidencia amostral e nao por construcao. SEGUNDO, a ",
        "volta a escala original e enviesada pela desigualdade de Jensen - o vies ",
        "medido foi de ", sprintf("%.1f%%", 100 * abs(vies) / mean(y)), " da media - e ",
        "os coeficientes estimados na escala transformada nao admitem leitura ",
        "agronomica. TERCEIRO, e decisivo para este trabalho, a transformacao ",
        "inviabiliza a Tarefa 3: nao se pode investigar se o parametro de precisao ",
        "varia dentro de um metodo cujo proposito e eliminar a variacao da dispersao. ",
        "A regressao beta trata a dispersao como parametro do modelo e permite ",
        "modela-la explicitamente.")

# --- Figura: as duas promessas, cumpridas
cmp_prom <- data.frame(
    criterio = rep(c("homocedasticidade\n(Breusch-Pagan)", "linearidade\n(RESET)"), each = 2),
    modelo   = rep(c("lm sobre y", "lm sobre asin(sqrt(y))"), 2),
    p        = c(bp_y$p.value, bp_z$p.value, rs_y$p.value, rs_z$p.value))
cmp_prom$modelo <- factor(cmp_prom$modelo,
                          levels = c("lm sobre y", "lm sobre asin(sqrt(y))"))

g3 <- ggplot(cmp_prom, aes(x = modelo, y = p, fill = modelo)) +
    geom_hline(yintercept = .05, colour = VERMELHO, linetype = "dashed") +
    geom_col(width = .55, colour = "black", linewidth = .3) +
    geom_text(aes(label = format.pval(p, digits = 2)), vjust = -0.6, size = 3.6) +
    facet_wrap(~ criterio) +
    scale_fill_manual(values = c("lm sobre y" = "grey80",
                                 "lm sobre asin(sqrt(y))" = VERDE), guide = "none") +
    scale_y_continuous(limits = c(0, max(cmp_prom$p) * 1.3)) +
    labs(x = NULL, y = "p-valor",
         title = "A transformacao cumpre o que promete nesta base",
         subtitle = "acima da linha tracejada em 0,05 o teste nao acusa problema",
         caption = "por isso o argumento contra ela e outro - ver etapas 2 a 4") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 12, hjust = 1))
salvar(g3, "B1_arco_cumpre_a_promessa.png", 9, 5.5)


# =====================================================================
# TAREFA 1.4 - COMPARACAO EMPIRICA DAS TRES ABORDAGENS
# =====================================================================
secao("TAREFA 1.4 - COMPARACAO EMPIRICA DAS TRES ABORDAGENS")

# O enunciado pede comparacao empirica. Uma advertencia metodologica:
# o AIC NAO pode ser usado aqui. lm(y), lm(asin(sqrt(y))) e betareg tem
# RESPOSTAS diferentes, logo verossimilhancas em escalas diferentes, e os
# valores nao sao comparaveis. A comparacao justa e o erro de predicao na
# escala ORIGINAL, comum as tres.

# mod_beta ja foi ajustado na Tarefa 1.2 como referencia de variancia.

erro <- function(pred, nome) data.frame(
    abordagem       = nome,
    RMSE            = round(sqrt(mean((y - pred)^2)), 5),
    MAE             = round(mean(abs(y - pred)), 5),
    vies            = round(mean(pred) - mean(y), 7),
    fora_do_dominio = sum(pred <= 0 | pred >= 1))

comp <- rbind(
    erro(pred_lm,          "linear sobre y"),
    erro(pred_arco,        "linear sobre asin(sqrt(y))"),
    erro(fitted(mod_beta), "regressao beta"))
print(comp, row.names = FALSE)

parecer("A regressao beta apresenta o menor RMSE (", comp$RMSE[3], " contra ",
        comp$RMSE[1], " do modelo linear), vies praticamente nulo e nenhuma ",
        "predicao fora do dominio. Cabe ser explicito quanto a magnitude: a ",
        "vantagem PREDITIVA e pequena, e nao e nela que o argumento se apoia. ",
        "A justificativa da beta e estrutural - ela e a unica das tres que ",
        "respeita o suporte da resposta por construcao, modela a variancia em ",
        "vez de escondê-la, e produz coeficientes interpretaveis na escala do ",
        "fenomeno.")

# --- Figura: faixa de predicoes contra o dominio valido
faixas <- data.frame(
    abordagem = factor(c("observado", "linear", "arco-seno", "beta"),
                       levels = c("beta", "arco-seno", "linear", "observado")),
    minimo = c(min(y), min(pred_lm), min(pred_arco), min(fitted(mod_beta))),
    maximo = c(max(y), max(pred_lm), max(pred_arco), max(fitted(mod_beta))))
faixas$invalida <- faixas$minimo < 0

g1 <- ggplot(faixas, aes(y = abordagem)) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
             fill = VERMELHO, alpha = .07) +
    geom_vline(xintercept = 0, colour = VERMELHO, linetype = "dashed") +
    geom_segment(aes(x = minimo, xend = maximo, yend = abordagem,
                     colour = invalida), linewidth = 5, lineend = "round") +
    scale_colour_manual(values = c(`FALSE` = VERDE, `TRUE` = VERMELHO),
                        guide = "none") +
    annotate("text", x = 0, y = Inf, label = "dominio impossivel", hjust = 1.05,
             vjust = 1.8, colour = VERMELHO, size = 3.3) +
    labs(x = "severidade predita", y = NULL,
         title = "So a beta nao consegue sair do intervalo valido",
         subtitle = paste0("o modelo linear prediz ", sum(pred_lm < 0),
                           " valores negativos, o menor deles ",
                           round(min(pred_lm), 4)),
         caption = "mesmas covariaveis nas tres abordagens") +
    theme_bw(base_size = 12)
salvar(g1, "B1_faixa_das_predicoes.png", 9, 5)

# --- Figura: residuos do lm contra o ajustado, mostrando o funil
res_lm <- data.frame(ajustado = pred_lm, residuo = residuals(mod_lm))
g2 <- ggplot(res_lm, aes(x = ajustado, y = residuo)) +
    geom_hline(yintercept = 0, colour = "grey55") +
    geom_vline(xintercept = 0, colour = VERMELHO, linetype = "dashed") +
    geom_point(pch = 21, fill = "white", colour = "black", size = 2, alpha = .65) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
                colour = VERDE, linewidth = .9) +
    labs(x = "valor ajustado", y = "residuo",
         title = "O funil dos residuos do modelo linear",
         subtitle = paste0("o espalhamento cresce com a predicao - Breusch-Pagan p = ",
                           format.pval(bp$p.value, digits = 3)),
         caption = "a linha vermelha marca o zero: pontos a esquerda dela sao predicoes impossiveis") +
    theme_bw(base_size = 12)
salvar(g2, "B1_funil_dos_residuos.png", 9, 6)


# =====================================================================
# TAREFA 2 - A REGRESSAO BETA E A LIGACAO
# =====================================================================
secao("TAREFA 2 - AJUSTE DA REGRESSAO BETA E JUSTIFICATIVA DA LIGACAO")

# Parametrizacao de Ferrari & Cribari-Neto (2004), citada pelo professor:
#
#     E[y] = mu          Var[y] = mu (1 - mu) / (1 + phi)
#
# A variancia ja depende da media por construcao - o problema estrutural
# do modelo linear desaparece. O phi e o parametro de precisao: quanto
# maior, mais concentrada a distribuicao em torno de mu.
#
# NOTA: a regressao beta nao e um MLG no sentido estrito de Nelder &
# Wedderburn - sua parametrizacao (mu, phi) nao se encaixa na forma
# canonica com phi separavel. E um MLG estendido, ajustado por maxima
# verossimilhanca direta. (03_MLGs.pdf, nota de rodape 2.) E por isso que
# o DHARMa nao a suporta nativamente, o que importara na Tarefa 5.

cat("LIGACAO DA MEDIA: logit\n\n")
cat("   mu precisa ficar em (0,1); o preditor linear varre toda a reta real.\n")
cat("   logit(mu) = log(mu/(1-mu)) leva (0,1) -> (-inf, +inf), e a inversa\n")
cat("   (a sigmoide) traz de volta. E a mesma logica da ligacao log do\n")
cat("   Desafio A, trocando a restricao 'positivo' por 'entre 0 e 1'.\n")
cat("   (03_MLGs.pdf secao 5.1: a ligacao resolve um descompasso de dominios.)\n\n")

cat("Resultado do ajuste:\n")
cat("   predicoes fora de (0,1):", sum(fitted(mod_beta) <= 0 | fitted(mod_beta) >= 1),
    " <- impossivel por construcao\n")
cat("   faixa das predicoes    :", round(min(fitted(mod_beta)), 4), "a",
    round(max(fitted(mod_beta)), 4), "\n")
cat("   phi estimado (constante):", round(mod_beta$coefficients$precision, 3), "\n")
cat("   log-verossimilhanca     :", round(as.numeric(logLik(mod_beta)), 2), "\n")
cat("   pseudo R2               :", round(mod_beta$pseudo.r.squared, 4), "\n\n")

print(summary(mod_beta))

parecer("Adotou-se a ligacao logit para a media, pelo mesmo criterio que ",
        "determinou a ligacao log no Desafio A: o dominio valido do parametro. ",
        "A media de uma proporcao precisa permanecer em (0,1), enquanto o ",
        "preditor linear percorre toda a reta real; o logit faz a traducao ",
        "entre os dois, de forma monotona e suave. A escolha e tambem o padrao ",
        "do pacote e permite leitura dos efeitos em razoes de chances. O ajuste ",
        "produziu phi constante estimado em ",
        round(mod_beta$coefficients$precision, 2),
        " e nenhuma predicao fora do dominio, o que nao e resultado do ajuste ",
        "e sim consequencia da estrutura do modelo.")

secao("FIM DAS TAREFAS 1 E 2")
cat("Figuras salvas em:", FIG, "\n")
cat("Proximo arquivo: dispersao variavel (Tarefa 3), efeitos marginais\n")
cat("(Tarefa 4) e diagnostico de residuos (Tarefa 5).\n")


# =====================================================================
#                  Creative Commons License 4.0
#                       (CC BY-NC-SA 4.0)
#
#  This is a human-readable summary of (and not a substitute for) the
#  license (https://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
# =====================================================================
