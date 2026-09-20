# =======================================================================================
# CONFERIDOR DO RELATORIO
#
#   Rscript relatorio/conferir.R
#
# Verifica, antes de compilar, o que o LaTeX so reclamaria no meio de um log de
# 400 linhas - ou pior, o que ele NAO reclama:
#
#   1. macro usada no texto que nao existe em numeros.tex
#   2. macro definida e nunca usada (candidata a remocao)
#   3. figura citada que nao esta em figs/
#   4. citacao \cite{} sem entrada no refs.bib
#   5. numero solto no texto que deveria ser macro
#   6. quantas marcacoes \pendente{} ainda restam
#
# O item 5 e o que mais importa: e a classe de erro que envelhece em silencio.
# =======================================================================================

localizar_script <- function() {
    a <- commandArgs(trailingOnly = FALSE)
    m <- grep("^--file=", a, value = TRUE)
    if (length(m)) return(normalizePath(sub("^--file=", "", m[1]), winslash = "/", mustWork = FALSE))
    normalizePath(getwd(), winslash = "/")
}
DIR <- dirname(localizar_script())
ler <- function(f) if (file.exists(file.path(DIR, f))) paste(readLines(file.path(DIR, f), warn = FALSE), collapse = "\n") else ""

conteudo <- ler("conteudo.tex")
principal <- ler("main.tex")
numeros  <- ler("numeros.tex")
bib      <- ler("refs.bib")
texto    <- paste(conteudo, principal)

# tira os comentarios LaTeX antes de analisar
texto_sem_com <- gsub("(?m)(?<!\\\\)%.*$", "", texto, perl = TRUE)

problemas <- 0
secao <- function(t) cat("\n", strrep("-", 72), "\n", t, "\n", strrep("-", 72), "\n", sep = "")

# --- 1 e 2: macros ---------------------------------------------------------------------
secao("1. MACROS")
definidas <- unique(regmatches(numeros, gregexpr("(?<=newcommand\\{\\\\)[a-zA-Z]+", numeros, perl = TRUE))[[1]])
usadas_todas <- unique(regmatches(texto_sem_com, gregexpr("\\\\[a-zA-Z]+", texto_sem_com))[[1]])
usadas_todas <- sub("^\\\\", "", usadas_todas)

# so interessam os nomes no padrao das nossas macros: base/a/b/c + maiuscula
padrao <- grepl("^(base|a|b|c)[A-Z]", usadas_todas)
usadas <- unique(usadas_todas[padrao])

faltando <- setdiff(usadas, definidas)
if (length(faltando)) {
    cat("  ERRO - usadas no texto mas NAO definidas:\n")
    for (m in faltando) cat("    \\", m, "\n", sep = "")
    problemas <- problemas + length(faltando)
} else cat("  ok - todas as", length(usadas), "macros usadas existem\n")

ociosas <- setdiff(definidas, usadas)
cat("  ", length(ociosas), "macros definidas e ainda nao usadas\n")

# --- 3: figuras ------------------------------------------------------------------------
secao("3. FIGURAS")
figs <- unique(regmatches(texto_sem_com,
        gregexpr("(?<=includegraphics)(\\[[^]]*\\])?\\{[^}]+\\}", texto_sem_com, perl = TRUE))[[1]])
figs <- gsub(".*\\{|\\}", "", figs)
if (!length(figs)) cat("  nenhuma figura citada\n")
for (f in figs) {
    alvo <- file.path(DIR, "figs", f)
    existe <- file.exists(alvo) ||
              length(Sys.glob(paste0(tools::file_path_sans_ext(alvo), ".*"))) > 0
    cat(if (existe) "  ok   " else "  ERRO ", f, "\n", sep = "")
    if (!existe) problemas <- problemas + 1
}

# --- 4: citacoes -----------------------------------------------------------------------
secao("4. CITACOES")
cits <- regmatches(texto_sem_com, gregexpr("(?<=cite\\{|citeonline\\{)[^}]+", texto_sem_com, perl = TRUE))[[1]]
cits <- unique(trimws(unlist(strsplit(cits, ","))))
chaves <- regmatches(bib, gregexpr("(?<=@)[a-zA-Z]+\\{[^,]+", bib, perl = TRUE))[[1]]
chaves <- sub(".*\\{", "", chaves)
orfas <- setdiff(cits, chaves)
if (length(orfas)) {
    cat("  ERRO - citadas sem entrada no refs.bib:\n")
    for (c in orfas) cat("    ", c, "\n")
    problemas <- problemas + length(orfas)
} else cat("  ok - todas as", length(cits), "citacoes tem entrada\n")
cat("  ", length(setdiff(chaves, cits)), "entradas no .bib ainda nao citadas\n")

# --- 5: numero solto -------------------------------------------------------------------
secao("5. NUMEROS DIGITADOS NO TEXTO")
cat("  Numeros com casa decimal fora de ambiente matematico sao suspeitos:\n")
cat("  ou viram macro, ou precisam de conferencia manual a cada mudanca.\n\n")
linhas <- strsplit(gsub("(?m)(?<!\\\\)%.*$", "", conteudo, perl = TRUE), "\n")[[1]]
susp <- 0
for (i in seq_along(linhas)) {
    l <- linhas[i]
    if (grepl("pendente|includegraphics", l)) next   # nao sao dados do trabalho
    l2 <- gsub("\\$[^$]*\\$", "", l)          # remove matematica inline
    l2 <- gsub("\\\\[a-zA-Z]+", "", l2)        # remove comandos
    m <- regmatches(l2, gregexpr("[0-9]+[,.][0-9]+", l2))[[1]]
    if (length(m)) {
        cat(sprintf("  linha %3d: %s\n", i, paste(m, collapse = ", ")))
        susp <- susp + length(m)
    }
}
if (!susp) cat("  ok - nenhum numero decimal digitado direto\n")

# --- 6: pendencias ---------------------------------------------------------------------
secao("6. PENDENCIAS")
np <- length(gregexpr("\\\\pendente\\{", conteudo)[[1]])
if (!grepl("\\\\pendente\\{", conteudo)) np <- 0
cat("  ", np, "marcacoes \\pendente{} no conteudo\n")
if (np > 0) cat("  (precisam sumir antes da entrega final)\n")

# --- veredicto -------------------------------------------------------------------------
cat("\n", strrep("=", 72), "\n", sep = "")
if (problemas == 0) {
    cat("SEM ERROS BLOQUEANTES.", susp, "numeros digitados e", np, "pendencias a resolver.\n")
} else {
    cat("!!", problemas, "PROBLEMA(S) que quebram a compilacao ou a referencia.\n")
}
cat(strrep("=", 72), "\n")
quit(status = if (problemas > 0) 1 else 0)
