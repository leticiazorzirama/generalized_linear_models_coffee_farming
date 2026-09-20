# =======================================================================================
# COMPILADOR DO RELATORIO - Windows, Linux e macOS
#
#   Rscript relatorio/build.R
#
# O que faz, em ordem:
#   1. regera relatorio/numeros.tex a partir dos modelos (numeros.R)
#   2. copia as figuras citadas de challenge_*/figuras/ para relatorio/figs/
#   3. compila main.tex -> main.pdf
#
# Por que em R e nao em shell:
#   todo mundo do grupo ja tem R, e o mesmo arquivo roda nos tres sistemas.
#   build.sh e build.ps1 sao atalhos que apenas chamam este aqui.
#
# Primeira vez? Veja relatorio/README.md - e um comando so.
# =======================================================================================

options(warn = 1)

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

DIR  <- dirname(localizar_script())
RAIZ <- dirname(DIR)

passo <- function(n, txt) cat(sprintf("\n[%d/3] %s\n%s\n", n, txt, strrep("-", 72)))
erro  <- function(...) stop("\n\n>> ", paste0(...), "\n", call. = FALSE)

# =======================================================================================
passo(1, "Regerando numeros.tex a partir dos modelos")

r <- system2(file.path(R.home("bin"), "Rscript"),
             shQuote(file.path(DIR, "numeros.R")), stdout = TRUE, stderr = TRUE)
cat(paste(utils::tail(r, 12), collapse = "\n"), "\n")
if (!file.exists(file.path(DIR, "numeros.tex")))
    erro("numeros.R nao produziu numeros.tex. Rode-o sozinho para ver o erro:\n",
         "   Rscript relatorio/numeros.R")

# =======================================================================================
passo(2, "Copiando figuras citadas")

FIGS <- file.path(DIR, "figs")
dir.create(FIGS, showWarnings = FALSE, recursive = TRUE)

# Toda figura que o relatorio cita entra aqui. Se uma figura nao existir, o
# build avisa e segue - o LaTeX e que vai reclamar, com o nome do arquivo.
origens <- c(file.path(RAIZ, "challenge_a", "figuras"),
             file.path(RAIZ, "challenge_b", "figuras"),
             file.path(RAIZ, "challenge_c", "figuras"),
             file.path(RAIZ, "mneis", "figuras_c"))

n_copiadas <- 0
for (o in origens) {
    if (!dir.exists(o)) next
    for (f in list.files(o, pattern = "[.](png|pdf|jpg)$", full.names = TRUE)) {
        file.copy(f, file.path(FIGS, basename(f)), overwrite = TRUE)
        n_copiadas <- n_copiadas + 1
    }
}
cat("   ", n_copiadas, "figuras disponiveis em relatorio/figs/\n")

# =======================================================================================
passo(3, "Compilando main.tex")

antigo <- setwd(DIR)
on.exit(setwd(antigo), add = TRUE)

compilou <- FALSE

# Caminho preferido: tinytex. Instala pacote LaTeX faltante sozinho, o que
# evita a via-crucis de "! LaTeX Error: File xxx.sty not found".
if (requireNamespace("tinytex", quietly = TRUE) && tinytex::is_tinytex()) {
    cat("   motor: TinyTeX (resolve pacotes faltantes automaticamente)\n")
    ok <- try(tinytex::latexmk("main.tex", engine = "pdflatex", bib_engine = "biber"),
              silent = TRUE)
    compilou <- !inherits(ok, "try-error")
}

# Alternativa: latexmk do sistema (TeX Live, MacTeX, MiKTeX com Perl)
if (!compilou && nzchar(Sys.which("latexmk"))) {
    cat("   motor: latexmk do sistema\n")
    compilou <- system2("latexmk",
        c("-pdf", "-interaction=nonstopmode", "-halt-on-error", "main.tex")) == 0
}

# Ultimo recurso: pdflatex + bibtex na mao, o ciclo classico de tres passadas
if (!compilou && nzchar(Sys.which("pdflatex"))) {
    cat("   motor: pdflatex + bibtex (manual)\n")
    system2("pdflatex", c("-interaction=nonstopmode", "-halt-on-error", "main.tex"))
    if (nzchar(Sys.which("biber"))) system2("biber", "main")
    system2("pdflatex", c("-interaction=nonstopmode", "-halt-on-error", "main.tex"))
    compilou <- system2("pdflatex",
        c("-interaction=nonstopmode", "-halt-on-error", "main.tex")) == 0
}

if (!compilou && !file.exists("main.pdf"))
    erro("Nenhum motor LaTeX encontrado ou a compilacao falhou.\n",
         "   Instalacao (uma vez so, igual nos tres sistemas):\n",
         "      Rscript -e \"install.packages('tinytex'); tinytex::install_tinytex()\"\n",
         "   Se compilou antes e quebrou agora, o erro esta em relatorio/main.log")

cat("\n", strrep("=", 72), "\n", sep = "")
if (file.exists("main.pdf")) {
    inf <- file.info("main.pdf")
    cat("OK -> ", normalizePath("main.pdf", winslash = "/"), "\n",
        "     ", round(inf$size / 1024), " KB, gerado em ",
        format(inf$mtime, "%d/%m/%Y %H:%M:%S"), "\n", sep = "")
    cat("\nLembrete: o limite do enunciado e 8 paginas.\n")
} else {
    cat("main.pdf NAO foi gerado. Veja relatorio/main.log\n")
}
cat(strrep("=", 72), "\n")
