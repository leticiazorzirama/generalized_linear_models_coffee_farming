# =======================================================================================
# PREPARACAO DO AMBIENTE PARA COMPILAR O RELATORIO
#
#   Rscript relatorio/setup.R
#
# Instala tudo que falta e nao mexe no que ja existe. Pode rodar quantas vezes
# quiser: e idempotente.
#
#   1. confere a versao do R
#   2. confere acesso a internet
#   3. instala os pacotes R usados pelos modelos e pelo build
#   4. instala o TinyTeX, se nao houver nenhum LaTeX na maquina
#   5. instala os pacotes LaTeX que o main.tex exige
#   6. compila um documento de teste, para provar que funciona de verdade
#
# Se voce nem tem o R instalado, use antes o bootstrap do seu sistema:
#   Linux/macOS : ./relatorio/setup.sh
#   Windows     : .\relatorio\setup.ps1
# =======================================================================================

options(warn = 1, timeout = 600)

OK    <- "  [ok]   "
ACAO  <- "  [novo] "
FALHA <- "  [FALHA]"
AVISO <- "  [!]    "

falhas <- character(0)
passo <- function(n, txt) cat(sprintf("\n[%d/6] %s\n%s\n", n, txt, strrep("-", 72)))

cat("\n", strrep("=", 72), "\n",
    "PREPARANDO O AMBIENTE - RELATORIO PBL CAFEICULTURA\n",
    strrep("=", 72), "\n", sep = "")
cat("  sistema : ", R.version$platform, "\n", sep = "")
cat("  R       : ", R.version.string, "\n", sep = "")

# =======================================================================================
passo(1, "Versao do R")

if (getRversion() < "4.0.0") {
    cat(FALHA, "R ", as.character(getRversion()), " e antigo demais.\n", sep = "")
    cat("         Instale a versao 4.x em https://cran.r-project.org\n")
    falhas <- c(falhas, "versao do R")
} else {
    cat(OK, "R ", as.character(getRversion()), "\n", sep = "")
}

# =======================================================================================
passo(2, "Acesso a internet")

tem_net <- FALSE
for (u in c("https://cloud.r-project.org", "https://cran.rstudio.com")) {
    r <- try(suppressWarnings(readLines(u, n = 1, warn = FALSE)), silent = TRUE)
    if (!inherits(r, "try-error")) { tem_net <- TRUE; break }
}
if (tem_net) {
    cat(OK, "CRAN acessivel\n", sep = "")
} else {
    cat(AVISO, "sem acesso ao CRAN.\n", sep = "")
    cat("         Se tudo ja estiver instalado, o script segue e apenas confere.\n")
    cat("         Se faltar alguma coisa, vai falhar - conecte e rode de novo.\n")
}

# espelho do CRAN, caso nao esteja configurado
if (is.null(getOption("repos")[["CRAN"]]) || getOption("repos")[["CRAN"]] == "@CRAN@")
    options(repos = c(CRAN = "https://cloud.r-project.org"))

# =======================================================================================
passo(3, "Pacotes R")

# usados pelos modelos (numeros.R) e pelo build
PACOTES <- c("tinytex", "MASS", "betareg", "AER", "pROC", "emmeans", "lmtest",
             "sandwich", "car", "ggplot2", "dplyr", "patchwork")

biblioteca_gravavel <- function() {
    for (l in .libPaths()) if (file.access(l, 2) == 0) return(l)
    # nenhuma biblioteca do sistema e gravavel: cria uma no perfil do usuario
    pessoal <- Sys.getenv("R_LIBS_USER")
    if (!nzchar(pessoal))
        pessoal <- file.path(path.expand("~"), "R",
                             paste0(R.version$platform, "-library"),
                             paste(R.version$major, substr(R.version$minor, 1, 1), sep = "."))
    dir.create(pessoal, recursive = TRUE, showWarnings = FALSE)
    .libPaths(c(pessoal, .libPaths()))
    pessoal
}
LIB <- biblioteca_gravavel()
cat("  biblioteca: ", LIB, "\n\n", sep = "")

for (p in PACOTES) {
    if (requireNamespace(p, quietly = TRUE)) {
        cat(OK, p, "\n", sep = "")
        next
    }
    if (!tem_net) {
        cat(FALHA, p, " ausente e sem internet\n", sep = "")
        falhas <- c(falhas, paste("pacote", p))
        next
    }
    cat(ACAO, p, " - instalando...\n", sep = "")
    r <- try(install.packages(p, lib = LIB, quiet = TRUE), silent = TRUE)
    if (requireNamespace(p, quietly = TRUE)) {
        cat(OK, p, " instalado\n", sep = "")
    } else {
        cat(FALHA, p, " nao instalou\n", sep = "")
        falhas <- c(falhas, paste("pacote", p))
    }
}

# =======================================================================================
passo(4, "Distribuicao LaTeX")

motor_do_sistema <- NULL
for (m in c("latexmk", "pdflatex", "xelatex")) {
    if (nzchar(Sys.which(m))) { motor_do_sistema <- m; break }
}

tem_tinytex <- requireNamespace("tinytex", quietly = TRUE) &&
               isTRUE(try(tinytex::is_tinytex(), silent = TRUE))

if (tem_tinytex) {
    cat(OK, "TinyTeX ja instalado\n", sep = "")
} else if (!is.null(motor_do_sistema)) {
    cat(OK, "LaTeX do sistema encontrado: ", Sys.which(motor_do_sistema), "\n", sep = "")
    cat(AVISO, "sera usado no lugar do TinyTeX.\n", sep = "")
    cat("         Se faltar algum pacote .sty, instale pelo gerenciador dele\n")
    cat("         (MiKTeX Console, tlmgr) ou apague-o e rode este setup de novo.\n")
} else if (!requireNamespace("tinytex", quietly = TRUE)) {
    cat(FALHA, "pacote tinytex ausente - veja as falhas do passo 3\n", sep = "")
    falhas <- c(falhas, "tinytex")
} else {
    cat(ACAO, "nenhum LaTeX na maquina. Instalando o TinyTeX...\n", sep = "")
    cat("         São cerca de 250 MB. Pode levar alguns minutos.\n")
    cat("         Nao pede senha de administrador.\n\n")
    r <- try(tinytex::install_tinytex(force = FALSE), silent = TRUE)
    if (inherits(r, "try-error")) {
        cat(FALHA, "instalacao do TinyTeX falhou:\n", sep = "")
        cat("         ", conditionMessage(attr(r, "condition")), "\n")
        falhas <- c(falhas, "TinyTeX")
    } else {
        tem_tinytex <- isTRUE(try(tinytex::is_tinytex(), silent = TRUE))
        if (tem_tinytex) {
            cat("\n", OK, "TinyTeX instalado\n", sep = "")
        } else {
            cat("\n", AVISO, "instalou, mas o R ainda nao enxerga.\n", sep = "")
            cat("         Feche e reabra o terminal, e rode este setup de novo.\n")
            falhas <- c(falhas, "TinyTeX (reabrir o terminal)")
        }
    }
}

# =======================================================================================
passo(5, "Pacotes LaTeX exigidos pelo main.tex")

# Nomes de PACOTE do TeX Live, que nem sempre coincidem com o do \usepackage.
# abntex2cite vem dentro de abntex2; babel-portuges e o portugues do babel.
TEX <- c("biblatex", "biblatex-abnt", "biber", "logreq",
         "babel-portuges", "hyphen-portuguese", "lm", "amsmath", "amsfonts",
         "geometry", "setspace", "booktabs", "multirow", "caption", "float",
         "graphics", "hyperref", "url", "xcolor", "etoolbox", "koma-script",
         "oberdiek", "tools", "latex-bin", "csquotes")

if (tem_tinytex) {
    instalados <- try(tinytex::tl_pkgs(), silent = TRUE)
    if (inherits(instalados, "try-error")) instalados <- character(0)
    faltam <- setdiff(TEX, instalados)
    if (!length(faltam)) {
        cat(OK, "todos os ", length(TEX), " pacotes ja presentes\n", sep = "")
    } else {
        cat(ACAO, length(faltam), " pacote(s) a instalar: ",
            paste(head(faltam, 8), collapse = ", "),
            if (length(faltam) > 8) ", ..." else "", "\n", sep = "")
        r <- try(tinytex::tlmgr_install(faltam), silent = TRUE)
        if (inherits(r, "try-error"))
            cat(AVISO, "alguns nao instalaram. O TinyTeX tenta de novo durante o build.\n", sep = "")
        else
            cat(OK, "instalados\n", sep = "")
    }
} else if (!is.null(motor_do_sistema)) {
    cat(AVISO, "LaTeX do sistema: nao da para instalar pacotes automaticamente.\n", sep = "")
    cat("         Se a compilacao acusar arquivo .sty faltando, instale-o pelo\n")
    cat("         gerenciador da sua distribuicao LaTeX. O provavel e 'abntex2'.\n")
} else {
    cat(FALHA, "sem LaTeX - passo pulado\n", sep = "")
}

# =======================================================================================
passo(6, "Teste real de compilacao")

# So conferir se o binario existe nao prova nada. Compila um documento minimo
# que usa exatamente os recursos do relatorio: ABNT, portugues, tabela, cor.
teste_dir <- file.path(tempdir(), "teste_latex")
dir.create(teste_dir, showWarnings = FALSE, recursive = TRUE)
writeLines(c(
    "\\documentclass[12pt,a4paper]{article}",
    "\\usepackage[utf8]{inputenc}",
    "\\usepackage[T1]{fontenc}",
    "\\usepackage[brazilian]{babel}",
    "\\usepackage{lmodern}",
    "\\usepackage[top=3cm,bottom=2cm,left=3cm,right=2cm]{geometry}",
    "\\usepackage{indentfirst,setspace,amsmath,booktabs,graphicx,float,xcolor}",
    "\\usepackage[backend=biber,style=abnt]{biblatex}",
    "\\usepackage[hidelinks]{hyperref}",
    "\\onehalfspacing",
    "\\begin{document}",
    "Acentuacao: coordenacao, severidade, precisao, cafe, Itajai.",
    "\\textcolor{red}{Cor funciona.}",
    "\\begin{table}[H]\\centering\\begin{tabular}{lr}\\toprule",
    "Cultivar & Severidade \\\\ \\midrule Icatu & 0,0665 \\\\ \\bottomrule",
    "\\end{tabular}\\end{table}",
    "$\\mathrm{Var}[y] = \\frac{\\mu(1-\\mu)}{1+\\phi}$",
    "\\end{document}"
), file.path(teste_dir, "teste.tex"))

antigo <- setwd(teste_dir); on.exit(setwd(antigo), add = TRUE)
compilou <- FALSE

if (tem_tinytex) {
    r <- try(tinytex::latexmk("teste.tex", engine = "pdflatex"), silent = TRUE)
    compilou <- !inherits(r, "try-error") && file.exists("teste.pdf")
} else if (!is.null(motor_do_sistema)) {
    if (motor_do_sistema == "latexmk")
        system2("latexmk", c("-pdf", "-interaction=nonstopmode", "teste.tex"),
                stdout = FALSE, stderr = FALSE)
    else
        system2("pdflatex", c("-interaction=nonstopmode", "teste.tex"),
                stdout = FALSE, stderr = FALSE)
    compilou <- file.exists("teste.pdf")
}

if (compilou) {
    cat(OK, "documento de teste compilou - ABNT, portugues, tabela e cor\n", sep = "")
} else {
    cat(FALHA, "o documento de teste NAO compilou\n", sep = "")
    if (file.exists("teste.log")) {
        linhas <- readLines("teste.log", warn = FALSE)
        erros <- grep("^!", linhas, value = TRUE)
        if (length(erros)) {
            cat("\n         Erros do LaTeX:\n")
            for (e in head(erros, 5)) cat("           ", e, "\n")
        }
    }
    falhas <- c(falhas, "compilacao de teste")
}
setwd(antigo)

# =======================================================================================
cat("\n", strrep("=", 72), "\n", sep = "")
if (!length(falhas)) {
    cat("AMBIENTE PRONTO.\n\n")
    cat("Compile o relatorio com:\n")
    cat("   Rscript relatorio/build.R\n\n")
    cat("E confira antes de entregar com:\n")
    cat("   Rscript relatorio/conferir.R\n")
} else {
    cat("AMBIENTE INCOMPLETO. Pendencias:\n\n")
    for (f in unique(falhas)) cat("   - ", f, "\n", sep = "")
    cat("\nCorrija e rode este setup de novo - ele nao refaz o que ja deu certo.\n")
}
cat(strrep("=", 72), "\n\n")
quit(status = if (length(falhas)) 1 else 0)
