# =====================================================================
# ENTENDENDO O OFFSET - PARTE 1: OS FUNDAMENTOS
#
# Responde, em ordem, quatro duvidas que precisam vir ANTES do offset
# fazer sentido:
#   (c) o que e mu, e por que log(0) nao ser definido nao e problema
#   (b) o que e uma funcao de ligacao (link)
#   (a) de onde sai a equacao log(mu/E) = x'beta
#   (d) o que quer dizer "link canonico"
#
# Rode por blocos, nao de uma vez. Cada secao imprime a evidencia.
# =====================================================================

localizar_script <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) return(normalizePath(sub("^--file=", "", m[1]), winslash = "/", mustWork = FALSE))
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
    if (identical(pai, d)) stop("Abra a PASTA glm-cafeicultura, nao o arquivo solto.", call. = FALSE)
    d <- pai
  }
}
RAIZ <- raiz_repo()
caf  <- read.csv(file.path(RAIZ, "data", "cafeicultura.csv"))
caf$esforco <- caf$n_armadilhas * caf$dias_exposicao

sec <- function(t) cat("\n", strrep("=", 70), "\n", t, "\n", strrep("=", 70), "\n", sep = "")


# =====================================================================
sec("(c) O QUE E 'mu'  -  NAO e o numero de brocas")
# =====================================================================
# y  = o que foi OBSERVADO no talhao. Inteiro. Pode ser 0.
# mu = o valor ESPERADO para um talhao com aquelas caracteristicas.
#      E um parametro do modelo, nao um dado. Nao precisa ser inteiro.
#
# Mesma distincao de "esse dado deu 3" versus "esse dado tem esperanca 3.5".

cat("Observado (y): inteiro, e ha zeros\n")
print(table(caf$n_brocas_capturadas)[1:8])
cat("\nzeros observados:", sum(caf$n_brocas_capturadas == 0), "de", nrow(caf), "\n")

m0 <- glm(n_brocas_capturadas ~ umidade_relativa_pct, family = poisson, data = caf)
mu <- fitted(m0)
cat("\nEstimado (mu): continuo, e NUNCA zero\n")
cat("  faixa de mu :", round(min(mu), 4), "a", round(max(mu), 4), "\n")
cat("  mu inteiro? :", all(mu == round(mu)), "\n")
cat("  algum mu=0? :", any(mu == 0), "\n")

cat("\n>> A duvida do log(0): mu nunca chega a zero, entao log(mu) sempre\n")
cat("   existe. E y = 0 continua perfeitamente possivel - uma Poisson com\n")
cat("   mu pequeno produz zeros o tempo todo:\n")
set.seed(1)
for (m_ in c(0.5, 2, 9)) {
  am <- rpois(10000, m_)
  cat(sprintf("     mu = %4.1f  ->  %5.1f%% das observacoes saem 0\n", m_, 100 * mean(am == 0)))
}
cat("   mu > 0 e restricao sobre a MEDIA, nao sobre os dados.\n")


# =====================================================================
sec("(b) O QUE E UMA FUNCAO DE LIGACAO (link)")
# =====================================================================
# O problema que o link resolve:
#
#   eta = x'beta   e o "preditor linear". Pode dar QUALQUER real:
#                  -8, 0, +1000. E uma soma de termos; nada o limita.
#
#   mu             para contagem, tem de ser > 0.
#
# Nao da para escrever mu = eta: um beta negativo grande daria media
# negativa de brocas. A ligacao g() e a traducao entre as duas escalas:
#
#        g(mu) = eta        leva mu para a escala sem limites
#        mu = g^-1(eta)     traz eta de volta para a escala valida

cat("Os dois links deste problema:\n\n")
cat("  identidade:  g(mu) = mu        ->  mu = eta        (sem traducao)\n")
cat("  log       :  g(mu) = log(mu)   ->  mu = exp(eta)   (sempre > 0)\n\n")
cat("O que cada um faz com valores de eta:\n\n")
eta <- c(-5, -2, 0, 1, 3)
print(data.frame(
  eta           = eta,
  mu_identidade = eta,
  situacao_1    = ifelse(eta > 0, "ok", "<< MEDIA NEGATIVA"),
  mu_log        = round(exp(eta), 4),
  situacao_2    = "ok"), row.names = FALSE)

cat("\n>> Ponte com rede neural: o link INVERSO e a funcao de ativacao da\n")
cat("   camada de saida. Regressao logistica = 1 neuronio com sigmoide, e a\n")
cat("   sigmoide E o inverso do link logit. Aqui a 'ativacao' e exp(). A\n")
cat("   parte linear calcula eta; a ativacao poe no dominio valido.\n")


# =====================================================================
sec("(a) DE ONDE SAI  log(mu/E) = x'beta")
# =====================================================================
# Nao e deducao de nada. E a DECLARACAO do modelo, montada em duas
# decisoes nossas:
#
#   DECISAO 1 - o que explicar.
#     Nao a contagem mu (que depende de quantas armadilhas instalaram),
#     e sim a TAXA:            lambda = mu / E
#
#   DECISAO 2 - em que escala essa coisa e linear.
#     lambda tambem e > 0, entao vale o mesmo argumento do link:
#                              log(lambda) = x'beta
#
#   Substituindo 1 em 2:       log(mu / E) = x'beta
#
#   E dai em diante e so algebra - a parte que voce ja tinha entendido:
#                              log(mu) - log(E) = x'beta
#                              log(mu) = log(E) + x'beta
#
# Em portugues: "o log da taxa de captura e uma combinacao linear das
# caracteristicas do talhao". Nao e verdade universal; e a hipotese que
# estamos assumindo, e que a analise de residuos vai testar depois.

cat("As duas decisoes, e o que a alternativa custaria:\n\n")
cat("  DECISAO 1: modelar mu (contagem) -> comparar talhoes de esforco\n")
cat("             diferente vira comparar coisas diferentes (FATO 1).\n")
cat("             modelar mu/E (taxa)   -> comparavel.          <-- escolhida\n\n")
cat("  DECISAO 2: lambda = x'beta       -> taxa negativa possivel.\n")
cat("             log(lambda) = x'beta  -> taxa sempre positiva. <-- escolhida\n")


# =====================================================================
sec("(d) O QUE QUER DIZER 'LINK CANONICO'")
# =====================================================================
# Toda distribuicao da familia exponencial (Poisson, binomial, gama,
# normal) pode ser escrita numa forma padrao, com um "parametro natural"
# chamado theta. Para a Poisson essa forma e:
#
#     f(y; mu) = exp( y * log(mu) - mu - log(fatorial de y) )
#                         ^^^^^^^
#                      theta = log(mu)
#
# O parametro natural da Poisson E log(mu). A ligacao que leva mu
# exatamente ate o parametro natural chama-se CANONICA. Para a Poisson,
# portanto, a canonica e a log - nao por convencao, mas porque cai da
# propria forma da distribuicao.
#
# Na pratica: com o link canonico as equacoes de estimacao simplificam e
# o algoritmo de ajuste (IRLS) se comporta bem. Nao e obrigatorio usar.
# E so o caminho sem pedra.

f <- n_brocas_capturadas ~ umidade_relativa_pct + densidade_plantio + altitude_m
cat("O mesmo modelo, com os dois links:\n\n")

cat("  link log (canonico), sem nenhuma ajuda:\n")
m_log <- glm(f, family = poisson(link = "log"), data = caf)
cat("     converged:", m_log$converged, "| iteracoes:", m_log$iter, "\n\n")

cat("  link identidade, sem nenhuma ajuda:\n")
r <- tryCatch({ glm(f, family = poisson(link = "identity"), data = caf); "ajustou" },
              error = function(e) paste("ERRO:", conditionMessage(e)))
cat("    ", r, "\n\n")

cat("  link identidade, com valores iniciais dados na mao:\n")
m_id <- suppressWarnings(glm(f, family = poisson(link = "identity"), data = caf,
                             start = c(1, 0, 0, 0)))
cat("     converged:", m_id$converged, "| iteracoes:", m_id$iter,
    "(teto padrao =", m_id$control$maxit, ")\n\n")

cat("  link identidade, dando 500 iteracoes de folga:\n")
m_id2 <- suppressWarnings(glm(f, family = poisson(link = "identity"), data = caf,
                              start = c(1, 0, 0, 0),
                              control = glm.control(maxit = 500)))
cat("     converged:", m_id2$converged, "| iteracoes:", m_id2$iter, "\n")

cat("\n>> ATENCAO - aqui esta o FATO 3 de novo, em outra roupa. O modelo com\n")
cat("   link identidade NAO CONVERGE, nem com 500 iteracoes. E mesmo assim\n")
cat("   devolve um objeto com AIC (", round(AIC(m_id2), 1), ") e p-valores.\n", sep = "")
cat("   Quem olhar so o AIC compara ", round(AIC(m_id2), 1), " contra ",
    round(AIC(m_log), 1), " e conclui\n", sep = "")
cat("   que o log e melhor - conclusao certa, por um caminho invalido:\n")
cat("   nao se compara AIC de um modelo que nao convergiu.\n")

cat("\n>> 'Canonico' nao quer dizer 'obrigatorio'. Quer dizer que sai da\n")
cat("   algebra da distribuicao, e que o ajuste nao briga com voce.\n")

sec("FIM DA PARTE 1")
cat("A PARTE 2 (mecanica do offset) entra depois que estes quatro pontos\n")
cat("estiverem de pe.\n")
