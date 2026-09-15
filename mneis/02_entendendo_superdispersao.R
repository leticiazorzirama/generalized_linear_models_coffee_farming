# =====================================================================
# ENTENDENDO A SUPERDISPERSÃO
#
# Cada número que apareceu na conversa, refeito passo a passo e na mão,
# ao lado do que o R devolve. Responde:
#   (1) o que é "gl"
#   (2) como se calcula Pearson/gl, e por que dividir por gl
#   (3) que simulação foi aquela
#   (4) de onde saíram o "mu médio" e o "observado" da tabela por quartis
#
# Rode por blocos. Tudo aqui é conferível.
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
y <- caf$n_brocas_capturadas

sec <- function(t) cat("\n", strrep("=", 70), "\n", t, "\n", strrep("=", 70), "\n", sep = "")


# =====================================================================
sec("(1) O QUE E 'gl' - GRAUS DE LIBERDADE")
# =====================================================================
# gl = quantas observacoes voce tem MENOS quantos parametros voce gastou
#      estimando. E quanta informacao independente sobra depois do ajuste.
#
# Intuicao: com 2 pontos voce ajusta uma reta PERFEITA - sobra zero
# informacao para avaliar se a reta presta. Com 100 pontos e 2 parametros,
# sobram 98 pedacos de informacao para julgar o ajuste.

f <- n_brocas_capturadas ~ regiao_produtora + cultivar + manejo + irrigacao +
  altitude_m + declividade_pct + idade_lavoura_anos + densidade_plantio +
  adubacao_n_kg_ha + precipitacao_safra_mm + umidade_relativa_pct +
  ph_solo + materia_organica_pct

mod <- glm(f, family = poisson(link = "log"), data = caf, offset = log(esforco))

cat("observacoes (n)              :", nrow(caf), "\n")
cat("parametros estimados (posto) :", mod$rank, "\n")
cat("gl residual = n - posto      :", nrow(caf) - mod$rank, "\n")
cat("o que o R diz  df.residual() :", df.residual(mod), "\n")
cat("\nPor que 19 parametros e nao 13? As categoricas viram varias colunas:\n")
cat("  regiao(4 niveis)->3 colunas + cultivar(4)->3 + manejo(3)->2 +\n")
cat("  irrigacao(2)->1 + 9 numericas + 1 intercepto = 19\n")
cat("  (cada fator gasta niveis-1: o primeiro nivel fica embutido no intercepto)\n")


# =====================================================================
sec("(2) COMO SE CALCULA Pearson/gl - NA MAO")
# =====================================================================
# PASSO 1 - o que o modelo previu para cada talhao
mu <- fitted(mod)

# PASSO 2 - o residuo de Pearson de cada talhao.
#   r_i = (observado - previsto) / desvio-padrao que a Poisson AFIRMA ter
#   Para a Poisson, Var = mu, entao o desvio-padrao e sqrt(mu).
r_mao <- (y - mu) / sqrt(mu)
r_R   <- residuals(mod, type = "pearson")

cat("PASSO 1 e 2 - residuos de Pearson dos 5 primeiros talhoes:\n\n")
print(data.frame(
  talhao   = caf$id_talhao[1:5],
  y        = y[1:5],
  mu       = round(mu[1:5], 4),
  sd_pois  = round(sqrt(mu[1:5]), 4),
  r_na_mao = round(r_mao[1:5], 4),
  r_do_R   = round(r_R[1:5], 4)), row.names = FALSE)
cat("\nmax |diferenca| entre a conta na mao e a do R:",
    format(max(abs(r_mao - r_R)), digits = 3), "\n")

cat("\nLeitura de um residuo: TAL001 observou", y[1], "broca onde o modelo\n")
cat("esperava", round(mu[1], 2), ". Isso e", round(r_mao[1], 2),
    "desvios-padrao abaixo - pela regua da Poisson.\n")

# PASSO 3 - somar os quadrados: a estatistica de Pearson
X2 <- sum(r_mao^2)
cat("\nPASSO 3 - estatistica de Pearson  X2 = soma dos r_i^2 =", round(X2, 3), "\n")

# PASSO 4 - dividir pelos graus de liberdade
cat("PASSO 4 - dispersao = X2 / gl =", round(X2, 3), "/", df.residual(mod),
    "=", round(X2 / df.residual(mod), 4), "\n")

cat("\nPOR QUE DIVIDIR POR gl?\n")
cat("Se a Poisson estiver certa, cada r_i tem variancia ~1, entao cada r_i^2\n")
cat("contribui ~1 para a soma. Com gl parcelas, a soma daria ~gl, e a razao ~1.\n")
cat("Ou seja: a razao e o QUADRADO MEDIO dos residuos padronizados.\n")
cat("  razao ~ 1 -> a Poisson descreve bem o espalhamento\n")
cat("  razao > 1 -> SUPERDISPERSAO: espalha mais do que a Poisson permite\n")
cat("  razao < 1 -> subdispersao\n")
cat("\nAqui deu", round(X2/df.residual(mod), 2),
    "-> os residuos espalham", round(X2/df.residual(mod), 1),
    "vezes mais do que a Poisson admite.\n")


# =====================================================================
sec("(2b) AS TRES LINHAS DA TABELA - cada uma e um modelo diferente")
# =====================================================================
razao <- function(m) sum(residuals(m, type = "pearson")^2) / df.residual(m)

m_a <- glm(n_brocas_capturadas ~ 1, family = poisson, data = caf)
m_b <- glm(n_brocas_capturadas ~ 1, family = poisson, data = caf, offset = log(esforco))
m_c <- mod

print(data.frame(
  modelo     = c("nulo, sem offset", "nulo, com offset", "completo, com offset"),
  formula    = c("y ~ 1", "y ~ 1 + offset", "y ~ 13 variaveis + offset"),
  parametros = c(m_a$rank, m_b$rank, m_c$rank),
  gl         = c(df.residual(m_a), df.residual(m_b), df.residual(m_c)),
  dispersao  = round(c(razao(m_a), razao(m_b), razao(m_c)), 3)), row.names = FALSE)

cat("\nA leitura da queda:\n")
cat("  7,60 -> nada explicado. Todo o espalhamento conta como dispersao.\n")
cat("  5,32 -> o offset explicou a parte devida ao esforco amostral.\n")
cat("  2,50 -> as covariaveis explicaram mais um tanto.\n")
cat("O que sobra em 2,50 e o que nenhuma variavel explicou: a superdispersao real.\n")

cat("\nConfere que a primeira linha e mesmo var/media dos dados crus:\n")
cat("  var(y)/mean(y) =", round(var(y)/mean(y), 3),
    " | Pearson/gl do modelo nulo =", round(razao(m_a), 3), "\n")
cat("  (batem porque no modelo nulo mu = media para todos, e gl = n-1)\n")


# =====================================================================
sec("(3) QUE SIMULACAO FOI AQUELA")
# =====================================================================
# rpois(n, mu) sorteia n numeros de uma Poisson de media mu.
# A ideia: pedir uma media, sortear muitas vezes, e MEDIR a variancia que
# saiu. Se a Poisson impoe Var = media, a variancia medida tem de bater
# com a media pedida - sem que ninguem tenha pedido isso.

set.seed(42)
cat("Exemplo pequeno primeiro - 10 sorteios de uma Poisson de media 9:\n")
cat(" ", rpois(10, 9), "\n\n")

cat("Agora 200.000 sorteios por linha, e medindo o que saiu:\n\n")
sim <- data.frame(mu_pedido = c(0.5, 2, 9, 20))
sim$media_medida <- sapply(sim$mu_pedido, function(m) round(mean(rpois(2e5, m)), 3))
sim$var_medida   <- sapply(sim$mu_pedido, function(m) round(var(rpois(2e5, m)), 3))
sim$diferenca    <- round(sim$var_medida - sim$mu_pedido, 3)
print(sim, row.names = FALSE)

cat("\nAs colunas 2 e 3 colam na 1 sempre. Nao e coincidencia nem escolha:\n")
cat("a Poisson tem UM parametro, e ele governa media e variancia juntos.\n")
cat("Nao existe argumento em rpois() para pedir media 9 com variancia 60.\n")


# =====================================================================
sec("(4) A TABELA POR QUARTIS - de onde saem as colunas")
# =====================================================================
# COLUNA 'mu medio'  : media dos mu PREVISTOS dos talhoes daquele grupo
# COLUNA 'observado' : variancia dos y REAIS daqueles mesmos talhoes
# A ideia: agrupar talhoes de infestacao parecida e comparar o espalhamento
# real contra o que a Poisson exigiria naquele nivel.

cortes <- quantile(mu, seq(0, 1, 0.25))
grupo  <- cut(mu, cortes, include.lowest = TRUE)

cat("PASSO 1 - os quartis de mu (os pontos de corte):\n")
print(round(cortes, 3))

cat("\nPASSO 2 - para cada grupo, na mao:\n\n")
res <- do.call(rbind, lapply(split(seq_along(y), grupo), function(ix)
  data.frame(
    n                 = length(ix),
    mu_medio          = round(mean(mu[ix]), 2),   # media dos PREVISTOS
    poisson_exigiria  = round(mean(mu[ix]), 2),   # ... que e a mesma coisa: Var = mu
    var_observada     = round(var(y[ix]), 2),     # variancia dos OBSERVADOS
    razao             = round(var(y[ix]) / mean(mu[ix]), 2))))
print(res)

cat("\nConferindo o primeiro grupo na unha:\n")
ix <- which(grupo == levels(grupo)[1])
cat("  talhoes no grupo     :", length(ix), "\n")
cat("  mean(mu[grupo])      :", round(mean(mu[ix]), 4), "<- coluna 'mu medio'\n")
cat("  var(y[grupo])        :", round(var(y[ix]), 4), "<- coluna 'observado'\n")
cat("  e var(y) na mao      :",
    round(sum((y[ix] - mean(y[ix]))^2) / (length(ix) - 1), 4), "\n")

cat("\nPOR QUE 'poisson exigiria' = 'mu medio'? Porque a Poisson diz Var = mu.\n")
cat("Entao a variancia exigida naquele nivel E o proprio mu. Mesma coluna,\n")
cat("nomes diferentes, para deixar a comparacao explicita.\n")

cat("\n>> RESSALVA: dentro de um grupo os mu nao sao identicos, so parecidos.\n")
cat("   Isso infla um pouco var(y) e a razao. E diagnostico visual, nao teste\n")
cat("   formal - o teste formal (AER::dispersiontest) vem na Tarefa 2.\n")

cat("\n>> O QUE IMPORTA NESSA TABELA: a razao CRESCE com mu (1,5 -> 6,1).\n")
cat("   Nao e um fator constante. Guardar isso: quasi-Poisson e binomial\n")
cat("   negativa preveem coisas diferentes sobre esse padrao, e e por ai que\n")
cat("   a Tarefa 3 se decide.\n")

sec("FIM")
