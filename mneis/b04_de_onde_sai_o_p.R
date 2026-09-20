# =====================================================================
# DESAFIO B - DE ONDE SAI O p: o nome, a estatistica, o embaralhamento
# =====================================================================
caf <- read.csv(file.path(
  if (file.exists("data/cafeicultura.csv")) "." else "..", "data/cafeicultura.csv"))
y <- caf$severidade_ferrugem
secao <- function(t) cat("\n", strrep("=", 66), "\n", t, "\n", strrep("=", 66), "\n", sep="")

secao("1. O QUE E 'ROTULO'")
cat("Rotulo = o valor da variavel categorica colado naquele talhao.\n")
cat("O TAL001 carrega o rotulo 'Catuai'. Sao duas colunas separadas:\n\n")
print(head(caf[, c("id_talhao","cultivar","severidade_ferrugem")], 6))
cat("\nEmbaralhar = manter a coluna de severidade PARADA e redistribuir\n")
cat("os 340 rotulos ao acaso entre os talhoes:\n\n")
set.seed(7)
demo <- head(data.frame(
  id = caf$id_talhao, severidade = caf$severidade_ferrugem,
  rotulo_real = caf$cultivar, rotulo_embaralhado = sample(caf$cultivar)), 6)
print(demo, row.names = FALSE)
cat("\nNAO sorteei um embaralhamento. Fiz 10.000 embaralhamentos COMPLETOS,\n")
cat("e de cada um extrai uma estatistica. Sao 10.000 numeros, nao 1.\n")

secao("2. DE ONDE SAI O 70,86 - a estatistica de Kruskal-Wallis, na mao")
cat("Ela trabalha com POSTOS, nao com os valores - igual ao Spearman.\n\n")
postos <- rank(y)                       # 1 = menor severidade ... 340 = maior
N <- length(y)
cat("PASSO 1 - ordenar os 340 talhoes por severidade e numerar de 1 a", N, "\n")
cat("  posto do menor:", min(postos), " | posto do maior:", max(postos), "\n")
cat("  posto medio geral = (N+1)/2 =", (N+1)/2, "\n\n")

cat("PASSO 2 - o posto medio DENTRO de cada cultivar:\n")
pm <- tapply(postos, caf$cultivar, mean)
ni <- table(caf$cultivar)
print(round(data.frame(n = as.vector(ni), posto_medio = round(as.vector(pm),1),
                       desvio_do_geral = round(as.vector(pm) - (N+1)/2, 1),
                       row.names = names(pm)), 1))
cat("\n  Se a cultivar nao importasse, as 4 linhas ficariam perto de", (N+1)/2, ".\n")

cat("\nPASSO 3 - somar os desvios ao quadrado, ponderados pelo n do grupo:\n")
H <- 12/(N*(N+1)) * sum(ni * (pm - (N+1)/2)^2)
cat("  H = 12/(N(N+1)) * soma[ n_i * (posto_medio_i - posto_medio_geral)^2 ]\n")
cat("  H =", round(H, 4), "\n")

cat("\nPASSO 4 - correcao para empates (ha severidades repetidas):\n")
t <- table(y); emp <- sum(t^3 - t)
Hc <- H / (1 - emp/(N^3 - N))
cat("  valores repetidos:", sum(t > 1), "| fator de correcao:",
    round(1 - emp/(N^3-N), 6), "\n")
cat("  H corrigido =", round(Hc, 4), "\n")
cat("  o que o R devolve =", round(unname(kruskal.test(y ~ factor(caf$cultivar))$statistic), 4), "\n")

secao("3. O NOME E A ORIGEM DO p")
cat("Chama-se p-valor. O 'p' e de PROBABILITY - e probabilidade mesmo,\n")
cat("um numero entre 0 e 1.\n\n")
cat("Mas probabilidade DE QUE, exatamente:\n\n")
cat("   p = P( estatistica tao ou mais extrema | H0 verdadeira )\n\n")
cat("H0 e a 'hipotese nula' - a afirmacao que se assume enquanto se calcula.\n")
cat("Aqui H0 = 'as 4 cultivares tem a mesma distribuicao de severidade'.\n\n")
cat("A barra | le-se 'dado que'. Entao o p e a probabilidade DOS DADOS\n")
cat("assumindo H0 - e nao a probabilidade de H0 ser verdadeira. Inverter\n")
cat("isso e o erro mais comum com p-valores.\n\n")
cat("Origem: Ronald Fisher, anos 1920. O corte 0,05 tambem e dele, e ele\n")
cat("mesmo tratava como referencia pratica, nao como lei.\n")

secao("4. OS DOIS CAMINHOS ATE O MESMO p")
cat("CAMINHO A - formula: assume-se que H segue uma distribuicao qui-quadrado\n")
cat("  com 3 graus de liberdade sob H0, e pergunta-se a area da cauda:\n")
cat("  p =", format(pchisq(Hc, df = 3, lower.tail = FALSE), digits = 3), "\n\n")
cat("CAMINHO B - embaralhamento: nao assume distribuicao nenhuma, so conta.\n")
set.seed(42)
acaso <- replicate(10000, unname(kruskal.test(y ~ factor(sample(caf$cultivar)))$statistic))
cat("  embaralhamentos que alcancaram", round(Hc,2), ":", sum(acaso >= Hc), "de 10.000\n")
cat("  p empirico =", sum(acaso >= Hc)/10000, "\n\n")
cat("Os dois respondem a mesma pergunta. O primeiro usa teoria, o segundo\n")
cat("usa forca bruta. Concordarem e sinal de que a teoria cabe aqui.\n")
