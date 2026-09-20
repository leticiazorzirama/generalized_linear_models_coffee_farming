# =====================================================================
# DESAFIO B - O CRITERIO DO EMBARALHAMENTO E QUANTOS EXISTEM
# =====================================================================
caf <- read.csv(file.path(
  if (file.exists("data/cafeicultura.csv")) "." else "..", "data/cafeicultura.csv"))
y <- caf$severidade_ferrugem
secao <- function(t) cat("\n", strrep("=", 64), "\n", t, "\n", strrep("=", 64), "\n", sep="")

secao("1. O n E QUANTIDADE DE TALHOES - sim")
ni <- table(caf$cultivar)
print(data.frame(cultivar = names(ni), n_talhoes = as.vector(ni)), row.names = FALSE)
cat("\nsoma =", sum(ni), "= N, o total de talhoes da base\n")

secao("2. O CRITERIO DO EMBARALHAMENTO: nenhum")
cat("Esse e o ponto. NAO ha criterio - e permutacao uniformemente aleatoria.\n")
cat("Todo rearranjo possivel tem exatamente a mesma chance de sair.\n\n")
cat("O sample() do R faz isso: sorteia uma ordem sem repetir ninguem.\n")
cat("Se houvesse qualquer regra, o embaralhamento deixaria de representar\n")
cat("'o mundo sem efeito' - que e justamente o que queremos simular.\n\n")
cat("A unica coisa preservada e o TAMANHO de cada grupo: continuam 54, 106,\n")
cat("62 e 118 talhoes. Muda so QUAIS talhoes caem em cada um.\n")

secao("3. QUANTOS EMBARALHAMENTOS EXISTEM")
cat("Nao e 340! - permutacoes que trocam talhoes DENTRO do mesmo grupo dao\n")
cat("o mesmo agrupamento e o mesmo H. O que conta e o coeficiente multinomial:\n\n")
cat("   340! / (54! * 106! * 62! * 118!)\n\n")
log10n <- (lfactorial(340) - sum(lfactorial(as.vector(ni)))) / log(10)
cat("   = 10 elevado a", round(log10n, 1), "\n")
cat("   ou seja, um numero com", round(log10n), "digitos.\n\n")
cat("Para comparacao: atomos no universo observavel ~ 10^80.\n")

secao("4. SUA PERGUNTA: ALGUM DELES ALCANCA O H REAL?")
H_real <- unname(kruskal.test(y ~ factor(caf$cultivar))$statistic)
cat("SIM - e a resposta e mais simples do que parece:\n\n")
cat("  O ARRANJO REAL E, ELE PROPRIO, UM DOS EMBARALHAMENTOS POSSIVEIS.\n\n")
cat("Entao pelo menos 1 alcanca, sempre: o proprio. Logo o p exato NUNCA\n")
cat("e zero - o minimo possivel e 1 dividido pelo total.\n\n")

cat("E qual o H MAXIMO que algum arranjo consegue? Aquele que separa\n")
cat("perfeitamente: os 62 menores postos todos no Icatu, e assim por diante.\n\n")
postos <- rank(y); N <- length(y)
ordem_ideal <- rep(names(sort(tapply(rank(y), caf$cultivar, mean))), times =
                   as.vector(ni[names(sort(tapply(rank(y), caf$cultivar, mean)))]))
extremo <- y[order(y)]
H_max <- unname(kruskal.test(extremo ~ factor(ordem_ideal))$statistic)
cat("  H maximo possivel (separacao perfeita):", round(H_max, 2), "\n")
cat("  H dos dados reais                     :", round(H_real, 2), "\n")
cat("  H maximo que o acaso alcancou em 10.000:  21.21\n\n")
cat("O real esta a", round(100*H_real/H_max), "% do maximo teorico - longe de ser perfeito.\n")
cat("Mas esta a mais de 3x do que o acaso consegue tentando 10.000 vezes.\n")

secao("5. ENTAO O 'p = 0' QUE EU REPORTEI ESTAVA ERRADO")
cat("Eu escrevi 'p empirico = 0'. Isso e impreciso, pelo motivo acima.\n\n")
cat("Com B = 10.000 embaralhamentos e 0 acertos, a forma correta e:\n\n")
cat("  p_empirico = (acertos + 1) / (B + 1) = 1/10001 =",
    format(1/10001, digits = 3), "\n\n")
cat("Le-se: 'p < 0,0001'. O embaralhamento nao consegue medir abaixo disso -\n")
cat("faltaria rodar mais embaralhamentos. Para chegar em 2,8e-15 pela forca\n")
cat("bruta seriam necessarios ~10^15 embaralhamentos.\n\n")
cat("E por isso que existe a formula (o qui-quadrado): ela extrapola a cauda\n")
cat("sem precisar simular. O embaralhamento serve para ENTENDER o p e para\n")
cat("conferir a formula na regiao onde os dois se encontram.\n")
