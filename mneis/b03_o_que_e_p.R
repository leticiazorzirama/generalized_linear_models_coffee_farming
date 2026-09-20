# =====================================================================
# DESAFIO B - O QUE E O p-VALOR
# Construido por embaralhamento, para nao ficar abstrato.
# =====================================================================
caf <- read.csv(file.path(
  if (file.exists("data/cafeicultura.csv")) "." else "..", "data/cafeicultura.csv"))
y <- caf$severidade_ferrugem
secao <- function(t) cat("\n", strrep("=", 66), "\n", t, "\n", strrep("=", 66), "\n", sep="")

secao("1. A PERGUNTA QUE O p RESPONDE")
cat("Observamos: Icatu 0,0665 e Bourbon 0,2122. Diferenca de 3,19x.\n\n")
cat("Mas uma diferenca dessas poderia aparecer POR ACASO, so porque\n")
cat("dividimos 340 talhoes em 4 grupos? O p responde exatamente isso:\n\n")
cat("  'SE nao houvesse diferenca nenhuma entre as cultivares,\n")
cat("   qual a chance de eu ver uma diferenca DESTE tamanho ou maior?'\n")

secao("2. VAMOS MEDIR ISSO EMBARALHANDO")
cat("Se a cultivar nao importasse, o rotulo de cada talhao seria irrelevante.\n")
cat("Entao: embaralho os rotulos 10.000 vezes e vejo que diferencas aparecem\n")
cat("por puro acaso.\n\n")

estat <- function(rotulos) unname(kruskal.test(y ~ factor(rotulos))$statistic)
real  <- estat(caf$cultivar)
set.seed(42)
acaso <- replicate(10000, estat(sample(caf$cultivar)))

cat("estatistica de Kruskal-Wallis nos dados REAIS :", round(real, 2), "\n\n")
cat("nos 10.000 embaralhamentos:\n")
cat("  minimo :", round(min(acaso), 2), "\n")
cat("  mediana:", round(median(acaso), 2), "\n")
cat("  maximo :", round(max(acaso), 2), "  <- o maior que o acaso conseguiu\n\n")
cat("quantos embaralhamentos chegaram a", round(real, 2), "?  ", sum(acaso >= real), "de 10.000\n")
cat("p empirico =", sum(acaso >= real) / 10000, "\n\n")
cat("Ou seja: embaralhei 10 mil vezes e o acaso NUNCA chegou perto.\n")
cat("O p do Kruskal-Wallis (2,8e-15) diz a mesma coisa por formula:\n")
cat("a chance disso ser acaso e de 3 em 1.000.000.000.000.000.\n")

secao("3. COMO LER")
cat("  p pequeno  -> o que vi seria MUITO improvavel se nao houvesse efeito\n")
cat("                -> ha evidencia de efeito\n")
cat("  p grande   -> o que vi cabe tranquilamente no acaso\n")
cat("                -> nao ha evidencia. NAO significa 'provei que nao ha efeito'\n\n")
cat("O corte 0,05 e convencao, nao lei. 0,049 e 0,051 nao sao mundos diferentes.\n\n")
cat("E o p NAO e a probabilidade de a hipotese ser verdadeira. E a probabilidade\n")
cat("dos DADOS, assumindo que a hipotese e verdadeira. Ordem invertida.\n")

secao("4. A ARMADILHA: p NAO MEDE TAMANHO DE EFEITO")
cat("Olhe sua propria tabela do passo 2:\n\n")
cmp <- data.frame(
  fator = c("cultivar","regiao_produtora","irrigacao","manejo"),
  razao_entre_extremos = c(3.19, 1.40, 1.25, 1.11),
  p = c("2.8e-15","0.025","0.016","0.889"))
print(cmp, row.names = FALSE)
cat("\nIRRIGACAO tem efeito MENOR que regiao (1,25x contra 1,40x),\n")
cat("mas p MENOR tambem (0,016 contra 0,025). Como?\n\n")
cat("Porque o p mistura tres coisas: tamanho do efeito, tamanho da amostra\n")
cat("e o desenho (irrigacao tem 2 niveis; regiao tem 4, e reparte a amostra).\n\n")
cat("Demonstrando - mesmo efeito, n diferente:\n\n")
set.seed(1)
for (n in c(20, 60, 200, 1000)) {
  a <- rnorm(n, 0.10, 0.05); b <- rnorm(n, 0.12, 0.05)   # efeito FIXO de 0,02
  cat(sprintf("  n = %4d por grupo -> diferenca %.4f, p = %s\n",
              n, mean(b) - mean(a),
              format.pval(wilcox.test(a, b)$p.value, digits = 2)))
}
cat("\nO efeito real e sempre 0,02. So o n muda, e o p despenca.\n")
cat("-> p responde SE ha efeito. Nunca QUANTO. Para o quanto, olhe a razao,\n")
cat("   o coeficiente, o intervalo de confianca.\n")
