# =====================================================================
# DESAFIO B - O QUI-QUADRADO E UMA DISTRIBUICAO, NAO UM TESTE
# =====================================================================
caf <- read.csv(file.path(
  if (file.exists("data/cafeicultura.csv")) "." else "..", "data/cafeicultura.csv"))
y <- caf$severidade_ferrugem
secao <- function(t) cat("\n", strrep("=", 64), "\n", t, "\n", strrep("=", 64), "\n", sep="")

secao("1. A CULPA DA CONFUSAO E DO R")
cat("Olhe o que ele imprime:\n\n")
print(kruskal.test(y ~ factor(caf$cultivar)))
cat("\nEle escreve 'Kruskal-Wallis chi-squared = 70.8566'. Isso junta duas\n")
cat("coisas diferentes na mesma linha, e e exatamente onde voce tropecou.\n")

secao("2. SAO DUAS COISAS SEPARADAS")
cat("  H = 70,86  ->  A ESTATISTICA. Calculada dos SEUS dados, pelos postos.\n")
cat("                 E o numero que mede o quanto os grupos se afastam.\n\n")
cat("  qui-quadrado -> A DISTRIBUICAO DE REFERENCIA. Nao tem nada dos seus\n")
cat("                 dados. E uma curva teorica, a REGUA usada para ler o H.\n\n")
cat("Analogia: H e a temperatura medida. O qui-quadrado e o termometro -\n")
cat("a escala que diz se 70,86 e frio ou quente.\n")

secao("3. POR QUE ESSA REGUA, E NAO OUTRA")
cat("Kruskal e Wallis PROVARAM que, se H0 for verdadeira e os grupos forem\n")
cat("grandes o bastante, o H se comporta como uma qui-quadrado com k-1\n")
cat("graus de liberdade, onde k = numero de grupos.\n\n")
cat("  k = 4 cultivares  ->  gl = 3\n\n")
cat("Como e essa curva - o que ela considera normal sob H0:\n\n")
q <- qchisq(c(.5,.9,.95,.99,.999), df=3)
for (i in seq_along(q))
  cat(sprintf("  %5.1f%% dos valores ficam abaixo de %6.2f\n",
              100*c(.5,.9,.95,.99,.999)[i], q[i]))
cat(sprintf("\n  seu H                            = %6.2f\n", 70.8566))
cat("  -> muito alem do 99,9%. Por isso o p e minusculo.\n\n")
cat("  p = area da cauda acima de 70,86 =",
    format(pchisq(70.8566, df=3, lower.tail=FALSE), digits=3), "\n")

secao("4. CONFERINDO A REGUA COM O EMBARALHAMENTO")
cat("O embaralhamento gera H sob H0 sem assumir curva nenhuma.\n")
cat("Se a teoria estiver certa, os H embaralhados devem seguir a qui-quadrado:\n\n")
set.seed(42)
acaso <- replicate(5000, unname(kruskal.test(y ~ factor(sample(caf$cultivar)))$statistic))
cmp <- data.frame(
  quantil = c("50%","90%","95%","99%"),
  embaralhamento = round(quantile(acaso, c(.5,.9,.95,.99)), 2),
  qui_quadrado_gl3 = round(qchisq(c(.5,.9,.95,.99), df=3), 2))
print(cmp, row.names = FALSE)
cat("\n-> batem bem. A aproximacao teorica cabe neste caso.\n")

secao("5. ONDE MAIS O QUI-QUADRADO APARECE NO SEU TRABALHO")
cat("Ele nao e 'o teste do Kruskal-Wallis'. E uma distribuicao que varios\n")
cat("testes diferentes usam como regua. No seu trabalho, tres lugares:\n\n")
cat("  (a) Kruskal-Wallis        - regua para o H dos postos\n")
cat("  (b) anova(m1, m2, test='Chisq') - regua para a diferenca de deviance\n")
cat("      entre dois modelos. E o que voces usaram no Desafio A.\n")
cat("  (c) teste qui-quadrado de Pearson - independencia em tabela de\n")
cat("      contingencia. Vai aparecer no Desafio C.\n\n")
cat("Sao perguntas diferentes, estatisticas diferentes, mesma regua.\n\n")
cat("Exemplo do (b), do Desafio A - repare no 'Pr(>Chi)':\n\n")
caf$esforco <- caf$n_armadilhas * caf$dias_exposicao
m0 <- glm(n_brocas_capturadas ~ 1, family=poisson, data=caf, offset=log(esforco))
m1 <- glm(n_brocas_capturadas ~ cultivar + umidade_relativa_pct,
          family=poisson, data=caf, offset=log(esforco))
print(anova(m0, m1, test="Chisq"))
cat("\nAqui a estatistica nao e H, e a DIFERENCA DE DEVIANCE. Outra conta,\n")
cat("mesma regua qui-quadrado.\n")
