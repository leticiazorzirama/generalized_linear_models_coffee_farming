# =====================================================================
# DESAFIO B - POSTOS, POSTO MEDIO E EMBARALHAMENTO
# Exemplo minusculo: 8 talhoes, 2 cultivares. Da para conferir a olho.
# =====================================================================
secao <- function(t) cat("\n", strrep("=", 64), "\n", t, "\n", strrep("=", 64), "\n", sep="")

secao("0. VOCABULARIO - cultivar NAO e rotulo")
cat("CULTIVAR  = variedade da planta de cafe. Termo de agronomia, nao de\n")
cat("            estatistica. Catuai, Bourbon, Icatu e Mundo Novo sao quatro\n")
cat("            variedades de Coffea arabica - como racas de cachorro.\n")
cat("            Na base, e o NOME DA COLUNA.\n\n")
cat("ROTULO    = palavra generica que eu usei para 'o valor daquela coluna\n")
cat("            naquela linha'. Foi confuso. Dali em diante: valor da coluna.\n\n")
cat("  coluna  -> cultivar\n")
cat("  valor   -> 'Catuai'\n")

# ---------------------------------------------------------------------
mini <- data.frame(
  talhao     = c("T1","T2","T3","T4","T5","T6","T7","T8"),
  cultivar   = c("Icatu","Bourbon","Icatu","Bourbon","Icatu","Bourbon","Icatu","Bourbon"),
  severidade = c(0.05, 0.22, 0.03, 0.18, 0.09, 0.31, 0.06, 0.25))

secao("1. O QUE E POSTO")
cat("Posto = a POSICAO do valor depois de ordenar tudo, do menor ao maior.\n")
cat("O menor recebe posto 1, o maior recebe posto 8.\n\n")
ord <- mini[order(mini$severidade), ]
ord$posto <- 1:8
print(ord[, c("talhao","cultivar","severidade","posto")], row.names = FALSE)
cat("\nRepare: o posto ignora A DISTANCIA. 0,03 e 0,05 viram 1 e 2;\n")
cat("0,25 e 0,31 viram 7 e 8. So a ordem conta.\n")

mini$posto <- rank(mini$severidade)

secao("2. O QUE E POSTO MEDIO, e como se calcula")
cat("Posto medio de uma cultivar = media dos postos dos talhoes dela.\n\n")
for (cv in c("Icatu","Bourbon")) {
  p <- mini$posto[mini$cultivar == cv]
  cat(" ", cv, ": postos", paste(p, collapse=" + "), "=", sum(p),
      " ->  /", length(p), " =  posto medio", mean(p), "\n")
}
cat("\nposto medio GERAL = (N+1)/2 = (8+1)/2 = 4.5\n")
cat("  (e sempre isso: a media de 1,2,...,N)\n\n")
cat("Icatu  ficou em 2.5 -> 2 postos ABAIXO do geral\n")
cat("Bourbon ficou em 6.5 -> 2 postos ACIMA do geral\n")
cat("Se a cultivar nao importasse, os dois estariam perto de 4.5.\n")

secao("3. OS SIMBOLOS DA FORMULA")
cat("  H   = a propria estatistica. E so o nome dela (Kruskal-Wallis).\n")
cat("  N   = total de observacoes. Aqui 8. Na base real, 340.\n")
cat("  n_i = quantas observacoes tem o grupo i. Aqui 4 e 4.\n")
cat("  R_i = posto medio do grupo i (a barra em cima significa 'media').\n\n")
N <- 8; ni <- c(4,4); Rm <- c(2.5, 6.5); Rgeral <- (N+1)/2
cat("  H = 12/(N(N+1)) * soma[ n_i * (R_i - R_geral)^2 ]\n\n")
cat("    = 12/(8*9)  * [ 4*(2.5-4.5)^2  +  4*(6.5-4.5)^2 ]\n")
cat("    = 12/72     * [ 4*4            +  4*4           ]\n")
cat("    = ", round(12/72,4), "     * ", 4*4+4*4, "\n")
H <- 12/(N*(N+1)) * sum(ni*(Rm-Rgeral)^2)
cat("    = ", H, "\n\n")
cat("  o que o R devolve:",
    round(unname(kruskal.test(severidade ~ factor(cultivar), data=mini)$statistic), 4), "\n")
cat("\n  Leitura: H grande = os grupos estao LONGE do posto medio geral.\n")
cat("           H = 0     = os grupos estao exatamente no geral.\n")

secao("4. COMO O EMBARALHAMENTO ACONTECE")
cat("A coluna severidade fica PARADA. So a coluna cultivar e redistribuida.\n")
cat("E sim: a cada embaralhamento, recalcula-se H do zero.\n\n")
set.seed(3)
cat("            severidade:", sprintf("%5.2f", mini$severidade), "\n")
cat("            real       :", sprintf("%7s", substr(mini$cultivar,1,5)), " -> H =",
    round(H,3), "\n")
for (k in 1:3) {
  emb <- sample(mini$cultivar)
  Hk <- unname(kruskal.test(mini$severidade ~ factor(emb))$statistic)
  cat("  embaralho", k, ":", sprintf("%7s", substr(emb,1,5)), " -> H =", round(Hk,3), "\n")
}
cat("\nFaz-se isso 10.000 vezes e junta-se os 10.000 H num monte. O p e a\n")
cat("fracao desse monte que alcancou o H real.\n")

secao("5. SUA PERGUNTA: SE NAO HOUVESSE EFEITO, p SERIA 1?")
cat("Nao. Essa e a parte contra-intuitiva.\n\n")
cat("Mesmo sem efeito nenhum, o acaso da diferencas pequenas entre grupos.\n")
cat("Entao H quase nunca da 0, e p quase nunca da 1.\n\n")
cat("Simulando 5.000 experimentos onde os grupos sao REALMENTE iguais\n")
cat("(mesma distribuicao, so rotulos sorteados):\n\n")
set.seed(9)
ps <- replicate(5000, {
  v <- rnorm(120); g <- rep(c("A","B","C"), each = 40)
  kruskal.test(v ~ factor(g))$p.value })
faixas <- table(cut(ps, seq(0,1,0.1)))
for (i in seq_along(faixas))
  cat(sprintf("  p em %-10s %4d  %s\n", names(faixas)[i], faixas[i],
              strrep("#", round(faixas[i]/12))))
cat("\n  media dos p:", round(mean(ps),3), " | mediana:", round(median(ps),3), "\n\n")
cat("O p sai UNIFORME entre 0 e 1 - qualquer valor e igualmente provavel.\n")
cat("Nao se acumula em 1. Em media da 0,5.\n\n")
cat("p = 1 exigiria os grupos IDENTICOS, H exatamente 0. Isso praticamente\n")
cat("nao acontece nem quando nao ha efeito algum.\n\n")
cat("E repare: ", sum(ps < 0.05), " dos 5.000 deram p < 0,05 (", 
    round(100*mean(ps < 0.05),1), "%) mesmo SEM efeito nenhum.\n", sep="")
cat("E o falso positivo que o corte de 5% aceita por construcao.\n")
