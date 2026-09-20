# =====================================================================
# DESAFIO B - ITEM 0.5: MULTICOLINEARIDADE
# Exigencia literal do §4b do enunciado.
# =====================================================================
if (!require(pacman)) { install.packages("pacman"); library(pacman) }
p_load(car, betareg, ggplot2)

caf <- read.csv(file.path(
  if (file.exists("data/cafeicultura.csv")) "." else "..", "data/cafeicultura.csv"),
  stringsAsFactors = TRUE)
secao <- function(t) cat("\n", strrep("=", 68), "\n", t, "\n", strrep("=", 68), "\n", sep="")

secao("1. O QUE E")
cat("Multicolinearidade = quando as COVARIAVEIS sao fortemente correlacionadas\n")
cat("ENTRE SI. Nao com a resposta - entre elas.\n\n")
cat("Exemplo caricato: se a base tivesse 'altitude_m' e 'altitude_pes', as duas\n")
cat("diriam a mesma coisa. O modelo nao teria como decidir qual merece o credito,\n")
cat("e distribuiria o efeito entre elas de forma arbitraria.\n")

secao("2. POR QUE E PROBLEMA - e por que NAO e")
cat("NAO quebra a predicao. O modelo continua acertando y.\n")
cat("QUEBRA a interpretacao de cada coeficiente isolado:\n\n")
cat("  - erros-padrao inflam -> intervalos largos, p-valores grandes\n")
cat("  - o sinal do coeficiente pode inverter\n")
cat("  - pequenas mudancas nos dados mudam muito as estimativas\n\n")
cat("Como o Desafio B e sobre INFERENCIA ('existe uma cultivar mais resistente?'),\n")
cat("e a interpretacao que importa. Por isso o enunciado cobra a verificacao.\n")

secao("3. DEMONSTRANDO O ESTRAGO")
cat("Vou criar uma copia quase identica da umidade e por as duas no modelo:\n\n")
set.seed(1)
caf$umidade_copia <- caf$umidade_relativa_pct + rnorm(nrow(caf), 0, 0.5)
cat("correlacao entre umidade e a copia:",
    round(cor(caf$umidade_relativa_pct, caf$umidade_copia), 4), "\n\n")

m_ok  <- lm(severidade_ferrugem ~ umidade_relativa_pct + densidade_plantio, data = caf)
m_bad <- lm(severidade_ferrugem ~ umidade_relativa_pct + umidade_copia + densidade_plantio,
            data = caf)
cmp <- rbind(
  data.frame(modelo="sem a copia", termo="umidade_relativa_pct",
             beta=coef(m_ok)[2], ep=summary(m_ok)$coefficients[2,2],
             p=summary(m_ok)$coefficients[2,4]),
  data.frame(modelo="COM a copia", termo="umidade_relativa_pct",
             beta=coef(m_bad)[2], ep=summary(m_bad)$coefficients[2,2],
             p=summary(m_bad)$coefficients[2,4]),
  data.frame(modelo="COM a copia", termo="umidade_copia",
             beta=coef(m_bad)[3], ep=summary(m_bad)$coefficients[3,2],
             p=summary(m_bad)$coefficients[3,4]))
cmp$beta <- round(cmp$beta, 5); cmp$ep <- round(cmp$ep, 5)
cmp$p <- format.pval(cmp$p, digits=2)
print(cmp, row.names=FALSE)
cat("\nO erro-padrao explodiu de", round(summary(m_ok)$coefficients[2,2],5), "para",
    round(summary(m_bad)$coefficients[2,2],5), "-",
    round(summary(m_bad)$coefficients[2,2]/summary(m_ok)$coefficients[2,2],1), "vezes maior.\n")
cat("A variavel que era altamente significativa virou nao significativa.\n")
cat("E a predicao? R2 sem a copia:", round(summary(m_ok)$r.squared,5),
    "| com a copia:", round(summary(m_bad)$r.squared,5), "- praticamente igual.\n")

secao("4. COMO SE MEDE: o VIF")
cat("VIF = Variance Inflation Factor (fator de inflacao da variancia).\n\n")
cat("Para cada covariavel j:\n")
cat("  1. ajusta-se uma regressao de j contra TODAS as outras covariaveis\n")
cat("  2. pega-se o R2 dessa regressao - quanto as outras explicam de j\n")
cat("  3. VIF_j = 1 / (1 - R2_j)\n\n")
cat("Leitura: VIF = quantas vezes a VARIANCIA do coeficiente foi inflada.\n")
cat("        sqrt(VIF) = quantas vezes o ERRO-PADRAO foi inflado (mais intuitivo).\n\n")
cat("  VIF = 1   -> covariavel independente das outras. Sem inflacao.\n")
cat("  VIF > 5   -> merece atencao\n")
cat("  VIF > 10  -> problema serio, convencao mais usada\n\n")
cat("Conferindo no exemplo estragado:\n")
print(round(vif(m_bad), 2))
cat("\n1/(1-R2) na mao para a umidade:",
    round(1/(1 - summary(lm(umidade_relativa_pct ~ umidade_copia + densidade_plantio,
                            data=caf))$r.squared), 2), "\n")

secao("5. O VIF DA NOSSA BASE")
cat("O VIF depende so das COVARIAVEIS (a matriz X), nao da resposta nem da\n")
cat("familia. Entao calcula-se num lm com as mesmas covariaveis - e valido\n")
cat("para a regressao beta tambem.\n\n")
f <- severidade_ferrugem ~ cultivar + manejo + irrigacao + regiao_produtora +
  umidade_relativa_pct + densidade_plantio + precipitacao_safra_mm + altitude_m +
  idade_lavoura_anos + adubacao_n_kg_ha + ph_solo + declividade_pct + materia_organica_pct
v <- vif(lm(f, data = caf))
tb <- data.frame(termo=rownames(v), GVIF=round(v[,1],3), gl=v[,2],
                 GVIF_ajustado=round(v[,3],3))
tb$equivale_a_VIF <- round(tb$GVIF_ajustado^2, 2)
tb <- tb[order(-tb$equivale_a_VIF), ]
print(tb, row.names=FALSE)

cat("\nNOTA sobre fatores: para variaveis categoricas o car devolve GVIF\n")
cat("(generalizado), porque elas ocupam varias colunas. A coluna comparavel\n")
cat("com o VIF classico e GVIF^(1/(2*gl)) elevado ao quadrado - a ultima acima.\n")

cat("\n# Parecer\n")
pior <- tb$equivale_a_VIF[1]
cat("O maior valor e", pior, "- muito abaixo de 5, quanto mais de 10.\n")
cat("Nao ha multicolinearidade preocupante nesta base. Nenhuma covariavel\n")
cat("precisa ser removida por esse motivo, e os coeficientes do modelo podem\n")
cat("ser interpretados individualmente sem ressalva de instabilidade.\n")


# =====================================================================
# GRAFICOS
# =====================================================================
p_load(ggplot2, patchwork, dplyr)
RAIZ <- if (file.exists("data/cafeicultura.csv")) "." else ".."
FIG  <- file.path(RAIZ, "mneis", "figuras")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
VERDE <- "#18675A"; VERM <- "#A3352A"; AMBAR <- "#9A5514"

secao("6. GRAFICO - o VIF de cada termo contra os limiares")

tbg <- tb[order(tb$equivale_a_VIF), ]
tbg$termo <- factor(tbg$termo, levels = tbg$termo)

g1 <- ggplot(tbg, aes(x = equivale_a_VIF, y = termo)) +
  annotate("rect", xmin = 5, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = VERM, alpha = .06) +
  geom_vline(xintercept = 1, colour = "grey55", linetype = "dashed") +
  geom_vline(xintercept = 5, colour = AMBAR, linewidth = .7) +
  geom_vline(xintercept = 10, colour = VERM, linewidth = .7) +
  annotate("text", x = 5, y = Inf, label = "5 = atencao", hjust = -.08, vjust = 1.6,
           colour = AMBAR, size = 3.2) +
  annotate("text", x = 10, y = Inf, label = "10 = severa", hjust = -.08, vjust = 1.6,
           colour = VERM, size = 3.2) +
  geom_segment(aes(x = 1, xend = equivale_a_VIF, yend = termo),
               colour = VERDE, linewidth = 1.2) +
  geom_point(size = 3.4, colour = VERDE) +
  geom_text(aes(label = sprintf("%.2f", equivale_a_VIF)), hjust = -.45, size = 3.2) +
  scale_x_continuous(limits = c(0.9, 12), breaks = c(1, 5, 10)) +
  labs(x = "VIF (para fatores: GVIF ajustado ao quadrado)", y = NULL,
       title = "Multicolinearidade: nenhum termo chega perto do limiar",
       subtitle = "maior valor 1,62 | limiares do script da aula de correlacao do professor",
       caption = "Docs/script_correlacao_aula_ver00.R secao 06.1") +
  theme_bw(base_size = 12)
print(g1); ggsave(file.path(FIG, "05_vif.png"), g1, width = 9, height = 6, dpi = 150, bg = "white")
cat("  figura salva: 05_vif.png\n")

cat("\n# Parecer\n")
cat("Todos os termos ficam na faixa entre 1,04 e 1,62 - colados no 1, que e a\n")
cat("ausencia total de inflacao. A cultivar, que e a variavel de interesse,\n")
cat("tem o MENOR VIF da tabela: praticamente ortogonal as demais.\n")

secao("7. GRAFICO - o estrago, quando ha multicolinearidade de verdade")

d <- data.frame(
  cenario = factor(c("sem a copia", "com a copia"), levels = c("sem a copia", "com a copia")),
  beta = c(coef(m_ok)[2], coef(m_bad)[2]),
  ep   = c(summary(m_ok)$coefficients[2,2], summary(m_bad)$coefficients[2,2]))
d$li <- d$beta - 1.96*d$ep; d$ls <- d$beta + 1.96*d$ep
v_ok <- vif(m_ok); v_bad <- vif(m_bad)
vif_ok_umidade <- if (is.matrix(v_ok)) v_ok["umidade_relativa_pct", 1] else v_ok["umidade_relativa_pct"]
vif_bad_umidade <- if (is.matrix(v_bad)) v_bad["umidade_relativa_pct", 1] else v_bad["umidade_relativa_pct"]
d$rot <- sprintf("VIF %.0f", c(vif_ok_umidade, vif_bad_umidade))

g2 <- ggplot(d, aes(x = cenario, y = beta)) +
  geom_hline(yintercept = 0, colour = "grey55", linetype = "dashed") +
  geom_errorbar(aes(ymin = li, ymax = ls), width = .12, linewidth = .9, colour = VERDE) +
  geom_point(size = 4, colour = VERDE) +
  geom_text(aes(label = rot), hjust = -0.35, size = 3.6, colour = "grey30") +
  labs(x = NULL, y = "coeficiente da umidade (IC 95%)",
       title = "O que a multicolinearidade faz com o coeficiente",
       subtitle = "a estimativa quase nao muda; o intervalo explode e cruza o zero",
       caption = "copia artificial da umidade com r = 0,9974 | R2 do modelo praticamente inalterado") +
  theme_bw(base_size = 12)
print(g2); ggsave(file.path(FIG, "05_estrago_multicolinearidade.png"), g2,
                  width = 8, height = 5.5, dpi = 150, bg = "white")
cat("  figura salva: 05_estrago_multicolinearidade.png\n")

cat("\n# Parecer\n")
cat("O ponto quase nao se move, mas a barra de incerteza passa a cruzar o zero.\n")
cat("E a assinatura da multicolinearidade: nao erra a estimativa, destroi a\n")
cat("precisao dela. Por isso e problema de INFERENCIA e nao de predicao.\n")
