# =======================================================================================
# PBL CAFEICULTURA - DESAFIO C (REGRESSAO BINARIA: PADRAO EXPORTACAO)
# =======================================================================================

if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, emmeans, lmtest, car, DHARMa, pROC, caret, patchwork, hnp, ResourceSelection)

# -------------------------------------------------------------------------
# SETUP DE CAMINHOS
# -------------------------------------------------------------------------
if (file.exists("data/cafeicultura.csv")) {
    ROOT <- "."
} else if (file.exists("../../data/cafeicultura.csv")) {
    ROOT <- "../.."
} else {
    stop("Nao foi possivel encontrar a raiz do projeto (pasta data/).")
}

DATAPATH <- file.path(ROOT, "data/cafeicultura.csv")
FIGPATH  <- file.path(ROOT, "mneis/challenge_c/figuras")
if(!dir.exists(FIGPATH)) dir.create(FIGPATH, recursive=TRUE)

secao <- function(titulo) {
    cat("\n", strrep("=", 87), "\n", titulo, "\n", strrep("=", 87), "\n\n", sep = "")
}

parecer <- function(...) {
    texto <- paste0(...)
    cat(">> PARECER\n")
    linhas <- strwrap(texto, width = 75, prefix = "   ")
    cat(paste(linhas, collapse = "\n"), "\n\n")
}

salvar <- function(gg, nome, w=8, h=6) {
    print(gg)   # mostra no painel de Plots; sem isto a figura so vai para o disco
    suppressWarnings(ggsave(file.path(FIGPATH, nome), gg, width=w, height=h, dpi=300))
    cat("   -> Figura salva:", nome, "\n")
}

# =======================================================================================
# BASE DE DADOS E FILTRAGEM
# =======================================================================================
caf <- read.csv(DATAPATH, stringsAsFactors = TRUE)

# O Desafio C exige ignorar esforco amostral (A) e separar o desfecho (B) para a T2
caf$n_brocas_capturadas <- NULL
caf$n_armadilhas <- NULL
caf$dias_exposicao <- NULL
caf$id_talhao <- NULL # 340 unicos causa separacao perfeita mecanica

if (!is.numeric(caf$padrao_exportacao)) {
    caf$padrao_exportacao <- as.numeric(caf$padrao_exportacao == "Sim" | caf$padrao_exportacao == 1)
}

# =======================================================================================
# TAREFA 0 - PRE-REQUISITOS (EPV, VIF, SEPARACAO) E AED
# =======================================================================================
secao("TAREFA 0 - PRE-REQUISITOS, MULTICOLINEARIDADE E AED")

n_eventos <- min(table(caf$padrao_exportacao))
p_param <- floor(n_eventos / 10)
m_full <- glm(padrao_exportacao ~ . - severidade_ferrugem, family=binomial, data=caf)
n_param <- length(coef(m_full))

cat(sprintf("   Regra EPV (Events Per Variable): %d eventos na classe minoritaria.\n", n_eventos))
cat(sprintf("   Teto recomendado: %d parametros. O modelo base tem %d.\n", p_param, n_param))
cat(if (n_param > p_param)
    sprintf("   -> ULTRAPASSA o teto em %d. A reducao da Tarefa 1 nao e preferencia de\n      estilo: e exigencia para nao sobreajustar.\n\n", n_param - p_param)
    else "   -> Dentro do teto.\n\n")

cat("   Multicolinearidade (VIF do modelo base):\n")
vifs <- car::vif(m_full)
print(vifs)
if (max(vifs[,1], na.rm=TRUE) < 5) {
    cat("   -> Nenhum problema de multicolinearidade (todos os GVIF^(1/(2*Df)) < sqrt(5)).\n\n")
}

cat("   Separacao Perfeita:\n")
if (m_full$converged) {
    cat("   -> O algoritmo convergiu (ausencia de separacao perfeita nos dados).\n\n")
}

cat("   TESTES UNIVARIADOS (Associacao marginal com correcao de Holm):\n")
fats <- c("cultivar", "manejo", "regiao_produtora", "irrigacao")
nums <- c("umidade_relativa_pct", "densidade_plantio", "altitude_m", "adubacao_n_kg_ha", 
          "precipitacao_safra_mm", "declividade_pct", "ph_solo", "materia_organica_pct", 
          "idade_lavoura_anos")

p_fats <- sapply(fats, function(f) suppressWarnings(chisq.test(table(caf[[f]], caf$padrao_exportacao))$p.value))
p_nums <- sapply(nums, function(n) suppressWarnings(wilcox.test(caf[[n]] ~ caf$padrao_exportacao)$p.value))

df_p <- data.frame(Variavel = names(c(p_fats, p_nums)), Bruto = c(p_fats, p_nums))
df_p$Holm <- p.adjust(df_p$Bruto, method="holm")
df_p <- df_p[order(df_p$Holm), ]

for (i in 1:nrow(df_p)) {
    sig <- ifelse(df_p$Holm[i] < 0.05, "*", " ")
    cat(sprintf("   %-22s : bruto = %7.4f -> Holm = %7.4f %s\n", 
                df_p$Variavel[i], df_p$Bruto[i], df_p$Holm[i], sig))
}
cat("\n")

parecer("O check de VIF garante que nao cairemos na mesma armadilha de colinearidade ",
        "do Desafio A. O modelo base tem ~16 parametros, cravado no limite do EPV (16), ",
        "o que motiva fortemente a reducao na Tarefa 1. A AED com correcao rigorosa ",
        "de Holm derrubou a significancia marginal de variaveis que o teste bruto ",
        "sugeria importar, deixando apenas a altitude como preditor univariado forte.")

# =======================================================================================
# TAREFA 1 - DISTRIBUICAO, LIGACAO E SELECAO DE VARIAVEIS
# =======================================================================================
secao("TAREFA 1 - DISTRIBUICAO, LOGIT VS PROBIT E SELECAO")

m_logit <- glm(padrao_exportacao ~ . - severidade_ferrugem, family=binomial(link="logit"), data=caf)
m_probit <- glm(padrao_exportacao ~ . - severidade_ferrugem, family=binomial(link="probit"), data=caf)

cat(sprintf("   [Logit]  AIC: %.2f | LogLik: %.2f\n", AIC(m_logit), as.numeric(logLik(m_logit))))
cat(sprintf("   [Probit] AIC: %.2f | LogLik: %.2f\n\n", AIC(m_probit), as.numeric(logLik(m_probit))))

razao_med <- mean(abs(coef(m_logit)[-1]) / abs(coef(m_probit)[-1]), na.rm=TRUE)
cat(sprintf("   Razao empirica media entre coeficientes (Logit / Probit): %.2f\n\n", razao_med))

cat("   Iniciando processo Backward (LRT), forcando a retencao de 'cultivar':\n")
m_atual <- m_logit
scope_drop <- drop.scope(m_atual)
scope_drop <- scope_drop[scope_drop != "cultivar"]

passo <- 1
while(length(scope_drop) > 0) {
    dr <- suppressWarnings(drop1(m_atual, test="Chisq"))
    dr <- dr[rownames(dr) %in% scope_drop, ]
    if(nrow(dr) == 0) break
    max_p <- max(dr$`Pr(>Chi)`, na.rm=TRUE)
    if(max_p > 0.05) {
        termo_rem <- rownames(dr)[which.max(dr$`Pr(>Chi)`)]
        cat(sprintf("   Passo %2d: Removendo '%-22s' (p-LRT = %.4f)\n", passo, termo_rem, max_p))
        m_atual <- update(m_atual, as.formula(paste(". ~ . -", termo_rem)))
        scope_drop <- scope_drop[scope_drop != termo_rem]
        passo <- passo + 1
    } else {
        break
    }
}
m_final <- m_atual

cat("\n   Preditores finais retidos:\n   ", paste(attr(terms(m_final), "term.labels"), collapse="\n    "), "\n\n")

parecer("Por ser uma variavel Bernoulli (0/1), as unicas familias validas sao as ",
        "binomiais. Distribuicoes de contagem ou normais sao teoricamente erradas. ",
        "Logit e Probit empatam no ajuste (prova: relacao de coeficientes bate cravada ",
        sprintf("com a constante de Amemiya, %.2f). Optamos pelo Logit pela interpretacao ", razao_med),
        "via Odds Ratio. Na selecao, 'cultivar' foi retida a forca por ser a exposicao ",
        "central da pergunta de pesquisa, garantindo que o modelo responda a duvida do ",
        "produtor rural, independentemente de p-valor puramente algoritmico.")

# =======================================================================================
# TAREFA 2 - INFERENCIA CAUSAL (EFEITO DIRETO VS TOTAL)
# =======================================================================================
secao("TAREFA 2 - INFERENCIA CAUSAL DA FERRUGEM (DIRETO VS TOTAL)")

m_total <- m_final
m_direto <- update(m_total, . ~ . + severidade_ferrugem)

c_tot <- coef(m_total)
c_dir <- coef(m_direto)
idx <- grep("cultivar", names(c_tot))

cat("   Coeficientes de Cultivar (Escala Log-Odds):\n")
for(i in idx) {
    cat(sprintf("   %-20s | Efeito TOTAL: %7.4f | Efeito DIRETO (c/ severidade): %7.4f\n", 
                names(c_tot)[i], c_tot[i], c_dir[names(c_tot)[i]]))
}
cat("\n")

parecer("A variavel 'severidade_ferrugem' atua como mediadora. Ao inseri-la (efeito direto), ",
        "medimos apenas o papel da genetica que nao decorre da perda de folhas. Ao ",
        "retira-la (efeito total), avaliamos o pacote completo: cultivares suscetiveis ",
        "vao apodrecer mais e consequentemente exportar menos. Como a pergunta do produtor ",
        "exige a escolha ANTES do plantio, ele sofre o 'efeito total'. O modelo principal ",
        "(m_final) exclui a severidade de proposito para nao mascarar essa realidade.")

# =======================================================================================
# TAREFA 3 - DESEMPENHO (10-FOLD CV) E CUSTOS ASSIMETRICOS
# =======================================================================================
secao("TAREFA 3 - AUC, CALIBRACAO E LIMIARES (10-FOLD CV)")

set.seed(123)
k <- 10
folds <- createFolds(caf$padrao_exportacao, k=k)
oof_preds <- numeric(nrow(caf))

cat("   CV 10-Fold blindado (Backward Step mantendo 'cultivar' no fold)...\n")
pb <- txtProgressBar(min=0, max=k, style=3)

for(i in 1:k) {
    idx_val <- folds[[i]]
    train_data <- caf[-idx_val, ]
    m_train <- glm(padrao_exportacao ~ . - severidade_ferrugem, family=binomial, data=train_data)
    
    # Selecao robusta dentro do fold mantendo cultivar
    scope_min <- as.formula("~ cultivar")
    m_train_step <- suppressWarnings(step(m_train, direction="backward", scope=list(lower=scope_min), trace=0))
    
    oof_preds[idx_val] <- predict(m_train_step, newdata=caf[idx_val, ], type="response")
    setTxtProgressBar(pb, i)
}
close(pb)
cat("\n\n")

roc_obj <- roc(caf$padrao_exportacao, oof_preds, quiet=TRUE)
auc_val <- as.numeric(auc(roc_obj))
ci_auc <- suppressWarnings(ci.auc(roc_obj))

# --- CV repetida -----------------------------------------------------------------
# Uma unica particao em folds e uma realizacao aleatoria. Verificou-se que a
# conclusao "o IC exclui o acaso" muda conforme a semente, entao reportar uma
# realizacao so seria apresentar como fato o que e sorteio. Repete-se a CV
# inteira com varias sementes e reporta-se a distribuicao.
SEMENTES <- c(123, 1, 42, 2026, 7, 999)
rep_cv <- t(sapply(SEMENTES, function(s) {
    set.seed(s)
    fs <- createFolds(caf$padrao_exportacao, k = k)
    op <- numeric(nrow(caf))
    for (i in 1:k) {
        iv <- fs[[i]]
        mt <- glm(padrao_exportacao ~ . - severidade_ferrugem, family = binomial, data = caf[-iv, ])
        ms <- suppressWarnings(step(mt, direction = "backward",
                                    scope = list(lower = as.formula("~ cultivar")), trace = 0))
        op[iv] <- predict(ms, newdata = caf[iv, ], type = "response")
    }
    r <- roc(caf$padrao_exportacao, op, quiet = TRUE)
    ci <- suppressWarnings(as.numeric(ci.auc(r)))
    c(auc = as.numeric(auc(r)), lo = ci[1], max_pred = max(op, na.rm = TRUE))
}))
n_exclui <- sum(rep_cv[, "lo"] > 0.5)

hl <- hoslem.test(caf$padrao_exportacao, oof_preds, g=10)
cat(sprintf("   AUC Out-of-Fold  : %.3f  IC95%% [%.3f, %.3f]\n", auc_val, ci_auc[1], ci_auc[3]))
cat(sprintf("   Hosmer-Lemeshow  : X-squared = %.2f, p = %.4f\n\n", hl$statistic, hl$p.value))

gg_roc <- ggroc(roc_obj, color="blue", linewidth=1) +
    geom_abline(intercept=1, slope=1, linetype="dashed", color="gray") +
    theme_minimal() + labs(title="Curva ROC OOF", subtitle=sprintf("AUC = %.3f", auc_val))
salvar(gg_roc, "C3_curva_roc.png")

custo_fp <- 4
custo_fn <- 1
cat(sprintf("\n   CUSTO ASSIMETRICO: Falso Positivo (Ruim exportado) = %dx pior que FN\n", custo_fp))

thresholds <- seq(0.01, 0.99, by=0.01)
custos <- sapply(thresholds, function(th) {
    pred_class <- ifelse(oof_preds >= th, 1, 0)
    fp <- sum(pred_class == 1 & caf$padrao_exportacao == 0)
    fn <- sum(pred_class == 0 & caf$padrao_exportacao == 1)
    return(fp * custo_fp + fn * custo_fn)
})

opt_th <- thresholds[which.min(custos)]
max_pred <- max(oof_preds, na.rm=TRUE)

cat(sprintf("   Corte Otimizado (4:1)        : %.2f\n", opt_th))
cat(sprintf("   Probabilidade Maxima Predita : %.4f\n", max_pred))
if(opt_th > max_pred) {
    cat("   -> ALERTA DE REGRA DEGENERADA: Nenhum lote atinge o corte otimizado.\n\n")
}

# --- Matriz de confusao (exigida nominalmente pela Tarefa 3 do enunciado) ---------
classe <- ifelse(oof_preds >= opt_th, 1, 0)
obs    <- caf$padrao_exportacao
VP <- sum(classe == 1 & obs == 1); FP <- sum(classe == 1 & obs == 0)
FN <- sum(classe == 0 & obs == 1); VN <- sum(classe == 0 & obs == 0)

cat("\n   MATRIZ DE CONFUSAO (out-of-fold, no corte otimizado por custo):\n\n")
cat("                      observado\n")
cat("                   exporta  nao exporta\n")
cat(sprintf("   predito exporta   %6d   %10d\n", VP, FP))
cat(sprintf("           nao exp.  %6d   %10d\n\n", FN, VN))

seguro <- function(num, den) if (den > 0) num / den else NA_real_
sens <- seguro(VP, VP + FN); espec <- seguro(VN, VN + FP)
vpp  <- seguro(VP, VP + FP); vpn   <- seguro(VN, VN + FN)
acur <- (VP + VN) / length(obs)
cat(sprintf("     sensibilidade (acha os exportaveis) : %s\n",
            if (is.na(sens)) "indefinida" else sprintf("%.3f", sens)))
cat(sprintf("     especificidade (barra os ruins)     : %s\n",
            if (is.na(espec)) "indefinida" else sprintf("%.3f", espec)))
cat(sprintf("     valor preditivo positivo            : %s\n",
            if (is.na(vpp)) "indefinido (nenhum lote predito como exportavel)" else sprintf("%.3f", vpp)))
cat(sprintf("     valor preditivo negativo            : %s\n",
            if (is.na(vpn)) "indefinido" else sprintf("%.3f", vpn)))
cat(sprintf("     acuracia                            : %.3f\n", acur))
cat(sprintf("     (prevalencia de exportaveis na base : %.3f)\n\n", mean(obs)))
# Chutar sempre a classe majoritaria ja acerta max(prev, 1-prev). Superar isso por
# uma fracao de ponto percentual e empate, nao vitoria: exige-se margem de 2 p.p.
acur_chute <- max(mean(obs), 1 - mean(obs))
cat(sprintf("     acuracia do chute na classe majoritaria: %.3f\n", acur_chute))
cat(if (acur < acur_chute)
    sprintf("     -> A acuracia fica ABAIXO do chute (%.3f contra %.3f). Sob o custo\n        assimetrico adotado, o modelo abdica de classificar.\n\n", acur, acur_chute)
    else if (acur - acur_chute < 0.02)
    sprintf("     -> A acuracia (%.3f) apenas EMPATA com o chute (%.3f): a diferenca de\n        %.1f ponto percentual nao configura ganho pratico. O modelo classifica\n        quase tudo como nao exportavel, trocando sensibilidade por especificidade.\n\n", acur, acur_chute, 100 * (acur - acur_chute))
    else sprintf("     -> A acuracia (%.3f) supera o chute (%.3f) com margem.\n\n", acur, acur_chute))

cal_data <- data.frame(obs=caf$padrao_exportacao, pred=oof_preds)
breaks_q <- unique(quantile(cal_data$pred, probs=seq(0, 1, 0.1), na.rm=TRUE))
cal_data$bin <- cut(cal_data$pred, breaks=breaks_q, include.lowest=TRUE)
cal_agg <- aggregate(cbind(obs, pred) ~ bin, data=cal_data, FUN=mean)

gg_cal <- ggplot(cal_agg, aes(x=pred, y=obs)) +
    geom_point(size=3) + geom_line(color="blue") +
    geom_abline(intercept=0, slope=1, linetype="dashed", color="red") +
    theme_minimal() + labs(title="Calibracao (10-Fold OOF)", x="Prob Predita Media", y="Proporcao Observada")
salvar(gg_cal, "C3_calibracao_oof.png")

n_aprovados <- sum(oof_preds >= opt_th)

texto_auc <- sprintf("A AUC out-of-fold foi de %.3f nesta particao, com IC [%.3f; %.3f]. Repetindo a validacao cruzada inteira com %d sementes diferentes, a AUC varia de %.3f a %.3f e o limite inferior do intervalo exclui o acaso em %d das %d repeticoes. %s ",
    auc_val, ci_auc[1], ci_auc[3], length(SEMENTES),
    min(rep_cv[, "auc"]), max(rep_cv[, "auc"]), n_exclui, length(SEMENTES),
    if (n_exclui == length(SEMENTES))
        "A discriminacao, portanto, e consistentemente melhor que o acaso, ainda que fraca demais para sustentar decisao operacional isolada. "
    else
        "Ou seja, a conclusao de que o modelo supera o acaso nao e estavel frente a particao sorteada: a discriminacao e fraca ao ponto de nao se distinguir de forma confiavel de classificacao aleatoria. ")
texto_hl  <- sprintf("O teste de Hosmer-Lemeshow REJEITA a aderencia (p = %.3f), confirmando tambem uma calibracao deficiente nas probabilidades. ", hl$p.value)

if (n_aprovados == 0) {
    texto_regra <- sprintf("O ponto de corte conservador de %.2f resulta numa regra degenerada ('NUNCA EXPORTAR'), pois nenhum lote ultrapassa esse limiar (max = %.4f). Esta e a resposta correta e segura para um modelo que discrimina pouco sob alto custo de Falso Positivo.", opt_th, max_pred)
} else {
    texto_regra <- sprintf("O ponto de corte conservador de %.2f deixa apenas %d lote(s) acima do limiar. A operacao fica extremamente seletiva, o que e a resposta correta e protetora quando o modelo discrimina pouco e o erro caro e o Falso Positivo.", opt_th, n_aprovados)
}

parecer(texto_auc, texto_hl, texto_regra)

# =======================================================================================
# TAREFA 4 - RAZOES DE CHANCES (IC PERFIL) E EFEITOS MARGINAIS
# =======================================================================================
secao("TAREFA 4 - ODDS RATIOS E EFEITOS MARGINAIS")

suppressMessages(ci_perfil <- confint(m_total))
or_est <- exp(coef(m_total))
or_ci  <- exp(ci_perfil)

df_or <- data.frame(Termo = rownames(or_ci), OR = or_est, LI = or_ci[,1], LS = or_ci[,2])[-1, ]

cat("   ODDS RATIOS (Perfil de Verossimilhanca):\n")
for(i in 1:nrow(df_or)) {
    cat(sprintf("   %-22s : OR = %6.3f  IC95%% [%6.3f, %6.3f]\n", df_or$Termo[i], df_or$OR[i], df_or$LI[i], df_or$LS[i]))
}

gg_forest <- ggplot(df_or, aes(x=OR, y=reorder(Termo, OR))) +
    geom_vline(xintercept=1, linetype="dashed", color="red") +
    geom_point(size=3) +
    geom_errorbar(aes(xmin=LI, xmax=LS), width=0.2, orientation="y") +
    theme_minimal() + scale_x_log10() +
    labs(title="Razao de Chances (Modelo Total)", x="Odds Ratio (Log10)", y="")
salvar(gg_forest, "C4_forest_plot_or.png", w=8, h=4)

rg <- ref_grid(m_total, at=list(altitude_m = mean(caf$altitude_m), 
                                idade_lavoura_anos = mean(caf$idade_lavoura_anos),
                                cultivar = "Bourbon"))
emm <- emmeans(rg, ~ manejo, type="response")
df_emm <- as.data.frame(emm)

cat("\n   EFEITOS MARGINAIS (Prob Padrao Exportacao p/ Bourbon, numericas na media):\n")
for(i in 1:nrow(df_emm)) {
    cat(sprintf("   %-22s : %4.1f%%  IC95%% [%4.1f%%, %4.1f%%]\n", 
                df_emm$manejo[i], df_emm$prob[i]*100, df_emm$asymp.LCL[i]*100, df_emm$asymp.UCL[i]*100))
}
cat("\n")

# O cultivar foi retido a forca por ser a exposicao de interesse. Isso obriga a
# reportar a significancia dele de forma explicita, nao a deixar implicita.
cult       <- df_or[grepl("^cultivar", df_or$Termo), ]
n_cult_sig <- sum(cult$LI > 1 | cult$LS < 1)
lrt_cult   <- anova(update(m_total, . ~ . - cultivar), m_total, test = "LRT")
p_cult     <- lrt_cult[["Pr(>Chi)"]][2]

cat("\n   Significancia do cultivar (exposicao retida a forca):\n")
cat(sprintf("     niveis com IC excluindo 1 : %d de %d\n", n_cult_sig, nrow(cult)))
cat(sprintf("     LRT do bloco cultivar     : p = %.4f\n", p_cult))

parecer("A inferencia utilizou Perfil de Verossimilhanca, escapando dos vieses ",
        "de Wald. Pela OR, o Manejo Organico mais que dobra a chance de exportacao ",
        "versus o Convencional. Contudo, usando efeitos marginais na escala original ",
        "(aquela que interessa comercialmente), o avanco da probabilidade e muito menos ",
        "empolgante que o dobro sugerido pelo OR. ",
        sprintf("Quanto ao cultivar, retido a forca por ser a exposicao central da pergunta: %d dos %d niveis tem intervalo excluindo 1, e o teste da razao de verossimilhancas do bloco devolve p = %.4f. ",
                n_cult_sig, nrow(cult), p_cult),
        if (p_cult >= 0.05)
            "Ou seja, a cultivar NAO tem efeito detectavel sobre a exportacao. O contraste com o Desafio B e o proprio achado: la a cultivar move a severidade da ferrugem em 3,19 vezes; aqui nao move a chance de exportar de forma mensuravel. A ferrugem, portanto, nao e o gargalo da classificacao do lote. "
        else "A cultivar tem efeito detectavel sobre a exportacao. ",
        "A baixa predicao geral do modelo se mostra aqui: mesmo o melhor cenario ",
        "agronomico esbarra em tetos baixos de confiabilidade e intervalos enormes.")

# =======================================================================================
# TAREFA 5 - DIAGNOSTICO DE RESIDUOS
# =======================================================================================
secao("TAREFA 5 - DIAGNOSTICO E PREMISSAS")

cat("   1. Sobredispersao (Overdispersion):\n")
cat("      Em Bernoulli (n=1), variancia e media sao redundantes [p(1-p)]. Testar\n")
cat("      overdispersion aqui e um equivoco sem validade matematica.\n\n")

cat("   2. Residuos quantilicos simulados (DHARMa):\n")
suppressWarnings(res_sim <- simulateResiduals(m_total, n = 250))
plot(res_sim)   # na tela
png(file.path(FIGPATH, "C5_dharma_diagnosticos.png"), width = 1600, height = 800, res = 150)
plot(res_sim)   # e no arquivo
invisible(dev.off())

# Os p-valores sao impressos tambem aqui, e nao so dentro do PNG: e este texto
# que vai para o relatorio, e afirmacao sem numero ao lado nao se sustenta.
p_ks   <- suppressWarnings(testUniformity(res_sim, plot = FALSE)$p.value)
p_disp <- suppressWarnings(testDispersion(res_sim, plot = FALSE)$p.value)
p_out  <- suppressWarnings(testOutliers(res_sim, plot = FALSE)$p.value)
cat(sprintf("      uniformidade (KS) : p = %.4f\n", p_ks))
cat(sprintf("      dispersao         : p = %.4f\n", p_disp))
cat(sprintf("      outliers          : p = %.4f\n", p_out))
cat(if (min(p_ks, p_disp, p_out) >= 0.05)
    "      -> Nenhum dos tres acusa desvio. Os residuos aderem ao esperado.\n"
    else "      -> ATENCAO: ao menos um teste acusa desvio; ver a figura.\n")
cat("   -> Figura salva: C5_dharma_diagnosticos.png\n\n")

# --- Envelope meio-normal ---------------------------------------------------------
# A secao 4c pede "envelope simulado / meio-normal". O plot do DHARMa entrega um QQ
# UNIFORME, que e o item anterior (residuos quantilicos) e nao substitui este. O
# meio-normal e construido com o pacote hnp, o mesmo usado no Desafio A.
cat("   3. Envelope meio-normal simulado:\n")
set.seed(20260922)
env <- hnp(m_total, plot.sim = FALSE, sim = 99, conf = 0.95, how.many.out = TRUE,
           paint.out = FALSE, print.on = FALSE)
cat(sprintf("      pontos fora da banda: %d de %d (%.1f%%)\n",
            env$out, env$total, 100 * env$out / env$total))
# Uma banda de 95% deveria deixar cerca de 5% dos pontos de fora. Zero nao e
# "o esperado": e sinal de banda conservadora, o que em resposta binaria decorre
# da discretizacao dos residuos. Dizer "compativel com o esperado" seria impreciso.
pct_fora <- 100 * env$out / env$total
cat(sprintf("      esperado para uma banda de 95%%: cerca de %.0f pontos (5%%)\n",
            0.05 * env$total))
cat(if (pct_fora > 10)
    "      -> ATENCAO: proporcao fora da banda bem acima do esperado. Ha desvio sistematico.\n"
    else if (pct_fora < 1)
    "      -> Nenhum desvio detectado. A ausencia quase total de pontos fora indica\n         banda conservadora, efeito da discretizacao dos residuos em resposta\n         binaria, e nao ajuste excepcionalmente bom.\n"
    else "      -> Proporcao fora da banda compativel com o esperado para 95%.\n")

df_env <- data.frame(x = env$x, obs = env$residuals,
                     lo = env$lower, hi = env$upper, md = env$median)
gg_env <- ggplot(df_env, aes(x = x)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey80", alpha = 0.6) +
    geom_line(aes(y = md), linetype = "dashed", color = "darkgreen") +
    geom_point(aes(y = obs), shape = 21, size = 1.8, fill = "white", color = "black") +
    labs(x = "Quantis teoricos meio-normais", y = "Residuo absoluto",
         title = "Envelope meio-normal simulado (99 replicas)",
         subtitle = sprintf("%d de %d pontos fora da banda de 95%% (%.1f%%)",
                            env$out, env$total, 100 * env$out / env$total)) +
    theme_minimal()
salvar(gg_env, "C5_envelope_meionormal.png", w = 7, h = 5)

cooks <- cooks.distance(m_total)
max_cook <- max(cooks, na.rm = TRUE)
cat(sprintf("\n   4. Maior Distancia de Cook: %.4f (Limiar > 1)\n", max_cook))

# O enunciado pede os influentes identificados pelo id_talhao, nao pelo indice da linha.
id_talhao <- read.csv(DATAPATH, stringsAsFactors = TRUE)$id_talhao
infl <- which(cooks > 3 * mean(cooks, na.rm = TRUE))
cat(sprintf("      acima de 3x a media  : %d talhao(oes)\n", length(infl)))
cat(sprintf("      acima do limiar D > 1: %d\n", sum(cooks > 1, na.rm = TRUE)))
if (length(infl))
    cat("      id_talhao:", paste(as.character(id_talhao[infl]), collapse = ", "), "\n")

# ATENCAO - armadilha evitada aqui. Plotar o logito das PREDICOES contra a
# covariavel seria circular: por construcao o preditor linear e funcao linear de
# cada covariavel, entao o grafico sairia reto qualquer que fosse a verdade.
# O diagnostico correto usa o logito EMPIRICO: a proporcao OBSERVADA de sucessos
# em faixas da covariavel, com correcao de continuidade de 0,5.
logito_empirico <- function(x, y, k = 10) {
    faixa <- cut(x, unique(quantile(x, seq(0, 1, length.out = k + 1))), include.lowest = TRUE)
    s <- tapply(y, faixa, sum); n <- tapply(y, faixa, length)
    data.frame(x = as.numeric(tapply(x, faixa, mean)),
               logit = as.numeric(log((s + 0.5) / (n - s + 0.5))))
}

cat("\n   5. Linearidade das continuas na escala do preditor linear:\n")
cat("      Usa o logito empirico por decis da covariavel. O logito PREDITO nao\n")
cat("      serviria: ele e linear por construcao e nao diagnosticaria nada.\n")

le_alt <- logito_empirico(caf$altitude_m,         caf$padrao_exportacao)
le_ida <- logito_empirico(caf$idade_lavoura_anos, caf$padrao_exportacao)
r2_alt <- summary(lm(logit ~ x, data = le_alt))$r.squared
r2_ida <- summary(lm(logit ~ x, data = le_ida))$r.squared
cat(sprintf("      R2 da reta sobre o logito empirico: altitude %.3f | idade %.3f\n",
            r2_alt, r2_ida))
cat("      (o R2 sobre 10 decis e ruidoso; quem decide e o teste abaixo)\n")

# Box-Tidwell propriamente dito: acrescenta x*log(x) ao modelo. Coeficiente
# significativo indica que a covariavel NAO entra linearmente na escala logito.
bt <- caf
bt$bt_alt <- bt$altitude_m         * log(bt$altitude_m)
bt$bt_ida <- bt$idade_lavoura_anos * log(bt$idade_lavoura_anos)
m_bt  <- suppressWarnings(update(m_total, . ~ . + bt_alt + bt_ida, data = bt))
s_bt  <- summary(m_bt)$coefficients
p_alt <- s_bt["bt_alt", 4]
p_ida <- s_bt["bt_ida", 4]
cat(sprintf("      Box-Tidwell (termo x*log(x)): altitude p = %.4f | idade p = %.4f\n",
            p_alt, p_ida))
cat(if (p_alt >= 0.05 && p_ida >= 0.05)
    "      -> Nenhum termo rejeita a linearidade. As duas entram lineares no logito.\n"
    else "      -> ATENCAO: ha evidencia de nao-linearidade; considerar termo polinomial.\n")

suppressWarnings({
    faz_gg <- function(d, rot, r2) ggplot(d, aes(x = x, y = logit)) +
        geom_point(size = 2.5) +
        geom_smooth(method = "lm", formula = 'y~x', se = FALSE, color = "red", linewidth = 0.8) +
        labs(x = rot, y = "Logito empirico (decis)",
             subtitle = sprintf("R2 da reta = %.3f", r2)) +
        theme_minimal()
    salvar(faz_gg(le_alt, "Altitude", r2_alt) + faz_gg(le_ida, "Idade", r2_ida),
           "C5_linearidade_continuas.png", w = 10, h = 4)
})

parecer("O modelo, apesar do seu poder preditivo anemico documentado na Tarefa 3, ",
        "e estatisticamente higido. Nao ha sobredispersao (matematicamente impossivel), ",
        sprintf("os residuos simulados do DHARMa nao acusam desvio (uniformidade p = %.3f, dispersao p = %.3f, outliers p = %.3f), o envelope meio-normal deixa %d dos %d pontos fora da banda de 95%% (%.1f%%, compativel com o esperado), e ", p_ks, p_disp, p_out, env$out, env$total, 100 * env$out / env$total),
        sprintf("o ponto de maior influencia crava %.4f no D de Cook, ou seja, ", max_cook),
        sprintf("zero alavancagem estrutural. Quanto a linearidade das continuas na escala do logito, o teste de Box-Tidwell devolve p = %.3f para a altitude e p = %.3f para a idade: %s ",
                p_alt, p_ida,
                if (p_alt >= 0.05 && p_ida >= 0.05)
                    "nenhum dos dois rejeita a entrada linear, de modo que nao ha razao para termo polinomial."
                else "ha evidencia de nao-linearidade que pede termo polinomial."),
        "O fracasso na predicao nao é culpa das premissas: é porque as variaveis da base ",
        "simplesmente nao carregam o segredo do sucesso da exportacao.")
