# =======================================================================================
# UNIVERSIDADE DO VALE DO ITAJAI - UNIVALI
# ESCOLA POLITECNICA
# PROGRAMA DE POS-GRADUACAO EM COMPUTACAO APLICADA - PPGCA
# MESTRADO EM COMPUTACAO APLICADA
# Disciplina: Modelagem estatistica
# Prof. Dr.: Rodrigo Sant'Ana
# Discentes: Andre Lucas Ribeiro, Leticia Zorzi Rama, Matheus Neis
# Itajai, Santa Catarina, Brasil
#
# =======================================================================================
# DESAFIO C - REGRESSAO LOGISTICA (BINOMIAL)
# =====================================================================
# Pergunta do corpo tecnico:
# O que determina a chance de um lote atingir o padrao de exportacao?
#
# Variavel resposta: `padrao_exportacao`.
# Para cada talhao existe apenas um resultado:
#   0 = o lote nao atingiu o padrao de exportacao;
#   1 = o lote atingiu o padrao de exportacao.
#
# Estatisticamente, esse resultado individual 0/1 e chamado de
# observacao de Bernoulli. (A diferenca entre esse caso e o caso
# "k sucessos em n tentativas" e tratada como aprofundamento no
# bloco do modelo nulo.)
#
# Variaveis candidatas:
# `cultivar`, `manejo`, `irrigacao`, `regiao_produtora`, `altitude_m`,
# `declividade_pct`, `idade_lavoura_anos`, `densidade_plantio`,
# `adubacao_n_kg_ha`, `precipitacao_safra_mm`, `umidade_relativa_pct`,
# `ph_solo`, `materia_organica_pct`
#
# Deixadas de fora do conjunto de candidatas, com justificativa:
# - `id_talhao`: identificador do talhao, nao uma caracteristica dele.
# - `n_armadilhas`, `dias_exposicao`: descrevem o esforco de monitoramento
#   da broca, nao o talhao; nao ha razao agronomica para relaciona-los
#   com a qualidade do lote.
# - `n_brocas_capturadas`: e desfecho do Desafio A e depende do esforco
#   amostral (armadilhas x dias). Sera usada apenas numa verificacao de
#   sensibilidade no Step 4.
# - `severidade_ferrugem`: e desfecho do Desafio B e sera analisada
#   como possivel mediadora da relacao entre cultivar e padrao de
#   exportacao. Por isso, entra somente na Tarefa 2 (Step 6), onde a
#   hipotese de mediacao sera investigada, nao assumida.
#
# Tarefas do enunciado:
# 1. Ajustem o modelo logistico e conduzam selecao de variaveis.
# 2. Questao de raciocinio causal: `severidade_ferrugem` e ela propria um
#    desfecho (Desafio B). Ajustem o modelo COM e SEM essa variavel e
#    comparem os coeficientes de `cultivar`. Qual dos dois responde melhor
#    a pergunta do produtor que precisa ESCOLHER a cultivar antes do plantio?
# 3. Avaliem a capacidade preditiva: matriz de confusao, curva ROC e AUC,
#    discutindo a escolha do ponto de corte em funcao do custo assimetrico.
# 4. Interpretem em razoes de chances com intervalos de confianca.
#
# O script esta estruturado em steps cumulativos, seguindo o ciclo de
# trabalho apresentado nos slides da disciplina (secao 12.1):
# especificar - ajustar - selecionar - diagnosticar - interpretar - decidir.
# - Configurar ambiente de trabalho e funcoes auxiliares
# - STEP 0: dados, natureza da resposta e modelo nulo
# - STEP 1: exploracao condicional e logit empirico
# - STEP 2: especificacao - por que binomial e por que logit
# - STEP 3: modelo amplo, deviance, LRT, multicolinearidade
# - STEP 4: selecao de variaveis
# - STEP 5: diagnostico de residuos, influencia, ligacao e escala
# - STEP 6: efeito total x direto (severidade como mediadora)
# - STEP 7: razoes de chances com IC
# - STEP 8: capacidade preditiva e ponto de corte
# - STEP 9: cenarios e recomendacao
#
# Base de dados: data/cafeicultura.csv
# Figuras geradas: challenge_c/figuras/
#
# Convencao dos comentarios de cada bloco:
#   OBJETIVO / O QUE O CODIGO FAZ / POR QUE ISSO E NECESSARIO / COMO INTERPRETAR
# A explicacao principal vem sempre em linguagem simples. Detalhes
# matematicos ficam em blocos marcados [APROFUNDAMENTO], que podem ser
# pulados na primeira leitura sem prejuizo para acompanhar a analise.
# =====================================================================


# =====================================================================
# CONFIGURAR AMBIENTE DE TRABALHO
# =====================================================================

# Pacotes - pacman verifica se ja estao instalados e instala se preciso
if (!require(pacman)) {
    print(paste0("Pacote ainda nao instalado. Instalando..."))
    install.packages("pacman", dependencies = TRUE)
    library(pacman)
} else {
    print(paste0("Pacote ja instalado e carregado"))
    library(pacman)
}

# Pacotes usados ao longo dos steps (cada um e explicado no step em que
# entra em uso):
# ggplot2, tidyr, patchwork, dplyr : manipulacao de dados e graficos
# car     : vif() - medida de multicolinearidade (Step 3)
# MuMIn   : dredge() - ajuste de todos os subconjuntos de modelos (Step 4)
# pROC    : roc(), auc() - capacidade preditiva (Step 8)
# hnp     : envelope simulado meio-normal (Step 5)
# statmod : qresiduals() - residuos quantilicos aleatorizados (Step 5)
# mgcv    : gam() - verificacao de linearidade (Step 5)
# lmtest  : lrtest() - teste da razao de verossimilhancas (Steps 3 e 4)
#
# dplyr e carregado por ultimo de proposito. Quando dois pacotes exportam
# uma funcao com o mesmo nome, prevalece a do pacote carregado por ultimo.
# O hnp carrega o MASS, que tambem possui select(); carregando o dplyr
# depois, garantimos que filter(), select() e mutate() sejam os do dplyr.
p_load(ggplot2, tidyr, patchwork, car, MuMIn, pROC, hnp, statmod, mgcv, lmtest, dplyr)

# Opcoes gerais
# scipen = 10 : evita notacao cientifica (0.0025 em vez de 2.5e-03)
# warn   = 1  : imprime cada aviso no momento em que ocorre, e nao no fim
#               da execucao; assim sabemos qual bloco gerou o aviso
options(scipen = 10, warn = 1)

# Semente de reprodutibilidade.
# Necessaria sempre que houver sorteio: residuos quantilicos aleatorizados,
# envelope simulado, validacao cruzada. Sem ela, o resultado muda a cada
# execucao, e o enunciado exige script reprodutivel.
set.seed(20260920)

# ---------------------------------------------------------------------
# Localizacao da raiz do repositorio.
# O script encontra os dados sozinho, nao importa de onde seja executado:
# pasta aberta como projeto, arquivo solto com Source, Ctrl+Enter linha a
# linha, ou Rscript. Nao depende de setwd().
# ---------------------------------------------------------------------
localizar_script <- function() {
    a <- commandArgs(trailingOnly = FALSE)
    m <- grep("^--file=", a, value = TRUE)
    if (length(m))
        return(normalizePath(sub("^--file=", "", m[1]), winslash = "/", mustWork = FALSE))
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
        if (identical(pai, d))
            stop("Nao encontrei '", marcador, "'.\n",
                 "  Abra a PASTA do repositorio, nao o arquivo solto.", call. = FALSE)
        d <- pai
    }
}

RAIZ <- raiz_repo()
FIG  <- file.path(RAIZ, "challenge_c", "figuras")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------
# Funcoes auxiliares
# ---------------------------------------------------------------------

# Cabecalho de secao
secao <- function(texto) {
    barra <- paste(rep("=", 72), collapse = "")
    cat("\n", barra, "\n", texto, "\n", barra, "\n\n", sep = "")
}

# Parecer em prosa. Marcado de forma uniforme para poder ser extraido de
# uma vez so na hora de montar o relatorio:
#   grep -A20 ">> PARECER" saida.txt
# Regra adotada: nenhum numero e digitado a mao no texto; usa-se sempre
# sprintf() com valores calculados, para que o parecer continue correto
# se os dados mudarem.
parecer <- function(...) {
    txt <- paste0(...)
    cat("\n>> PARECER\n")
    cat(strwrap(txt, width = 72, prefix = "   "), sep = "\n")
    cat("\n")
}

# Salvar figura em disco E mostrar no painel de Plots
salvar <- function(g, nome, w = 9, h = 6) {
    print(g)
    ggsave(file.path(FIG, nome), g, width = w, height = h, dpi = 150, bg = "white")
    cat("   figura salva: ", nome, "\n", sep = "")
}

# Esta funcao impede a interpretacao de modelos que nao convergiram
# e destaca sinais de possivel separacao completa.
#
# Nao convergencia: o ajuste de um glm() e iterativo e pode nao chegar a
# uma resposta estavel. Mesmo assim, o R devolve um objeto cujo summary()
# imprime coeficientes e p-valores normalmente. Aqui, a nao convergencia
# interrompe o script com erro.
#
# Possivel separacao completa (slides da disciplina, secao 9.1.4): ocorre
# quando alguma combinacao de covariaveis separa perfeitamente os
# talhoes com 0 dos talhoes com 1. Nesse caso a estimativa de maxima
# verossimilhanca nao existe como numero finito.
#
# Sinais observados nesta funcao: probabilidades ajustadas numericamente
# iguais a 0 ou 1, e erros-padrao muito grandes (acima de 10). Sao sinais
# de alerta para investigacao, nao uma prova formal de separacao.
#
# Por que nao interpretar os coeficientes nessa situacao: as estimativas
# e os erros-padrao deixam de ter significado numerico, e as conclusoes
# extraidas deles seriam arbitrarias.
#
# [APROFUNDAMENTO] O limiar de probabilidade usado abaixo e o mesmo que
# o proprio glm.fit() usa para emitir o aviso "fitted probabilities
# numerically 0 or 1 occurred". O limiar de 10 para o erro-padrao e uma
# regra pratica deste projeto. Se o alerta aparecer, os slides indicam
# como alternativa a regressao logistica penalizada de Firth (pacote
# logistf) ou a remocao da covariavel responsavel.
verificar <- function(modelo, nome) {
    conv <- if (!is.null(modelo$converged)) modelo$converged else TRUE
    if (!isTRUE(conv))
        stop("O modelo '", nome, "' NAO convergiu. Nao interprete nada dele.",
             call. = FALSE)
    cat("   [ok] ", nome, " convergiu", sep = "")
    if (!is.null(modelo$iter)) cat(" em ", modelo$iter, " iteracoes", sep = "")
    cat("\n")

    # Sinais de possivel separacao completa: verificados apenas em glm binomial
    if (inherits(modelo, "glm") && identical(modelo$family$family, "binomial")) {
        eps      <- 10 * .Machine$double.eps          # mesmo limiar do glm.fit
        p_ajust  <- fitted(modelo)
        extremos <- sum(p_ajust < eps | p_ajust > 1 - eps)
        ep       <- summary(modelo)$coefficients[, "Std. Error"]
        ep_gig   <- sum(ep > 10)
        if (extremos > 0 || ep_gig > 0) {
            cat("   [ATENCAO] possivel SEPARACAO COMPLETA em '", nome, "':\n", sep = "")
            if (extremos > 0)
                cat("      ", extremos, " probabilidade(s) ajustada(s) exatamente 0 ou 1\n", sep = "")
            if (ep_gig > 0)
                cat("      ", ep_gig, " erro(s)-padrao acima de 10\n", sep = "")
            cat("      Nao interprete os coeficientes. Ver slides 9.1.4 (Firth / logistf).\n")
        }
    }
    invisible(modelo)
}

# Tabela de razoes de chances com intervalo de confianca (usada nos
# Steps 6 e 7). O intervalo e o de verossimilhanca perfilada, que e o
# que confint() devolve para um glm; e mais confiavel do que o de Wald
# em amostras moderadas e fica assimetrico na escala da razao de chances.
tabela_or <- function(modelo, nivel = 0.95) {
    ic <- suppressMessages(confint(modelo, level = nivel))
    b  <- coef(modelo)
    ep <- summary(modelo)$coefficients[, "Std. Error"]
    data.frame(termo  = names(b),
               beta   = unname(b),
               ep     = unname(ep),
               OR     = unname(exp(b)),
               IC_inf = unname(exp(ic[, 1])),
               IC_sup = unname(exp(ic[, 2])))
}

# Cores usadas nos graficos (mesmas do Desafio B)
VERDE <- "#18675A"; VERMELHO <- "#A3352A"; CINZA <- "grey55"


# =====================================================================
# ---- STEP 0 ---- BASE DE DADOS E VALIDACAO DA RESPOSTA
# =====================================================================
secao("STEP 0 - BASE DE DADOS E NATUREZA DA RESPOSTA")

# OBJETIVO: carregar a base com as variaveis categoricas ja reconhecidas
#   como fatores.
# O QUE O CODIGO FAZ: le o CSV. O argumento stringsAsFactors = TRUE converte
#   cada coluna de texto (cultivar, manejo, irrigacao, regiao) em fator,
#   com os niveis em ordem alfabetica.
# POR QUE ISSO E NECESSARIO: o glm() precisa tratar "Bourbon" como uma
#   categoria, e nao como texto livre. Alem disso, a ordem alfabetica dos
#   niveis define qual deles sera a categoria de referencia (ver adiante).
dados_cafe <- read.csv(file.path(RAIZ, "data", "cafeicultura.csv"),
                       stringsAsFactors = TRUE)
resposta_exportacao <- dados_cafe$padrao_exportacao

cat("observacoes:", nrow(dados_cafe), "| variaveis:", ncol(dados_cafe), "\n")

# ---------------------------------------------------------------------
# Validacao da resposta e dos fatores
# ---------------------------------------------------------------------
# OBJETIVO: confirmar que a resposta tem o formato que o modelo logistico
#   espera, antes de calcular qualquer coisa a partir dela.
# O QUE O CODIGO FAZ: confere, em ordem, (1) que a coluna existe, (2) que
#   so contem 0 e 1, (3) que nao ha valores ausentes (NA), (4) que as duas
#   classes aparecem, (5) que cada classe tem um minimo de observacoes,
#   (6) que nenhum fator tem nivel sem observacao. Se algo falhar, o
#   script para com uma mensagem explicando o problema.
# POR QUE ISSO E NECESSARIO: um valor fora de {0, 1}, uma classe ausente
#   ou um nivel de fator sem observacoes nao interrompem o glm(). Eles
#   produzem estimativas sem sentido ou coeficientes NA, que poderiam ser
#   interpretados por engano.
cat("\nverificacao da resposta:\n")

if (!"padrao_exportacao" %in% names(dados_cafe))
    stop("A coluna 'padrao_exportacao' nao existe na base.", call. = FALSE)
cat("   [ok] coluna padrao_exportacao existe\n")

valores_invalidos <- setdiff(unique(na.omit(resposta_exportacao)), c(0, 1))
if (length(valores_invalidos) > 0)
    stop("A resposta contem valores fora de {0, 1}: ",
         paste(valores_invalidos, collapse = ", "), call. = FALSE)
cat("   [ok] valores restritos a 0 e 1\n")

n_na_resposta <- sum(is.na(resposta_exportacao))
if (n_na_resposta > 0)
    stop("A resposta tem ", n_na_resposta, " NA. Decida o tratamento antes de seguir.",
         call. = FALSE)
cat("   [ok] nenhum NA na resposta\n")

if (!all(c(0, 1) %in% resposta_exportacao))
    stop("So uma das classes esta presente. Nao ha o que modelar.", call. = FALSE)
cat("   [ok] as duas classes (0 e 1) estao presentes\n")

# Minimo por classe. Uma classe com poucas observacoes pode produzir
# estimativas mais instaveis e intervalos de confianca mais largos.
# Neste projeto, 30 observacoes e apenas um limite pratico para emitir
# um alerta, nao uma garantia de que o modelo sera adequado: atingir
# esse numero nao dispensa os diagnosticos dos steps seguintes.
MINIMO_POR_CLASSE <- 30
menor_classe <- min(table(resposta_exportacao))
if (menor_classe < MINIMO_POR_CLASSE)
    warning("A classe menor tem so ", menor_classe, " observacoes (< ",
            MINIMO_POR_CLASSE, "). Inferencia fragil.", call. = FALSE)
cat("   [ok] classe menor com ", menor_classe, " observacoes (minimo declarado: ",
    MINIMO_POR_CLASSE, ")\n", sep = "")

# Fatores. Um nivel registrado no fator mas sem nenhuma observacao nao
# gera erro no glm(); gera um coeficiente NA no summary(). Se isso
# ocorrer, droplevels() remove o nivel vazio e o script avisa.
fatores <- names(dados_cafe)[sapply(dados_cafe, is.factor) & names(dados_cafe) != "id_talhao"]
for (f in fatores) {
    contagem <- table(dados_cafe[[f]])
    if (any(contagem == 0)) {
        warning("O fator '", f, "' tem nivel sem observacoes: ",
                paste(names(contagem)[contagem == 0], collapse = ", "),
                ". Removido com droplevels().", call. = FALSE)
        dados_cafe[[f]] <- droplevels(dados_cafe[[f]])
    }
}
cat("   [ok] nenhum fator com nivel vazio (", length(fatores), " fatores conferidos)\n", sep = "")


# =====================================================================
# ---- STEP 0 ---- QUANTOS LOTES ATINGEM O PADRAO?
# =====================================================================
# OBJETIVO: responder a pergunta mais simples possivel: de cada N talhoes,
#   quantos produziram lote com padrao de exportacao?
# O QUE O CODIGO FAZ: conta os 0s e os 1s e divide o numero de 1s pelo
#   total. Como a resposta so vale 0 ou 1, mean() faz exatamente essa
#   divisao: a soma dos valores conta os 1s, e dividir por n da a fracao.
# POR QUE ISSO E NECESSARIO: essa proporcao e o ponto de partida. Os
#   proximos modelos verificarao como as caracteristicas dos talhoes
#   alteram a probabilidade estimada de atingir o padrao de exportacao.
# COMO INTERPRETAR: se a proporcao for 0,40, 40 de cada 100 talhoes
#   observados atingiram o padrao. Neste momento, apresentamos apenas as
#   quantidades e porcentagens, sem classificar as classes como
#   balanceadas ou desbalanceadas. O possivel desbalanceamento sera
#   retomado no Step 8, principalmente durante a avaliacao da matriz de
#   confusao, da sensibilidade e da especificidade.
frequencia_exportacao <- table(resposta_exportacao)
proporcao_exportacao  <- mean(resposta_exportacao)

cat("\nresposta: padrao_exportacao\n")
cat(sprintf("  nao atingiu (0) : %3d  (%.1f%%)\n",
            frequencia_exportacao[["0"]], 100 * (1 - proporcao_exportacao)))
cat(sprintf("  atingiu     (1) : %3d  (%.1f%%)\n",
            frequencia_exportacao[["1"]], 100 * proporcao_exportacao))
cat(sprintf("  proporcao de 1  : %.4f\n", proporcao_exportacao))

# [APROFUNDAMENTO] Variancia de uma resposta 0/1.
# Em linguagem simples: quanto mais equilibrados os 0s e os 1s, maior a
# "incerteza" sobre o resultado de um talhao qualquer; quando quase todos
# sao 0 ou quase todos sao 1, essa incerteza e pequena.
# Em termos estatisticos: uma variavel de Bernoulli com probabilidade p
# tem variancia p(1-p). Ela e maxima em p = 0,5 (vale 0,25) e tende a
# zero quando p se aproxima de 0 ou de 1. A variancia fica, portanto,
# determinada pela media; nao existe um parametro de variancia separado.
# Essa relacao entre media e variancia e um dos motivos pelos quais o
# modelo linear normal nao e adequado como modelo principal para esta
# resposta. O Step 2 apresentara as demais limitacoes e justificara a
# escolha do modelo binomial.
cat(sprintf("  Var(Y) = p(1-p) : %.4f  (maximo possivel: 0,25, em p = 0,5)\n",
            proporcao_exportacao * (1 - proporcao_exportacao)))

# ---------------------------------------------------------------------
# Figura C0: distribuicao da resposta
# ---------------------------------------------------------------------
# Pergunta que responde: quantos lotes atingem o padrao?
# Como ler: a altura de cada barra e a contagem; o rotulo e a porcentagem.
# Limitacao: o grafico descreve a resposta, mas nao informa nada sobre
# as caracteristicas associadas a ela. Ele mostra o que sera modelado.
df_resposta <- data.frame(
    classe     = factor(c("Nao atingiu (0)", "Atingiu (1)"),
                        levels = c("Nao atingiu (0)", "Atingiu (1)")),
    n          = as.vector(frequencia_exportacao),
    porcentagem = 100 * as.vector(frequencia_exportacao) / nrow(dados_cafe))

g_resposta <- ggplot(df_resposta, aes(x = classe, y = n, fill = classe)) +
    geom_col(width = .6, colour = "black") +
    geom_text(aes(label = sprintf("%d  (%.1f%%)", n, porcentagem)),
              vjust = -.5, size = 4) +
    scale_fill_manual(values = c(VERMELHO, VERDE), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, .12))) +
    labs(x = NULL, y = "numero de talhoes",
         title = "Quantos lotes atingem o padrao de exportacao?",
         subtitle = sprintf("%d talhoes; proporcao observada de 1 = %.3f",
                            nrow(dados_cafe), proporcao_exportacao)) +
    theme_bw(base_size = 12)
salvar(g_resposta, "C0_distribuicao_resposta.png", 7, 5)

# ---------------------------------------------------------------------
# Niveis de referencia dos fatores
# ---------------------------------------------------------------------
# OBJETIVO: saber, antes de qualquer modelo, em relacao a qual nivel cada
#   coeficiente de fator sera lido.
# O QUE O CODIGO FAZ: imprime os niveis de cada fator. O primeiro nivel
#   e a categoria de referencia (ordem alfabetica, definida na leitura).
# POR QUE ISSO E NECESSARIO: nos proximos steps, cada coeficiente de um
#   fator representa uma comparacao entre um nivel e a referencia. Um
#   coeficiente para o nivel Icatu, por exemplo, so tem significado
#   quando se sabe que a comparacao e "Icatu em relacao a Bourbon".
#   Registrar as referencias aqui evita leituras erradas mais adiante.
#
# Referencias resultantes:
#   cultivar         -> Bourbon
#   manejo           -> Convencional
#   irrigacao        -> Nao
#   regiao_produtora -> Cerrado Mineiro
cat("\nniveis dos fatores (o primeiro e a referencia):\n")
for (f in fatores)
    cat(sprintf("  %-17s: %s\n", f, paste(levels(dados_cafe[[f]]), collapse = " | ")))

# ---------------------------------------------------------------------
# Formula das covariaveis candidatas
# ---------------------------------------------------------------------
# OBJETIVO: fixar, uma unica vez, o conjunto de covariaveis candidatas.
# O QUE O CODIGO FAZ: guarda o lado direito da formula (sem a resposta),
#   que sera reutilizado por todos os modelos dos proximos steps.
# POR QUE ISSO E NECESSARIO: definir o conjunto num unico lugar garante
#   que todos os modelos partam das mesmas candidatas e que nenhuma
#   covariavel seja adicionada sem registro. As exclusoes (broca, esforco
#   amostral, severidade) estao justificadas no cabecalho.
formula_candidatas <- ~ cultivar + manejo + irrigacao + regiao_produtora +
    altitude_m + declividade_pct + idade_lavoura_anos + densidade_plantio +
    adubacao_n_kg_ha + precipitacao_safra_mm + umidade_relativa_pct +
    ph_solo + materia_organica_pct

cat("\ncovariaveis candidatas:", length(attr(terms(formula_candidatas), "term.labels")), "\n")


# =====================================================================
# ---- STEP 0 ---- MODELO NULO
# =====================================================================
secao("STEP 0 - MODELO NULO")

# Modelo nulo, em uma frase: ele ignora as caracteristicas dos talhoes e
# atribui a mesma probabilidade de exportacao a todos eles. Essa
# probabilidade corresponde a proporcao geral observada na base. Ele
# sera usado como referencia para avaliar se as variaveis explicativas
# melhoram o modelo.
#
# OBJETIVO: ajustar esse modelo e mostrar que o unico numero que ele
#   estima (o intercepto) e a proporcao geral, escrita em outra escala.
# O QUE O CODIGO FAZ: glm() com a formula "resposta ~ 1" (so intercepto),
#   familia binomial (resposta 0/1) e funcao de ligacao logit.
# POR QUE ISSO E NECESSARIO: e a linha de base. Cada modelo com
#   covariaveis sera comparado com este para responder se as covariaveis
#   melhoram o ajuste.
#
# [APROFUNDAMENTO] Bernoulli e binomial agrupada.
# A familia binomial cobre dois formatos de dados. No primeiro, cada
# linha e um unico resultado 0/1 (uma observacao de Bernoulli), como
# aqui: um talhao, um lote, um resultado. No segundo, cada linha resume
# k sucessos em n tentativas (por exemplo, k lotes aprovados entre n
# lotes de um mesmo talhao), e o n entra no ajuste como peso. O glm()
# aceita os dois; neste trabalho todos os n valem 1.
modelo_nulo <- glm(padrao_exportacao ~ 1,
                   family = binomial(link = "logit"),
                   data   = dados_cafe)
verificar(modelo_nulo, "modelo_nulo")

# ---------------------------------------------------------------------
# O intercepto e a proporcao, em outra escala
# ---------------------------------------------------------------------
# Tres escalas, uma de cada vez:
#
# Probabilidade (p): de cada 100 talhoes, quantos atingem o padrao.
#   Varia entre 0 e 1. E a escala que interpretamos.
#
# Chance (em ingles, odds): quantos talhoes atingem o padrao para cada um
#   que nao atinge. Formula: p / (1 - p).
#   Exemplos: p = 0,5 da chance 1 (um para um); p = 0,8 da chance 4
#   (quatro para um). Varia de 0 a infinito.
#
# Logit: o logaritmo da chance. Formula: log(p / (1 - p)).
#   Pode assumir qualquer valor real, negativo ou positivo. E a escala em
#   que o modelo logistico combina os efeitos das covariaveis; so no fim
#   o resultado e convertido de volta para probabilidade.
#
# No R: qlogis(p) calcula o logit de p; plogis(x) devolve a probabilidade
# cujo logit e x.
#
# O QUE O CODIGO FAZ: compara o intercepto do modelo nulo com qlogis() da
#   proporcao observada e, em seguida, faz o caminho de volta com plogis().
# COMO INTERPRETAR: se os dois numeros coincidem, o intercepto do modelo
#   nulo e a proporcao observada escrita na escala do logit. Dentro do
#   modelo nulo, que nao utiliza nenhuma caracteristica dos talhoes, essa
#   e a estimativa da probabilidade comum atribuida a todos eles.
intercepto_nulo  <- unname(coef(modelo_nulo))
logit_proporcao  <- qlogis(proporcao_exportacao)

cat(sprintf("   proporcao observada          : %.6f\n", proporcao_exportacao))
cat(sprintf("   chance observada  p/(1-p)    : %.6f\n",
            proporcao_exportacao / (1 - proporcao_exportacao)))
cat(sprintf("   logit observado   qlogis(p)  : %.6f\n", logit_proporcao))
cat(sprintf("   intercepto do modelo nulo    : %.6f\n", intercepto_nulo))
cat("   coincidem (all.equal)        :",
    isTRUE(all.equal(intercepto_nulo, logit_proporcao, tolerance = 1e-6)), "\n")
cat(sprintf("   plogis(intercepto)           : %.4f  <- de volta a proporcao\n",
            plogis(intercepto_nulo)))

# [APROFUNDAMENTO] Por que all.equal() e nao ==.
# O glm() e ajustado por um algoritmo iterativo, chamado minimos quadrados
# reponderados iterativamente (IRLS, na sigla em ingles; slides, secao 6).
# O algoritmo para quando a melhoria entre iteracoes fica abaixo de uma
# tolerancia. Por isso os dois numeros coincidem ate a sexta ou oitava
# casa decimal, e nao bit a bit. A comparacao com == poderia dar FALSE
# por uma diferenca da ordem de 1e-10; all.equal() tolera essa diferenca.

# [APROFUNDAMENTO] Maxima verossimilhanca e ligacao canonica.
# O glm() escolhe o intercepto que torna os dados observados os mais
# provaveis possiveis; esse criterio e chamado de maxima verossimilhanca.
# A ligacao logit e a ligacao canonica da familia binomial (slides,
# secoes 4.1 e 5). Uma consequencia da ligacao canonica e que a media
# ajustada reproduz a media observada. E por isso que o intercepto do
# modelo nulo coincide com o logit da proporcao, e nao apenas se aproxima
# dele.

# ---------------------------------------------------------------------
# Figura C0: probabilidade, chance e logit
# ---------------------------------------------------------------------
# Pergunta que responde: como probabilidade, chance e logit se relacionam?
# Como ler o painel da esquerda: a chance cresce devagar para p pequeno e
# cresce muito rapido quando p se aproxima de 1.
# Como ler o painel da direita: o logit vale zero quando p = 0,5.
# Probabilidades igualmente afastadas de 0,5, como 0,2 e 0,8, produzem
# logits com o mesmo tamanho e sinais opostos. O logit pode assumir
# qualquer valor real.
# Em ambos os paineis, o ponto marca a proporcao observada na base.
# A probabilidade e a escala interpretada; o logit e a escala em que o
# modelo combina os efeitos.
p_grade <- seq(0.02, 0.98, by = 0.005)
df_escalas <- data.frame(p = p_grade,
                         chance = p_grade / (1 - p_grade),
                         logit  = qlogis(p_grade))
ponto <- data.frame(p = proporcao_exportacao,
                    chance = proporcao_exportacao / (1 - proporcao_exportacao),
                    logit  = logit_proporcao)

g_chance <- ggplot(df_escalas, aes(x = p, y = chance)) +
    geom_hline(yintercept = 1, colour = CINZA, linetype = "dashed") +
    geom_line(colour = VERDE, linewidth = 1) +
    geom_point(data = ponto, size = 3.5, colour = VERMELHO) +
    annotate("text", x = ponto$p, y = ponto$chance, hjust = 1.12, vjust = -.5,
             size = 3.3, colour = VERMELHO,
             label = sprintf("p = %.3f\nchance = %.3f", ponto$p, ponto$chance)) +
    coord_cartesian(ylim = c(0, 12)) +
    labs(x = "probabilidade  p", y = "chance  p / (1 - p)",
         title = "Chance: quantos exportam para cada um que nao exporta",
         subtitle = "linha tracejada: chance = 1 (p = 0,5)") +
    theme_bw(base_size = 11)

g_logit <- ggplot(df_escalas, aes(x = p, y = logit)) +
    geom_hline(yintercept = 0, colour = CINZA, linetype = "dashed") +
    geom_line(colour = VERDE, linewidth = 1) +
    geom_point(data = ponto, size = 3.5, colour = VERMELHO) +
    annotate("text", x = ponto$p, y = ponto$logit, hjust = -.12, vjust = 1.25,
             size = 3.3, colour = VERMELHO,
             label = sprintf("p = %.3f\nlogit = %.3f\n= intercepto do nulo",
                             ponto$p, ponto$logit)) +
    labs(x = "probabilidade  p", y = "logit  log(p / (1 - p))",
         title = "Logit: a escala em que o modelo faz contas",
         subtitle = "simetrico em torno de p = 0,5; cobre toda a reta real") +
    theme_bw(base_size = 11)

salvar(g_chance | g_logit, "C0_prob_chance_logit.png", 11, 5)

# ---------------------------------------------------------------------
# [APROFUNDAMENTO] Deviance nula
# ---------------------------------------------------------------------
# Em linguagem simples: a deviance mede a distancia entre o modelo
# ajustado e um modelo de referencia que reproduz os dados perfeitamente.
#
# Em termos estatisticos (slides, secao 7.1): esse modelo de referencia
# e o modelo saturado, que tem um parametro por observacao. A deviance e
# a diferenca de log-verossimilhanca entre o modelo saturado e o modelo
# ajustado, multiplicada por 2. Entre modelos ajustados a mesma resposta
# e as mesmas observacoes, uma deviance menor indica maior proximidade
# do modelo saturado. A comparacao so faz sentido nessas condicoes:
# deviances de modelos com respostas ou bases diferentes nao sao
# comparaveis. No modelo linear normal, a deviance coincide com a soma
# de quadrados dos residuos.
#
# A deviance nula e a deviance do modelo so com intercepto. Seus graus de
# liberdade sao n - 1, porque um parametro (o intercepto) foi estimado.
# No Step 3, a queda da deviance ao incluir covariaveis sera usada para
# testar se elas melhoram o ajuste (analise de deviance).
#
# Sobre a segunda linha impressa abaixo: neste caso de observacoes
# individuais 0/1, a deviance coincide com -2 vezes a log-verossimilhanca
# porque a log-verossimilhanca do modelo saturado e zero. Essa igualdade
# nao vale para qualquer modelo; e uma particularidade dos dados binarios
# nao agrupados.
deviance_nula <- deviance(modelo_nulo)
gl_nulo       <- df.residual(modelo_nulo)
cat(sprintf("\n   deviance nula : %.4f  com %d graus de liberdade (n - 1)\n",
            deviance_nula, gl_nulo))
cat(sprintf("   -2 * logLik   : %.4f  (deve ser identico: deviance = -2 log-verossimilhanca)\n",
            -2 * as.numeric(logLik(modelo_nulo))))

# ---------------------------------------------------------------------
# Parecer do Step 0
# ---------------------------------------------------------------------
# O parecer abaixo esta organizado em partes, para separar o que e
# descricao dos dados do que e leitura estatistica:
#
# CONCLUSAO LITERAL:
# Resume o resultado em quantidades e porcentagens.
#
# CONCLUSAO ESTATISTICA:
# Relaciona a resposta Bernoulli, o modelo nulo e a escala logit.
#
# REFERENCIA PARA OS PROXIMOS STEPS:
# Registra a deviance nula que sera usada nas comparacoes futuras.
#
# PROXIMO PASSO:
# Investigar quais caracteristicas alteram a probabilidade estimada.
#
# Observacao: a expressao "melhor estimativa", no texto do parecer, deve
# ser lida dentro da estrutura do modelo nulo, que nao utiliza nenhuma
# caracteristica dos talhoes.
parecer("CONCLUSAO LITERAL. De cada ", nrow(dados_cafe), " talhoes da cooperativa, ",
        frequencia_exportacao[["1"]], " produziram lote com padrao de exportacao e ",
        frequencia_exportacao[["0"]], " nao (", sprintf("%.1f%%", 100 * proporcao_exportacao),
        " contra ", sprintf("%.1f%%", 100 * (1 - proporcao_exportacao)),
        "). Sem saber nada sobre um talhao, a melhor estimativa da probabilidade ",
        "de o lote dele exportar e essa proporcao - e e exatamente isso que o ",
        "modelo nulo diz para todos os talhoes. ",
        "CONCLUSAO ESTATISTICA. A resposta e binaria: cada talhao contribui com uma ",
        "observacao de Bernoulli, com variancia p(1-p) = ",
        sprintf("%.3f", proporcao_exportacao * (1 - proporcao_exportacao)),
        " determinada pela propria media. O modelo nulo com ligacao logit devolve ",
        "intercepto ", sprintf("%.4f", intercepto_nulo), ", identico a log(p/(1-p)) da ",
        "proporcao observada - consequencia da ligacao canonica, pela qual a media ",
        "ajustada reproduz a media amostral. A deviance nula de ",
        sprintf("%.2f", deviance_nula), " com ", gl_nulo, " graus de liberdade e a ",
        "referencia contra a qual todo modelo com covariaveis sera testado. ",
        "O QUE FICA PARA OS PROXIMOS STEPS: quais caracteristicas dos talhoes movem ",
        "essa probabilidade para cima ou para baixo (Steps 1 a 4) e se vale a pena ",
        "modelar (Step 3).")


# =====================================================================
# ---- STEP 1 ---- EXPLORACAO CONDICIONAL: FATORES
# =====================================================================
secao("STEP 1 - EXPLORACAO CONDICIONAL: FATORES")

# Pergunta do step: olhando uma caracteristica de cada vez, quais parecem
# mudar a proporcao de lotes com padrao de exportacao? E, para as
# caracteristicas numericas, essa mudanca parece "em linha reta" na
# escala em que o modelo trabalha (o logit)?
#
# Aviso que vale para todo o step: cada tabela e cada grafico daqui olha
# UMA caracteristica de cada vez, sem levar em conta as demais. Isso se
# chama analise marginal. Uma diferenca vista aqui pode encolher, crescer
# ou mudar de sinal quando as outras caracteristicas entrarem no modelo
# (Step 3). Por isso tudo neste step e hipotese, nao conclusao.

# ---------------------------------------------------------------------
# Proporcao condicional por nivel de cada fator
# ---------------------------------------------------------------------
# OBJETIVO: ver a proporcao de lotes com padrao de exportacao dentro de
#   cada grupo (so os talhoes Bourbon, so os organicos, e assim por diante).
# O QUE O CODIGO FAZ: para cada fator, agrupa os talhoes por nivel e
#   calcula, em cada nivel, o numero de talhoes (n), quantos atingiram o
#   padrao e a proporcao. Calcula tambem a chance, p / (1 - p), ja
#   apresentada no Step 0.
# POR QUE ISSO E NECESSARIO: a proporcao dentro de um grupo e chamada de
#   proporcao condicional. Comparar as proporcoes condicionais entre os
#   niveis e a forma mais simples de perceber se um fator parece estar
#   associado a resposta.
# COMO INTERPRETAR: um nivel com proporcao acima da geral parece
#   favorecer o padrao de exportacao; abaixo, parece desfavorecer. A
#   coluna n importa: uma proporcao alta num grupo pequeno e menos
#   confiavel do que num grupo grande.
#
# [APROFUNDAMENTO] Intervalo de confianca de Wilson para uma proporcao.
# A proporcao de cada grupo e uma estimativa e carrega incerteza. O
# intervalo de Wilson e uma forma de expressar essa incerteza que se
# comporta bem mesmo em grupos pequenos ou com proporcoes perto de 0 ou
# de 1. No R, prop.test() sem correcao de continuidade devolve esse
# intervalo. Ele aparece como barra de erro na figura C1.
tabela_por_fator <- function(f) {
    tab <- dados_cafe %>%
        group_by(nivel = .data[[f]]) %>%
        summarise(n         = n(),
                  atingiu   = sum(padrao_exportacao),
                  proporcao = mean(padrao_exportacao),
                  .groups   = "drop") %>%
        mutate(fator  = f,
               chance = proporcao / (1 - proporcao))
    ic <- mapply(function(x, n) prop.test(x, n, correct = FALSE)$conf.int,
                 tab$atingiu, tab$n)
    tab$ic_inf <- ic[1, ]
    tab$ic_sup <- ic[2, ]
    tab
}

tab_fatores <- bind_rows(lapply(fatores, tabela_por_fator))

cat(sprintf("proporcao geral de referencia: %.3f\n", proporcao_exportacao))
for (f in fatores) {
    cat("\n", f, "\n", sep = "")
    t_f <- tab_fatores %>% filter(fator == f)
    print(data.frame(nivel     = as.character(t_f$nivel),
                     n         = t_f$n,
                     atingiu   = t_f$atingiu,
                     proporcao = round(t_f$proporcao, 3),
                     chance    = round(t_f$chance, 3),
                     IC95_inf  = round(t_f$ic_inf, 3),
                     IC95_sup  = round(t_f$ic_sup, 3)),
          row.names = FALSE)
}

# Verificacao: grupos pequenos. Abaixo de ~15 observacoes, a proporcao
# de um nivel e pouco confiavel. O limite e uma regra pratica deste
# projeto para emitir um alerta, nao um criterio estatistico formal.
MINIMO_POR_NIVEL <- 15
niveis_pequenos <- tab_fatores %>% filter(n < MINIMO_POR_NIVEL)
if (nrow(niveis_pequenos) > 0) {
    cat("\n   [ATENCAO] niveis com menos de ", MINIMO_POR_NIVEL, " observacoes:\n", sep = "")
    print(niveis_pequenos[, c("fator", "nivel", "n")], row.names = FALSE)
} else {
    cat(sprintf("\n   [ok] todos os niveis tem pelo menos %d observacoes (menor: %d)\n",
                MINIMO_POR_NIVEL, min(tab_fatores$n)))
}

# Amplitude de cada fator: diferenca entre a maior e a menor proporcao
# entre os niveis. E uma medida informal de "quanto o fator parece
# importar" na analise marginal. Sera usada no parecer.
amplitude_fatores <- tab_fatores %>%
    group_by(fator) %>%
    summarise(nivel_max = as.character(nivel[which.max(proporcao)]),
              prop_max  = max(proporcao),
              nivel_min = as.character(nivel[which.min(proporcao)]),
              prop_min  = min(proporcao),
              amplitude = prop_max - prop_min,
              .groups   = "drop") %>%
    arrange(desc(amplitude))

cat("\namplitude da proporcao entre niveis (maior - menor), por fator:\n")
print(data.frame(fator     = amplitude_fatores$fator,
                 maior     = sprintf("%s (%.3f)", amplitude_fatores$nivel_max, amplitude_fatores$prop_max),
                 menor     = sprintf("%s (%.3f)", amplitude_fatores$nivel_min, amplitude_fatores$prop_min),
                 amplitude = round(amplitude_fatores$amplitude, 3)),
      row.names = FALSE)

# ---------------------------------------------------------------------
# Figura C1: proporcao por nivel de cada fator
# ---------------------------------------------------------------------
# Pergunta que responde: qual nivel de cada fator tem maior proporcao de
# lotes com padrao de exportacao?
# Como ler: a altura da barra e a proporcao do nivel; a barra de erro e o
# intervalo de confianca de Wilson (95%); o n esta escrito na base da
# barra; a linha tracejada e a proporcao geral da base. Barras cujo
# intervalo cruza a linha tracejada nao se distinguem claramente da media.
# Limitacao: analise marginal. As diferencas podem mudar quando os fatores
# forem considerados em conjunto.
tab_fatores_plot <- tab_fatores %>%
    mutate(nivel = factor(nivel, levels = unique(nivel)),
           fator = factor(fator, levels = fatores))

g_fatores <- ggplot(tab_fatores_plot, aes(x = nivel, y = proporcao)) +
    geom_hline(yintercept = proporcao_exportacao, linetype = "dashed", colour = VERMELHO) +
    geom_col(fill = VERDE, colour = "black", width = .65) +
    geom_errorbar(aes(ymin = ic_inf, ymax = ic_sup), width = .2, colour = "black") +
    geom_text(aes(y = 0.03, label = sprintf("n = %d", n)), colour = "white", size = 3.2) +
    facet_wrap(~ fator, scales = "free_x") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
    labs(x = NULL, y = "proporcao com padrao de exportacao",
         title = "Proporcao de lotes com padrao de exportacao, por nivel de cada fator",
         subtitle = sprintf("barras de erro: IC de Wilson 95%%; linha tracejada: proporcao geral = %.3f",
                            proporcao_exportacao),
         caption = "analise marginal: uma caracteristica de cada vez, sem ajuste pelas demais") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 15, hjust = 1))
salvar(g_fatores, "C1_proporcao_por_fator.png", 11, 8)


# =====================================================================
# ---- STEP 1 ---- EXPLORACAO CONDICIONAL: NUMERICAS
# =====================================================================
secao("STEP 1 - EXPLORACAO CONDICIONAL: NUMERICAS")

# OBJETIVO: ver se as caracteristicas numericas (altitude, idade da
#   lavoura, ...) diferem entre os talhoes que atingiram o padrao e os
#   que nao atingiram.
# O QUE O CODIGO FAZ: identifica as candidatas numericas a partir da
#   formula do Step 0 (tudo o que nao e fator), reorganiza a base no
#   formato longo (uma linha por talhao e por variavel) e calcula a
#   mediana de cada variavel dentro de cada classe da resposta.
# POR QUE ISSO E NECESSARIO: e a versao, para variaveis numericas, da
#   pergunta feita aos fatores. Se a mediana de altitude dos talhoes que
#   atingiram o padrao e maior do que a dos que nao atingiram, a altitude
#   parece estar associada a resposta.
# COMO INTERPRETAR: diferencas de mediana grandes em relacao a dispersao
#   da variavel sugerem associacao; medianas quase iguais sugerem pouca
#   associacao. Nada aqui diz qual e a causa e qual e o efeito.
numericas <- setdiff(attr(terms(formula_candidatas), "term.labels"), fatores)
cat("candidatas numericas:", length(numericas), "\n")

dados_longos <- dados_cafe %>%
    select(all_of(numericas), padrao_exportacao) %>%
    pivot_longer(all_of(numericas), names_to = "variavel", values_to = "valor") %>%
    mutate(classe = factor(ifelse(padrao_exportacao == 1, "Atingiu (1)", "Nao atingiu (0)"),
                           levels = c("Nao atingiu (0)", "Atingiu (1)")),
           variavel = factor(variavel, levels = numericas))

tab_medianas <- dados_longos %>%
    group_by(variavel, classe) %>%
    summarise(mediana = median(valor), .groups = "drop") %>%
    pivot_wider(names_from = classe, values_from = mediana) %>%
    mutate(diferenca = `Atingiu (1)` - `Nao atingiu (0)`)

cat("\nmediana de cada variavel numerica, por classe da resposta:\n")
print(data.frame(variavel          = as.character(tab_medianas$variavel),
                 mediana_nao       = round(tab_medianas$`Nao atingiu (0)`, 2),
                 mediana_atingiu   = round(tab_medianas$`Atingiu (1)`, 2),
                 diferenca         = round(tab_medianas$diferenca, 2)),
      row.names = FALSE)

# ---------------------------------------------------------------------
# Figura C1: distribuicao das numericas por classe da resposta
# ---------------------------------------------------------------------
# Pergunta que responde: as caracteristicas numericas diferem entre quem
# atingiu e quem nao atingiu o padrao?
# Como ler: cada painel e uma variavel; cada caixa resume a distribuicao
# em uma classe (linha central = mediana; caixa = metade central dos
# talhoes). Caixas deslocadas uma da outra sugerem associacao; caixas
# sobrepostas sugerem pouca associacao.
# Limitacao: este e um painel de exploracao. So as variaveis com
# diferenca visivel merecem ir para o relatorio. O grafico tambem nao
# mostra a forma da relacao na escala do modelo; isso e feito a seguir.
g_numericas <- ggplot(dados_longos, aes(x = classe, y = valor, fill = classe)) +
    geom_boxplot(width = .6, outlier.size = 1, outlier.alpha = .5) +
    scale_fill_manual(values = c(VERMELHO, VERDE), guide = "none") +
    facet_wrap(~ variavel, scales = "free_y", ncol = 3) +
    labs(x = NULL, y = NULL,
         title = "Distribuicao de cada caracteristica numerica, por classe da resposta",
         subtitle = "caixas deslocadas sugerem associacao; caixas sobrepostas sugerem pouca associacao",
         caption = "painel de exploracao; analise marginal") +
    theme_bw(base_size = 11)
salvar(g_numericas, "C1_numericas_por_resposta.png", 11, 9)


# =====================================================================
# ---- STEP 1 ---- LOGIT EMPIRICO: A RELACAO PARECE LINEAR NA ESCALA DO MODELO?
# =====================================================================
secao("STEP 1 - LOGIT EMPIRICO")

# Por que olhar o logit, e nao a probabilidade:
# O modelo logistico e linear na escala do logit (Step 0). Isso quer
# dizer que, se a altitude entrar no modelo "como esta", o modelo assume
# que cada metro a mais soma a mesma quantidade ao logit, em qualquer
# ponto da faixa de altitude. Uma relacao que parece reta na escala da
# probabilidade seria curva no logit, e vice-versa. Por isso a pergunta
# "parece linear?" tem de ser feita na escala do logit.
#
# Logit empirico, em linguagem simples: dividir os talhoes em faixas de
# uma variavel (por exemplo, cinco faixas de altitude com o mesmo numero
# de talhoes), calcular a proporcao de sucesso em cada faixa e converter
# essa proporcao em logit. Cinco pontos quase alinhados sugerem que o
# termo linear e adequado; uma curva sugere que a variavel pode precisar
# de transformacao. O nome tecnico e "logit empirico" (Agresti, Categorical
# Data Analysis, cap. 4); e um metodo informal, de exploracao.
#
# OBJETIVO: ver, para cada numerica, se a relacao com a resposta parece
#   aproximadamente linear na escala do logit.
# O QUE O CODIGO FAZ: a funcao logit_empirico() corta a variavel em k
#   faixas de tamanho parecido (quantis), calcula em cada faixa o n, o
#   numero de sucessos, a proporcao, o valor medio da variavel e o logit.
# POR QUE ISSO E NECESSARIO: a linearidade no logit e uma suposicao do
#   modelo. Ver os dados agora evita descobrir tarde (no diagnostico do
#   Step 5) que uma variavel precisava de outra forma.
# COMO INTERPRETAR: pontos alinhados = termo linear plausivel; pontos em
#   curva = hipotese de nao linearidade, a ser testada formalmente no
#   Step 5. Pontos muito espalhados, com barras de erro largas, indicam
#   que os dados nao permitem distinguir uma tendencia do ruido.
#
# [APROFUNDAMENTO] Correcao de continuidade e erro-padrao do logit.
# Se uma faixa tiver proporcao 0 ou 1, o logit e -infinito ou +infinito e
# nao pode ser desenhado. A correcao usada aqui soma 0,5 aos sucessos e
# aos fracassos: logit = log((s + 0,5) / (n - s + 0,5)), com s = sucessos
# e n = talhoes na faixa. Quando a proporcao esta longe de 0 e de 1, a
# correcao quase nao muda o valor. O erro-padrao aproximado desse logit e
# sqrt(1/(s + 0,5) + 1/(n - s + 0,5)); ele define as barras de erro da
# figura (logit +- 1,96 erro-padrao).
logit_empirico <- function(x, y, k = 5) {
    cortes <- unique(quantile(x, probs = seq(0, 1, length.out = k + 1)))
    faixa  <- cut(x, breaks = cortes, include.lowest = TRUE)
    data.frame(faixa    = levels(faixa),
               n        = as.vector(table(faixa)),
               sucessos = as.vector(tapply(y, faixa, sum)),
               x_medio  = as.vector(tapply(x, faixa, mean))) %>%
        mutate(proporcao       = sucessos / n,
               logit_bruto     = qlogis(proporcao),
               logit_corrigido = log((sucessos + 0.5) / (n - sucessos + 0.5)),
               ep_logit        = sqrt(1 / (sucessos + 0.5) + 1 / (n - sucessos + 0.5)))
}

K_FAIXAS <- 5
tab_logit_empirico <- bind_rows(lapply(numericas, function(v) {
    cbind(variavel = v, logit_empirico(dados_cafe[[v]], resposta_exportacao, k = K_FAIXAS))
})) %>% mutate(variavel = factor(variavel, levels = numericas))

for (v in numericas) {
    cat("\n", v, "\n", sep = "")
    t_v <- tab_logit_empirico %>% filter(variavel == v)
    print(data.frame(faixa           = t_v$faixa,
                     n               = t_v$n,
                     sucessos        = t_v$sucessos,
                     proporcao       = round(t_v$proporcao, 3),
                     logit_bruto     = round(t_v$logit_bruto, 3),
                     logit_corrigido = round(t_v$logit_corrigido, 3)),
          row.names = FALSE)
}

# Verificacoes: numero de faixas efetivo (quantis repetidos reduzem o
# numero de faixas) e faixas em que a correcao de continuidade foi
# necessaria (proporcao 0 ou 1).
faixas_por_variavel <- tab_logit_empirico %>% count(variavel, name = "faixas")
if (any(faixas_por_variavel$faixas < K_FAIXAS)) {
    cat("\n   [ATENCAO] variaveis com menos de ", K_FAIXAS, " faixas (quantis repetidos):\n", sep = "")
    print(faixas_por_variavel %>% filter(faixas < K_FAIXAS), row.names = FALSE)
} else {
    cat(sprintf("\n   [ok] todas as %d variaveis foram divididas em %d faixas\n",
                length(numericas), K_FAIXAS))
}
faixas_extremas <- sum(tab_logit_empirico$proporcao %in% c(0, 1))
cat(sprintf("   [ok] faixas com proporcao exatamente 0 ou 1 (correcao necessaria): %d\n",
            faixas_extremas))
cat(sprintf("   [ok] menor faixa: %d talhoes; erro-padrao tipico do logit: %.2f\n",
            min(tab_logit_empirico$n), median(tab_logit_empirico$ep_logit)))

# Resumo por variavel: amplitude do logit entre a primeira e a ultima
# faixa (sinal e tamanho da tendencia) e o R2 de uma reta passada pelos
# k pontos (quanto da variacao entre faixas uma reta explica). Sao
# medidas informais para ordenar as hipoteses; nao sao testes.
resumo_logit <- tab_logit_empirico %>%
    group_by(variavel) %>%
    summarise(logit_faixa1   = first(logit_corrigido),
              logit_faixaK   = last(logit_corrigido),
              amplitude      = logit_faixaK - logit_faixa1,
              r2_reta        = cor(x_medio, logit_corrigido)^2,
              .groups        = "drop") %>%
    arrange(desc(abs(amplitude)))

cat("\nresumo por variavel (ordenado pela amplitude absoluta do logit):\n")
print(data.frame(variavel     = as.character(resumo_logit$variavel),
                 logit_faixa1 = round(resumo_logit$logit_faixa1, 3),
                 logit_faixaK = round(resumo_logit$logit_faixaK, 3),
                 amplitude    = round(resumo_logit$amplitude, 3),
                 r2_reta      = round(resumo_logit$r2_reta, 2)),
      row.names = FALSE)

# ---------------------------------------------------------------------
# Figura C1: logit empirico por faixa
# ---------------------------------------------------------------------
# Pergunta que responde: a relacao de cada numerica com a resposta parece
# linear na escala do logit?
# Como ler: cada ponto e uma faixa (eixo x = valor medio da variavel na
# faixa; eixo y = logit da proporcao de sucesso); a barra de erro e
# logit +- 1,96 erro-padrao; a linha pontilhada vermelha e uma reta
# ajustada aos pontos, so como referencia visual; a linha tracejada cinza
# e o logit da proporcao geral. O eixo y e o mesmo em todos os paineis,
# para que as inclinacoes sejam comparaveis.
# Limitacao: poucos pontos por variavel; as faixas sao arbitrarias; e um
# metodo informal. A verificacao formal da escala fica para o Step 5.
g_logit_emp <- ggplot(tab_logit_empirico, aes(x = x_medio, y = logit_corrigido)) +
    geom_hline(yintercept = logit_proporcao, linetype = "dashed", colour = CINZA) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                colour = VERMELHO, linewidth = .6, linetype = "dotted") +
    geom_errorbar(aes(ymin = logit_corrigido - 1.96 * ep_logit,
                      ymax = logit_corrigido + 1.96 * ep_logit),
                  width = 0, colour = CINZA) +
    geom_line(colour = VERDE, linewidth = .5) +
    geom_point(colour = VERDE, size = 2.5) +
    facet_wrap(~ variavel, scales = "free_x", ncol = 3) +
    labs(x = "valor medio da variavel na faixa", y = "logit empirico da proporcao de sucesso",
         title = "Logit empirico por faixa: a relacao parece linear na escala do modelo?",
         subtitle = sprintf("%d faixas de tamanho parecido por variavel; barras: logit +- 1,96 EP; pontilhado: reta de referencia",
                            K_FAIXAS),
         caption = "metodo informal de exploracao; a verificacao formal da escala e feita no Step 5") +
    theme_bw(base_size = 11)
salvar(g_logit_emp, "C1_logit_empirico.png", 11, 9)

# ---------------------------------------------------------------------
# Parecer do Step 1
# ---------------------------------------------------------------------
# O parecer lista hipoteses, nao conclusoes. Os criterios abaixo sao
# regras praticas declaradas para ordenar o que parece importar:
#   - fatores: ordenados pela amplitude da proporcao entre niveis;
#   - numericas: ordenadas pela amplitude absoluta do logit entre a
#     primeira e a ultima faixa; uma amplitude menor do que o dobro do
#     erro-padrao tipico e tratada como "nao distinguivel do ruido";
#   - possivel nao linearidade: amplitude acima desse limite, mas reta
#     explicando menos da metade da variacao entre faixas (R2 < 0,5).
limite_ruido <- 2 * median(tab_logit_empirico$ep_logit)
num_tendencia <- resumo_logit %>% filter(abs(amplitude) >= limite_ruido)
num_ruido     <- resumo_logit %>% filter(abs(amplitude) <  limite_ruido)
num_curva     <- num_tendencia %>% filter(r2_reta < 0.5)

descreve_tendencia <- function(t) {
    if (nrow(t) == 0) return("nenhuma")
    paste(sprintf("%s (%s, amplitude %.2f)", t$variavel,
                  ifelse(t$amplitude > 0, "crescente", "decrescente"), t$amplitude),
          collapse = "; ")
}

parecer("CONCLUSAO LITERAL. Entre os fatores, a maior diferenca de proporcao ",
        "entre niveis aparece em ", amplitude_fatores$fator[1], ": de ",
        sprintf("%.1f%%", 100 * amplitude_fatores$prop_min[1]), " em ",
        amplitude_fatores$nivel_min[1], " a ",
        sprintf("%.1f%%", 100 * amplitude_fatores$prop_max[1]), " em ",
        amplitude_fatores$nivel_max[1], ". A menor aparece em ",
        amplitude_fatores$fator[nrow(amplitude_fatores)], " (",
        sprintf("%.1f", 100 * amplitude_fatores$amplitude[nrow(amplitude_fatores)]),
        " pontos percentuais entre o maior e o menor nivel). Entre as numericas, ",
        "as que mostram tendencia acima do ruido na escala do logit sao: ",
        descreve_tendencia(num_tendencia), ". As que nao se distinguem do ruido ",
        "com ", K_FAIXAS, " faixas sao: ",
        if (nrow(num_ruido) > 0) paste(num_ruido$variavel, collapse = ", ") else "nenhuma", ". ",
        "CONCLUSAO ESTATISTICA. Todas as diferencas acima sao marginais: cada ",
        "caracteristica foi olhada isoladamente, sem ajuste pelas demais. Para as ",
        "numericas com tendencia, a reta ajustada aos ", K_FAIXAS, " pontos explica ",
        "de ", sprintf("%.0f%%", 100 * min(num_tendencia$r2_reta)), " a ",
        sprintf("%.0f%%", 100 * max(num_tendencia$r2_reta)), " da variacao entre ",
        "faixas, o que e compativel com o termo linear no logit para ",
        if (nrow(num_curva) == 0) "todas elas" else
            paste0("a maioria; ", paste(num_curva$variavel, collapse = ", "),
                   " tem sinal de possivel nao linearidade a testar"), ". ",
        "RESPOSTA PARCIAL. Ha caracteristicas candidatas fortes e fracas, mas ",
        "nenhuma delas pode ainda ser chamada de determinante: a associacao marginal ",
        "nao diz se o efeito se mantem quando as outras caracteristicas sao ",
        "consideradas, nem qual e a direcao causal. ",
        "O QUE FICA PARA OS PROXIMOS STEPS: os efeitos condicionais, com todas as ",
        "candidatas juntas (Step 3); a verificacao formal da escala de cada numerica ",
        "no preditor linear (Step 5).")


# =====================================================================
# ---- STEP 2 ---- POR QUE NAO A REGRESSAO LINEAR
# =====================================================================
secao("STEP 2 - POR QUE NAO A REGRESSAO LINEAR")

# Pergunta do step: por que nao usar a regressao linear que ja conhecemos?
# E, entre as funcoes de ligacao possiveis para uma resposta 0/1, por que
# a logit?
#
# A regressao linear (modelo linear normal) e o ponto de partida natural,
# porque e o modelo mais conhecido. Aqui ela e ajustada de proposito, com
# todas as candidatas, para que suas limitacoes aparecam nos proprios
# dados e nao apenas na teoria. Os slides da disciplina (secao 2.1)
# chamam isso de "os fracassos do modelo linear" para uma resposta
# binaria; sao dois, mostrados a seguir.

# ---------------------------------------------------------------------
# Limitacao 1: predicoes fora do intervalo possivel
# ---------------------------------------------------------------------
# OBJETIVO: verificar se a regressao linear produz "probabilidades"
#   menores que 0 ou maiores que 1 para algum talhao.
# O QUE O CODIGO FAZ: ajusta lm() com a resposta 0/1 e todas as
#   candidatas (o mesmo conjunto do Step 0), extrai os valores ajustados
#   e conta quantos ficam fora de [0, 1].
# POR QUE ISSO E NECESSARIO: uma reta nao tem piso nem teto. A resposta
#   e uma probabilidade, que so existe entre 0 e 1. Se a reta sai desse
#   intervalo para algum talhao, o modelo esta prevendo algo impossivel.
# COMO INTERPRETAR: um unico valor fora de [0, 1] ja mostra o problema
#   estrutural; nao e a quantidade que importa, e o fato de o modelo nao
#   ter como impedir que isso ocorra.
formula_completa <- update(formula_candidatas, padrao_exportacao ~ .)

modelo_linear_amplo <- lm(formula_completa, data = dados_cafe)
pred_linear <- fitted(modelo_linear_amplo)
fora_dominio <- sum(pred_linear < 0 | pred_linear > 1)

cat(sprintf("regressao linear com as %d candidatas:\n", length(numericas) + length(fatores)))
cat(sprintf("   faixa das predicoes     : %.4f a %.4f\n", min(pred_linear), max(pred_linear)))
cat(sprintf("   predicoes fora de [0, 1]: %d de %d talhoes\n", fora_dominio, nrow(dados_cafe)))
if (fora_dominio > 0) {
    piores <- dados_cafe$id_talhao[order(pmax(pred_linear - 1, -pred_linear), decreasing = TRUE)][seq_len(min(3, fora_dominio))]
    cat("   talhoes com predicao impossivel (ate 3):", paste(piores, collapse = ", "), "\n")
}

# ---------------------------------------------------------------------
# Limitacao 2: a variancia nao e constante
# ---------------------------------------------------------------------
# OBJETIVO: mostrar que o espalhamento dos erros da regressao linear muda
#   conforme a probabilidade prevista, em vez de ser constante.
# O QUE O CODIGO FAZ: divide os talhoes em cinco faixas pelo valor
#   ajustado da regressao linear e calcula, em cada faixa, a variancia
#   dos residuos observados e a variancia que uma resposta 0/1 teria
#   naquela probabilidade, p(1-p) (Step 0).
# POR QUE ISSO E NECESSARIO: a regressao linear supoe que os erros tem a
#   mesma variancia em toda a faixa de predicao (homocedasticidade). Para
#   uma resposta 0/1 isso nao pode ser verdade: a variancia e p(1-p), que
#   depende da propria probabilidade. Nao e um defeito dos dados; e uma
#   propriedade da resposta.
# COMO INTERPRETAR: se a variancia observada dos residuos acompanha
#   p(1-p) de faixa em faixa, a suposicao de variancia constante esta
#   sendo violada de forma sistematica, e os erros-padrao da regressao
#   linear nao sao confiaveis.
faixa_pred <- cut(pred_linear, quantile(pred_linear, seq(0, 1, .2)), include.lowest = TRUE)
tab_variancia <- data.frame(
    faixa        = levels(faixa_pred),
    n            = as.vector(table(faixa_pred)),
    p_medio      = as.vector(tapply(pred_linear, faixa_pred, mean)),
    var_residuos = as.vector(tapply(residuals(modelo_linear_amplo), faixa_pred, var)))
tab_variancia$var_teorica <- tab_variancia$p_medio * (1 - tab_variancia$p_medio)

cat("\nvariancia dos residuos da regressao linear, por faixa de predicao:\n")
print(data.frame(faixa        = tab_variancia$faixa,
                 n            = tab_variancia$n,
                 p_medio      = round(tab_variancia$p_medio, 3),
                 var_residuos = round(tab_variancia$var_residuos, 4),
                 var_teorica  = round(tab_variancia$var_teorica, 4)),
      row.names = FALSE)
cat(sprintf("   razao entre a maior e a menor variancia observada: %.2f\n",
            max(tab_variancia$var_residuos) / min(tab_variancia$var_residuos)))

# ---------------------------------------------------------------------
# Figura C2: faixa das predicoes
# ---------------------------------------------------------------------
# Pergunta que responde: a regressao linear sai do intervalo possivel?
# Como ler: cada segmento vai da menor a maior predicao de um modelo; a
# faixa sombreada marca o que e impossivel para uma probabilidade
# (abaixo de 0 e acima de 1). O modelo logistico e ajustado aqui apenas
# para a comparacao; ele e o assunto dos proximos blocos.
# Limitacao: mostra o alcance das predicoes, nao a qualidade do ajuste.
modelo_logit_amplo <- glm(formula_completa, family = binomial(link = "logit"), data = dados_cafe)
pred_logit <- fitted(modelo_logit_amplo)

df_faixas <- data.frame(
    modelo = factor(c("observado (0 ou 1)", "regressao linear", "logistico (logit)"),
                    levels = c("logistico (logit)", "regressao linear", "observado (0 ou 1)")),
    minimo = c(0, min(pred_linear), min(pred_logit)),
    maximo = c(1, max(pred_linear), max(pred_logit)))
df_faixas$impossivel <- df_faixas$minimo < 0 | df_faixas$maximo > 1

g_faixas <- ggplot(df_faixas, aes(y = modelo)) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = VERMELHO, alpha = .07) +
    annotate("rect", xmin = 1, xmax = Inf, ymin = -Inf, ymax = Inf, fill = VERMELHO, alpha = .07) +
    geom_vline(xintercept = c(0, 1), colour = VERMELHO, linetype = "dashed") +
    geom_segment(aes(x = minimo, xend = maximo, yend = modelo, colour = impossivel),
                 linewidth = 5, lineend = "round") +
    scale_colour_manual(values = c(`FALSE` = VERDE, `TRUE` = VERMELHO), guide = "none") +
    annotate("text", x = 1, y = Inf, label = "acima de 1: impossivel", hjust = -.05,
             vjust = 1.8, colour = VERMELHO, size = 3.3) +
    coord_cartesian(xlim = c(-0.06, 1.34)) +
    labs(x = "probabilidade predita", y = NULL,
         title = "Faixa das predicoes: a regressao linear sai do intervalo possivel",
         subtitle = sprintf("regressao linear: %d predicao(oes) fora de [0, 1]; a maior vale %.3f",
                            fora_dominio, max(pred_linear)),
         caption = "mesmas 13 candidatas nos dois modelos") +
    theme_bw(base_size = 12)
salvar(g_faixas, "C2_faixa_predicoes.png", 9, 4.5)


# =====================================================================
# ---- STEP 2 ---- O MODELO LINEAR GENERALIZADO E A FUNCAO DE LIGACAO
# =====================================================================
secao("STEP 2 - QUAL FUNCAO DE LIGACAO?")

# As duas limitacoes acima tem a mesma origem: a regressao linear liga a
# probabilidade diretamente a uma reta. O modelo linear generalizado
# (MLG) resolve isso separando o problema em tres pecas (slides, secao
# 4.1). Uma de cada vez:
#
# 1. Componente aleatorio: qual distribuicao a resposta segue. Aqui,
#    cada talhao e um resultado 0/1 - familia binomial. Isso ja embute a
#    variancia p(1-p) no modelo; a limitacao 2 deixa de existir.
#
# 2. Componente sistematico: a soma ponderada das caracteristicas do
#    talhao, chamada preditor linear e escrita como eta:
#       eta = beta0 + beta1 * x1 + beta2 * x2 + ...
#    onde x1, x2, ... sao as covariaveis (altitude, cultivar, ...) e os
#    betas sao os coeficientes a estimar. O eta pode ser qualquer numero.
#
# 3. Funcao de ligacao: a regra que conecta a probabilidade ao preditor
#    linear, g(p) = eta. Ela e escolhida de modo que, qualquer que seja o
#    valor de eta, a probabilidade fique entre 0 e 1; a limitacao 1 deixa
#    de existir.
#
# Para a familia binomial, as ligacoes disponiveis no glm() sao logit,
# probit, complemento log-log (cloglog) e cauchit. Todas levam a reta
# real para o intervalo (0, 1); diferem na forma da curva. A logit e a
# ligacao canonica da binomial (Step 0). Os slides (secao 9.1.2) mostram
# que so na logit o coeficiente exponenciado, exp(beta), tem um nome: a
# razao de chances. E essa a escala pedida pela Tarefa 4.
#
# OBJETIVO: comparar as quatro ligacoes nos dados e decidir qual adotar.
# O QUE O CODIGO FAZ: ajusta o mesmo modelo (todas as candidatas) com
#   cada ligacao, confere convergencia e possivel separacao, e compara a
#   log-verossimilhanca e o AIC.
# POR QUE ISSO E NECESSARIO: o requisito (a) do enunciado pede que a
#   familia e a ligacao sejam justificadas. A justificativa tem uma parte
#   estrutural (natureza da resposta) e uma parte empirica (os dados
#   preferem alguma ligacao?).
# COMO INTERPRETAR: diferencas de AIC menores que cerca de 2 sao
#   convencionalmente tratadas como equivalencia - os dados nao
#   distinguem as ligacoes. Nesse caso, a escolha se apoia na
#   interpretabilidade e na canonicidade, nao no ajuste.
#
# [APROFUNDAMENTO] Quando o AIC pode ser comparado.
# O AIC = -2 log-verossimilhanca + 2k (k = numero de parametros) so e
# comparavel entre modelos ajustados a mesma resposta, com a mesma
# verossimilhanca, nas mesmas observacoes. Aqui isso vale: a resposta e a
# mesma, a familia e a mesma, so a ligacao muda. No Desafio B (Tarefa 1.4)
# o AIC nao podia ser usado, porque as abordagens comparadas tinham
# respostas diferentes.
ligacoes <- c("logit", "probit", "cloglog", "cauchit")
modelos_por_ligacao <- lapply(ligacoes, function(L) {
    m <- glm(formula_completa, family = binomial(link = L), data = dados_cafe)
    verificar(m, paste0("ligacao ", L))
    m
})
names(modelos_por_ligacao) <- ligacoes

tab_ligacoes <- data.frame(
    ligacao = ligacoes,
    logLik  = sapply(modelos_por_ligacao, function(m) as.numeric(logLik(m))),
    AIC     = sapply(modelos_por_ligacao, AIC),
    fora    = sapply(modelos_por_ligacao, function(m) sum(fitted(m) <= 0 | fitted(m) >= 1)))
tab_ligacoes$dAIC <- tab_ligacoes$AIC - min(tab_ligacoes$AIC)
tab_ligacoes <- tab_ligacoes[order(tab_ligacoes$AIC), ]

cat("\ncomparacao das ligacoes (mesma resposta, mesmas candidatas):\n")
print(data.frame(ligacao = tab_ligacoes$ligacao,
                 logLik  = round(tab_ligacoes$logLik, 2),
                 AIC     = round(tab_ligacoes$AIC, 2),
                 dAIC    = round(tab_ligacoes$dAIC, 2),
                 fora_de_0_1 = tab_ligacoes$fora),
      row.names = FALSE)

n_equivalentes <- sum(tab_ligacoes$dAIC < 2)
dAIC_logit     <- tab_ligacoes$dAIC[tab_ligacoes$ligacao == "logit"]
cat(sprintf("\n   melhor por AIC: %s | ligacoes com dAIC < 2: %d de %d | dAIC da logit: %.2f\n",
            tab_ligacoes$ligacao[1], n_equivalentes, length(ligacoes), dAIC_logit))

# ---------------------------------------------------------------------
# As ligacoes dao respostas diferentes para o mesmo talhao?
# ---------------------------------------------------------------------
# OBJETIVO: ver se a escolha da ligacao muda a probabilidade predita de
#   um talhao concreto.
# O QUE O CODIGO FAZ: monta um talhao de referencia (medianas das
#   numericas e categoria de referencia dos fatores) e mais dois talhoes
#   reais - o de menor e o de maior probabilidade predita pela logit -
#   e calcula a probabilidade de cada um sob cada ligacao.
# POR QUE ISSO E NECESSARIO: as curvas das ligacoes sao quase iguais
#   perto de p = 0,5 e diferem nas pontas. Se os talhoes da base ficam
#   longe das pontas, a ligacao pouco importa para as predicoes.
# COMO INTERPRETAR: probabilidades parecidas entre ligacoes para o mesmo
#   talhao significam que a conclusao pratica nao depende da escolha.
talhao_referencia <- dados_cafe[1, c(fatores, numericas)]
for (v in numericas) talhao_referencia[[v]] <- median(dados_cafe[[v]])
for (f in fatores)   talhao_referencia[[f]] <- factor(levels(dados_cafe[[f]])[1], levels = levels(dados_cafe[[f]]))

id_min <- which.min(pred_logit)
id_max <- which.max(pred_logit)
talhoes_comparacao <- rbind(talhao_referencia,
                            dados_cafe[id_min, c(fatores, numericas)],
                            dados_cafe[id_max, c(fatores, numericas)])
rotulos_comparacao <- c("referencia (medianas + niveis de referencia)",
                        sprintf("%s (menor p pela logit)", dados_cafe$id_talhao[id_min]),
                        sprintf("%s (maior p pela logit)", dados_cafe$id_talhao[id_max]))

tab_pred_ligacoes <- sapply(modelos_por_ligacao, function(m)
    predict(m, newdata = talhoes_comparacao, type = "response"))
tab_pred_ligacoes <- data.frame(talhao = rotulos_comparacao, round(tab_pred_ligacoes, 3))
cat("\nprobabilidade predita para o mesmo talhao, por ligacao:\n")
print(tab_pred_ligacoes, row.names = FALSE)

maior_divergencia <- max(apply(tab_pred_ligacoes[, ligacoes], 1, function(r) max(r) - min(r)))
cat(sprintf("   maior diferenca entre ligacoes para um mesmo talhao: %.3f\n", maior_divergencia))

# ---------------------------------------------------------------------
# Figura C2: as curvas de ligacao
# ---------------------------------------------------------------------
# Pergunta que responde: o que cada funcao de ligacao faz com o preditor
# linear?
# Como ler: o eixo x e o preditor linear eta (qualquer numero); o eixo y
# e a probabilidade que cada ligacao devolve. Todas ficam entre 0 e 1.
# Logit e probit sao simetricas em torno de p = 0,5. A probit parece
# mais inclinada porque o eta dela esta em outra escala; ao ajustar o
# modelo, os coeficientes compensam essa diferenca, e as duas curvas
# ficam quase indistinguiveis (slides, figura 8). A cloglog e
# assimetrica (sobe mais rapido de um lado do que do outro); a cauchit
# tem caudas mais pesadas (demora mais para chegar a 0 e a 1).
# Limitacao: figura conceitual; nao usa os dados.
eta_grade <- seq(-5, 5, length.out = 400)
df_curvas <- rbind(
    data.frame(eta = eta_grade, p = plogis(eta_grade),           ligacao = "logit"),
    data.frame(eta = eta_grade, p = pnorm(eta_grade),            ligacao = "probit"),
    data.frame(eta = eta_grade, p = 1 - exp(-exp(eta_grade)),    ligacao = "cloglog"),
    data.frame(eta = eta_grade, p = pcauchy(eta_grade),          ligacao = "cauchit"))
df_curvas$ligacao <- factor(df_curvas$ligacao, levels = ligacoes)

g_curvas <- ggplot(df_curvas, aes(x = eta, y = p, colour = ligacao, linetype = ligacao)) +
    geom_hline(yintercept = c(0, 1), colour = CINZA, linewidth = .4) +
    geom_hline(yintercept = .5, colour = CINZA, linetype = "dotted") +
    geom_line(linewidth = 1) +
    scale_colour_manual(values = c(logit = VERDE, probit = "#4A7FA5",
                                   cloglog = "#8A6D3B", cauchit = VERMELHO)) +
    scale_linetype_manual(values = c(logit = "solid", probit = "solid",
                                     cloglog = "dashed", cauchit = "dotdash")) +
    labs(x = "preditor linear  eta", y = "probabilidade  p", colour = NULL, linetype = NULL,
         title = "O que cada funcao de ligacao faz com o preditor linear",
         subtitle = "todas mantem a probabilidade entre 0 e 1; diferem na forma da curva",
         caption = "logit e probit sao simetricas em torno de p = 0,5; cloglog e assimetrica; cauchit tem caudas pesadas") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
salvar(g_curvas, "C2_curvas_de_ligacao.png", 9, 6)

# ---------------------------------------------------------------------
# Parecer do Step 2
# ---------------------------------------------------------------------
# CONCLUSAO LITERAL:
# Resume o que a regressao linear fez de errado nestes dados.
#
# CONCLUSAO ESTATISTICA:
# Relaciona as duas limitacoes aos tres componentes do MLG e registra o
# resultado da comparacao empirica das ligacoes.
#
# DECISAO:
# Familia binomial com ligacao logit, com os motivos.
#
# O QUE FICA PARA OS PROXIMOS STEPS:
# A ligacao escolhida sera verificada com os dados no diagnostico.
parecer("CONCLUSAO LITERAL. A regressao linear ajustada com as ",
        length(numericas) + length(fatores), " candidatas produziu ", fora_dominio,
        " predicao(oes) fora do intervalo [0, 1] (a maior vale ",
        sprintf("%.3f", max(pred_linear)), "), e a variancia dos seus residuos ",
        "varia ", sprintf("%.2f", max(tab_variancia$var_residuos) / min(tab_variancia$var_residuos)),
        " vezes entre as faixas de predicao, acompanhando p(1-p) em vez de ser ",
        "constante. ",
        "CONCLUSAO ESTATISTICA. As duas limitacoes tem a mesma origem: a regressao ",
        "linear conecta a probabilidade diretamente a uma reta e supoe variancia ",
        "constante. O modelo linear generalizado separa o problema em tres ",
        "componentes: a familia binomial embute a variancia p(1-p); a funcao de ",
        "ligacao mantem a probabilidade entre 0 e 1 para qualquer valor do preditor ",
        "linear. Entre as ", length(ligacoes), " ligacoes disponiveis, ",
        n_equivalentes, " ficam dentro de dAIC < 2 (a melhor por AIC e ",
        tab_ligacoes$ligacao[1], "; a logit fica a ", sprintf("%.2f", dAIC_logit),
        " unidades), e a maior diferenca de probabilidade predita entre ligacoes ",
        "para um mesmo talhao e de ", sprintf("%.3f", maior_divergencia), ". ",
        if (n_equivalentes == length(ligacoes))
            "Os dados, portanto, nao distinguem as ligacoes. " else
            "Os dados distinguem parcialmente as ligacoes; a diferenca precisa ser discutida. ",
        "DECISAO. Adota-se a familia binomial com ligacao logit. A familia decorre ",
        "da natureza da resposta (um resultado 0/1 por talhao). A ligacao logit e ",
        "adotada por ser a ligacao canonica da binomial e por ser a unica em que ",
        "exp(beta) e uma razao de chances, a escala de interpretacao pedida pela ",
        "Tarefa 4; nao e adotada por vantagem de ajuste, que nao existe aqui. ",
        "O QUE FICA PARA OS PROXIMOS STEPS: a adequacao da ligacao logit sera ",
        "verificada com os proprios dados no diagnostico (Step 5), pelo metodo da ",
        "apostila (secao 3.5).")


# =====================================================================
# ---- STEP 3 ---- MODELO AMPLO: TODAS AS CANDIDATAS JUNTAS
# =====================================================================
secao("STEP 3 - MODELO AMPLO")

# Pergunta do step: colocando todas as caracteristicas candidatas de uma
# vez, o modelo explica algo alem da proporcao geral? Quais termos
# parecem importar? Ha caracteristicas que dizem a mesma coisa?
#
# Modelo amplo, em uma frase: e o modelo logistico com todas as
# candidatas ao mesmo tempo. E o ponto de partida da selecao de
# variaveis (requisito b do enunciado). Diferente do Step 1, aqui cada
# efeito e estimado levando em conta as demais caracteristicas: e o
# efeito condicional, nao o marginal.
#
# OBJETIVO: ajustar o modelo amplo e conferir que ele pode ser lido.
# O QUE O CODIGO FAZ: glm() binomial com ligacao logit e todas as
#   candidatas. O argumento na.action = na.fail faz o ajuste parar se
#   houver valor ausente, em vez de descartar linhas em silencio; isso
#   tambem e exigido pelo dredge() do Step 4.
# POR QUE ISSO E NECESSARIO: e o mesmo modelo logistico do Step 2, agora
#   ajustado como objeto de trabalho dos proximos steps.
# COMO INTERPRETAR: nada ainda; primeiro conferimos convergencia e
#   sinais de separacao (verificar()), depois lemos o summary().
modelo_amplo <- glm(formula_completa,
                    family    = binomial(link = "logit"),
                    data      = dados_cafe,
                    na.action = na.fail)
verificar(modelo_amplo, "modelo_amplo")

# ---------------------------------------------------------------------
# Como ler o summary()
# ---------------------------------------------------------------------
# Cada linha e um coeficiente. As colunas:
#   Estimate   : o efeito estimado na escala do logit (beta). Para um
#                fator, e a diferenca em relacao ao nivel de referencia;
#                para uma numerica, e a mudanca no logit por unidade.
#   Std. Error : a incerteza dessa estimativa (erro-padrao).
#   z value    : Estimate dividido por Std. Error.
#   Pr(>|z|)   : p-valor do teste "esse coeficiente e zero?".
# O teste da ultima coluna e chamado de teste de Wald (slides, secao
# 7.2). Ele e rapido, mas os slides avisam que nao e confiavel quando o
# coeficiente e grande em modulo (efeito Hauck-Donner); a comparacao com
# o teste da razao de verossimilhancas vem logo abaixo.
#
# Numero de coeficientes: 1 intercepto + 3 (cultivar tem 4 niveis) +
# 2 (manejo) + 1 (irrigacao) + 3 (regiao) + 9 numericas = 19. Cada fator
# com k niveis gasta k - 1 coeficientes, porque o nivel de referencia e
# a base de comparacao (ver mneis/03_por_que_19_parametros.R).
print(summary(modelo_amplo))
cat(sprintf("   coeficientes estimados: %d (1 + 3 + 2 + 1 + 3 + 9)\n", length(coef(modelo_amplo))))

# ---------------------------------------------------------------------
# O modelo amplo explica algo alem do modelo nulo?
# ---------------------------------------------------------------------
# Em linguagem simples: o modelo nulo da a mesma probabilidade a todos os
# talhoes; o amplo usa 13 caracteristicas. A pergunta e se o ganho de
# ajuste do amplo e grande demais para ser obra do acaso.
#
# Em termos estatisticos: analise de deviance (slides, secao 7.3;
# apostila, secao 2.9). A diferenca entre a deviance do modelo nulo e a
# do modelo amplo, D(nulo) - D(amplo), segue aproximadamente uma
# distribuicao qui-quadrado com q graus de liberdade, onde q e o numero
# de parametros a mais no modelo amplo. Esse e o teste da razao de
# verossimilhancas (LRT), que os slides chamam de padrao-ouro para
# comparar modelos aninhados (um dentro do outro). A aproximacao vale
# porque, na binomial, o parametro de dispersao e conhecido (phi = 1).
#
# OBJETIVO: testar se as 13 candidatas, juntas, melhoram o ajuste.
# O QUE O CODIGO FAZ: calcula a queda da deviance e os graus de liberdade
#   a mao, obtem o p-valor da qui-quadrado, e confere o mesmo teste com
#   anova(modelo_nulo, modelo_amplo, test = "Chisq").
# COMO INTERPRETAR: p pequeno indica que o conjunto de caracteristicas
#   explica parte de quem atinge o padrao. Nao diz quais delas importam.
deviance_amplo <- deviance(modelo_amplo)
gl_amplo       <- df.residual(modelo_amplo)
lrt_global     <- deviance_nula - deviance_amplo
gl_global      <- gl_nulo - gl_amplo
p_global       <- pchisq(lrt_global, gl_global, lower.tail = FALSE)

cat(sprintf("\n   deviance nula  : %.2f com %d gl\n", deviance_nula, gl_nulo))
cat(sprintf("   deviance amplo : %.2f com %d gl\n", deviance_amplo, gl_amplo))
cat(sprintf("   queda          : %.2f com %d gl  ->  p = %s\n",
            lrt_global, gl_global, format.pval(p_global, digits = 3)))
cat("\n   o mesmo teste pelo anova():\n")
print(anova(modelo_nulo, modelo_amplo, test = "Chisq"))

# [APROFUNDAMENTO] Deviance residual nao mede ajuste global em dados 0/1.
# E comum ver a regra "se a deviance residual for proxima dos graus de
# liberdade, o ajuste esta bom". Os slides (secao 7.4) advertem que ela
# NAO vale para dados binarios: com uma observacao por talhao, a
# deviance residual nao tem distribuicao qui-quadrado nem aproximada.
# Por isso a razao deviance/gl nao e usada aqui como medida de ajuste.
# O ajuste global sera avaliado no Step 5 por calibracao
# (Hosmer-Lemeshow) e por analise grafica de residuos, como os slides
# recomendam.

# ---------------------------------------------------------------------
# Qual termo faz falta? O teste de cada termo pelo LRT
# ---------------------------------------------------------------------
# Em linguagem simples: retirar uma caracteristica de cada vez e medir
# quanto o ajuste piora. Se piora mais do que o acaso explicaria, a
# caracteristica faz falta.
#
# OBJETIVO: obter, para cada termo, o teste da razao de verossimilhancas
#   condicional aos demais termos.
# O QUE O CODIGO FAZ: drop1(modelo_amplo, test = "Chisq") ajusta 13
#   modelos, cada um sem um dos termos, e compara cada um com o amplo
#   pela queda da deviance. Para um fator, o teste retira todos os
#   niveis de uma vez (um p por termo, com k - 1 graus de liberdade).
# POR QUE ISSO E NECESSARIO: o summary() da um p por coeficiente (Wald,
#   um por nivel); o drop1() da um p por termo (LRT). Sao perguntas
#   diferentes, e o LRT e o teste recomendado.
# COMO INTERPRETAR: p pequeno = o termo faz falta ao modelo amplo; p
#   grande = o termo pode sair sem perda perceptivel. E a base da
#   selecao do Step 4.
tab_drop1 <- drop1(modelo_amplo, test = "Chisq")
cat("\nteste da razao de verossimilhancas de cada termo (drop1):\n")
print(tab_drop1)

# ---------------------------------------------------------------------
# Wald e LRT: as duas perguntas
# ---------------------------------------------------------------------
# OBJETIVO: colocar lado a lado o p de Wald do summary() e o p do LRT do
#   drop1(), para ver onde concordam e onde nao.
# O QUE O CODIGO FAZ: para cada termo, recolhe o menor p de Wald entre os
#   seus coeficientes (um fator tem varios) e o p do LRT do termo.
# COMO INTERPRETAR: para numericas (1 coeficiente), os dois p sao quase
#   iguais. Para fatores, podem diferir: o Wald de um nivel pergunta "este
#   nivel difere da referencia?"; o LRT do termo pergunta "o fator inteiro
#   faz falta?". Um nivel pode diferir sem que o fator, como um todo,
#   melhore o ajuste de forma clara - e vice-versa.
coefs_amplo   <- summary(modelo_amplo)$coefficients
termos_amplo  <- attr(terms(modelo_amplo), "term.labels")
assign_amplo  <- attr(model.matrix(modelo_amplo), "assign")   # qual termo gera cada coluna

tab_wald_lrt <- data.frame(
    termo      = termos_amplo,
    gl         = tab_drop1$Df[-1],
    LRT        = tab_drop1$LRT[-1],
    p_LRT      = tab_drop1$`Pr(>Chi)`[-1],
    p_Wald_min = sapply(seq_along(termos_amplo), function(i)
        min(coefs_amplo[assign_amplo == i, "Pr(>|z|)"])))
tab_wald_lrt <- tab_wald_lrt[order(tab_wald_lrt$p_LRT), ]

cat("\np do LRT (por termo) e menor p de Wald (entre os coeficientes do termo):\n")
print(data.frame(termo      = tab_wald_lrt$termo,
                 gl         = tab_wald_lrt$gl,
                 LRT        = round(tab_wald_lrt$LRT, 3),
                 p_LRT      = round(tab_wald_lrt$p_LRT, 4),
                 p_Wald_min = round(tab_wald_lrt$p_Wald_min, 4)),
      row.names = FALSE)

# ---------------------------------------------------------------------
# Figura C3: quanto cada termo reduz a deviance
# ---------------------------------------------------------------------
# Pergunta que responde: quais termos mais fazem falta ao modelo amplo?
# Como ler: a barra e a estatistica do LRT de cada termo (queda da
# deviance ao retira-lo); quanto maior, mais o termo contribui. Barras
# de fatores tem mais graus de liberdade, por isso a cor marca o p-valor,
# que ja leva os graus de liberdade em conta.
# Limitacao: cada barra e condicional aos demais termos do amplo; ao
# retirar termos (Step 4), as barras dos que ficam podem mudar.
tab_wald_lrt$significativo <- ifelse(tab_wald_lrt$p_LRT < 0.05, "p < 0,05", "p >= 0,05")
g_lrt <- ggplot(tab_wald_lrt, aes(x = reorder(termo, LRT), y = LRT, fill = significativo)) +
    geom_col(colour = "black", width = .7) +
    geom_text(aes(label = sprintf("p = %.3f (gl %d)", p_LRT, gl)), hjust = -.1, size = 3.2) +
    coord_flip() +
    scale_fill_manual(values = c(`p < 0,05` = VERDE, `p >= 0,05` = "grey80"), name = "LRT") +
    scale_y_continuous(expand = expansion(mult = c(0, .35))) +
    labs(x = NULL, y = "estatistica do LRT (queda da deviance ao retirar o termo)",
         title = "Quanto cada termo faz falta ao modelo amplo",
         subtitle = "teste da razao de verossimilhancas de cada termo, condicional aos demais",
         caption = "fatores tem mais graus de liberdade; o p-valor ja considera isso") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
salvar(g_lrt, "C3_lrt_por_termo.png", 9, 6)

# ---------------------------------------------------------------------
# Multicolinearidade: caracteristicas que dizem a mesma coisa
# ---------------------------------------------------------------------
# Em linguagem simples: se duas caracteristicas andam juntas (por
# exemplo, umidade e precipitacao), o modelo tem dificuldade em separar
# o efeito de cada uma; as estimativas ficam instaveis e os erros-padrao
# crescem. O enunciado pede que isso seja verificado entre as variaveis
# de solo e clima.
#
# Em termos estatisticos: fator de inflacao de variancia (VIF) - quanto
# a variancia de um coeficiente esta inflada pela correlacao com as
# demais covariaveis. Para fatores com mais de um coeficiente, usa-se o
# VIF generalizado (GVIF; Fox e Monette, 1992, implementado em
# car::vif()). A coluna GVIF^(1/(2*Df)) coloca fatores e numericas na
# mesma escala: e comparavel a raiz quadrada de um VIF comum.
#
# OBJETIVO: verificar se alguma candidata e redundante com as demais.
# O QUE O CODIGO FAZ: car::vif() no modelo amplo; marca os termos cujo
#   GVIF ajustado passa de 2 (equivale a VIF > 4).
# COMO INTERPRETAR: o limite 2 e uma regra pratica, nao um criterio
#   formal. Valores bem abaixo dele indicam que os efeitos podem ser
#   estimados separadamente.
tab_vif <- as.data.frame(car::vif(modelo_amplo))
names(tab_vif) <- c("GVIF", "Df", "GVIF_ajustado")
tab_vif$termo <- rownames(tab_vif)
tab_vif$atencao <- ifelse(tab_vif$GVIF_ajustado > 2, "sim", "")
cat("\nfator de inflacao de variancia generalizado (GVIF):\n")
print(data.frame(termo         = tab_vif$termo,
                 GVIF          = round(tab_vif$GVIF, 3),
                 Df            = tab_vif$Df,
                 GVIF_ajustado = round(tab_vif$GVIF_ajustado, 3),
                 acima_de_2    = tab_vif$atencao),
      row.names = FALSE)
gvif_max <- tab_vif[which.max(tab_vif$GVIF_ajustado), ]
cat(sprintf("   maior GVIF ajustado: %.3f (%s); limite pratico: 2\n",
            gvif_max$GVIF_ajustado, gvif_max$termo))

# ---------------------------------------------------------------------
# Dispersao: o que se pode e o que nao se pode dizer com dados 0/1
# ---------------------------------------------------------------------
# No Desafio A, a razao entre a estatistica de Pearson e os graus de
# liberdade foi central: valores muito acima de 1 indicavam
# superdispersao. Aqui ela e calculada por exigencia do requisito (c),
# mas com uma ressalva dos slides (secao 9.1.4): com uma observacao 0/1
# por talhao, a variancia e sempre p(1-p), sem folga; a superdispersao
# nao e identificavel. A razao abaixo, portanto, nao pode acusar nem
# descartar superdispersao. A verificacao util da dispersao, para dados
# binarios, e a do envelope simulado (Step 5).
pearson_gl_amplo <- sum(residuals(modelo_amplo, type = "pearson")^2) / df.residual(modelo_amplo)
cat(sprintf("\n   estatistica de Pearson / gl (modelo amplo): %.3f\n", pearson_gl_amplo))
cat("   (nao interpretavel como dispersao em dados binarios nao agrupados)\n")

# ---------------------------------------------------------------------
# Parecer do Step 3
# ---------------------------------------------------------------------
# CONCLUSAO LITERAL:
# As caracteristicas, juntas, ajudam? Quais parecem fazer falta?
#
# CONCLUSAO ESTATISTICA:
# LRT global, LRT por termo, colinearidade e a ressalva da dispersao.
#
# RESPOSTA PARCIAL:
# Primeira lista de determinantes, agora ajustada mutuamente.
#
# O QUE FICA PARA OS PROXIMOS STEPS:
# Reduzir o modelo (Step 4) e diagnostica-lo (Step 5).
termos_fazem_falta  <- tab_wald_lrt$termo[tab_wald_lrt$p_LRT < 0.05]
termos_dispensaveis <- tab_wald_lrt$termo[tab_wald_lrt$p_LRT > 0.5]
# Numericas com tendencia no Step 1 (analise marginal) que nao fazem
# falta no modelo amplo (analise condicional), e vice-versa.
perdem_forca <- setdiff(as.character(num_tendencia$variavel), termos_fazem_falta)
ganham_forca <- setdiff(intersect(termos_fazem_falta, numericas), as.character(num_tendencia$variavel))

parecer("CONCLUSAO LITERAL. Com as ", length(termos_amplo), " caracteristicas juntas, o ",
        "modelo explica parte de quem atinge o padrao de exportacao: a queda da ",
        "deviance em relacao ao modelo nulo e de ", sprintf("%.1f", lrt_global),
        " com ", gl_global, " graus de liberdade (p = ", format.pval(p_global, digits = 2),
        "). Retirando um termo de cada vez, os que fazem falta ao ajuste (LRT com ",
        "p < 0,05) sao: ", paste(termos_fazem_falta, collapse = ", "),
        ". Os que podem sair sem perda perceptivel (p > 0,5) sao: ",
        paste(termos_dispensaveis, collapse = ", "), ". ",
        "CONCLUSAO ESTATISTICA. O teste global e o teste de cada termo sao testes ",
        "da razao de verossimilhancas (analise de deviance), validos porque o ",
        "parametro de dispersao da binomial e conhecido. Os p de Wald do summary() e ",
        "os p do LRT concordam para as numericas e divergem para alguns fatores, ",
        "porque respondem a perguntas diferentes (um nivel contra a referencia; o ",
        "fator inteiro contra o modelo sem ele). ",
        if (gvif_max$GVIF_ajustado > 2)
            paste0("Ha sinal de multicolinearidade: o maior GVIF ajustado e ",
                   sprintf("%.2f", gvif_max$GVIF_ajustado), " (", gvif_max$termo,
                   "), acima do limite pratico de 2; a selecao do Step 4 deve considerar isso. ")
        else
            paste0("Nao ha multicolinearidade relevante: o maior GVIF ajustado e ",
                   sprintf("%.2f", gvif_max$GVIF_ajustado), " (", gvif_max$termo,
                   "), abaixo do limite pratico de 2, de modo que os efeitos de solo e ",
                   "clima podem ser estimados separadamente. "),
        "A razao de Pearson sobre os graus de liberdade vale ",
        sprintf("%.2f", pearson_gl_amplo),
        ", mas nao e interpretavel como medida de dispersao em dados 0/1 nao agrupados. ",
        "RESPOSTA PARCIAL. Ajustadas umas pelas outras, ", length(termos_fazem_falta),
        " caracteristica(s) sustenta(m) evidencia clara. Das numericas com tendencia ",
        "marginal no Step 1, ",
        if (length(perdem_forca) > 0) paste0(paste(perdem_forca, collapse = ", "),
                                             " deixam de fazer falta no modelo conjunto")
        else "todas continuam a fazer falta no modelo conjunto",
        if (length(ganham_forca) > 0) paste0("; e ", paste(ganham_forca, collapse = ", "),
                                             " passam a fazer falta sem tendencia marginal clara")
        else "", ": a associacao marginal e o efeito condicional sao coisas diferentes. ",
        "O QUE FICA PARA OS PROXIMOS STEPS: reduzir o modelo amplo com criterio ",
        "explicito e caminho documentado (Step 4) e diagnosticar o modelo final antes ",
        "de interpreta-lo (Step 5).")


# =====================================================================
# ---- STEP 4 ---- SELECAO DE VARIAVEIS
# =====================================================================
secao("STEP 4 - SELECAO DE VARIAVEIS")

# Pergunta do step: quais caracteristicas podem sair do modelo sem perda
# de explicacao, e como registrar o caminho de modo que outra pessoa
# possa critica-lo?
#
# Parcimonia, em uma frase: um modelo com menos termos e mais facil de
# interpretar e menos sujeito a captar ruido; so vale manter o que paga
# o seu custo em ajuste.
#
# ---------------------------------------------------------------------
# Regra de decisao, escrita ANTES de olhar os resultados
# ---------------------------------------------------------------------
# 1. Criterio principal: eliminacao backward pelo teste da razao de
#    verossimilhancas (LRT) a 5%, como na Tarefa 6 do Desafio B. A cada
#    passo, retira-se o termo cuja retirada menos prejudica o modelo
#    (maior p no drop1); para-se quando toda retirada prejudicar (todos
#    os p < 0,05).
# 2. `cultivar` e mantida por desenho, sem competir na selecao: e a
#    caracteristica que o produtor escolhe e a que a Tarefa 2 compara.
#    E o mesmo tratamento dado a uma "exposicao de interesse" em
#    epidemiologia: o modelo precisa dela para responder a pergunta.
# 3. Criterio de conferencia: todos os subconjuntos por AICc e por BIC
#    (dredge). Um termo que o backward retirou so volta se (a) o modelo
#    com ele estiver dentro de dAICc < 2 do melhor E (b) o seu LRT for
#    significativo E (c) houver sentido agronomico. Sem os tres, e ruido.
# 4. O modelo final e testado contra o amplo por LRT: a simplificacao nao
#    pode ter custado ajuste.
# 5. "Irrelevante por construcao": termo com p > 0,5 no LRT do modelo
#    amplo E soma de pesos de Akaike abaixo de 0,30.
# Os limiares (5%, 2, 0,5, 0,30) sao convencoes praticas declaradas, nao
# teoremas.

# ---------------------------------------------------------------------
# Eliminacao backward por LRT
# ---------------------------------------------------------------------
# OBJETIVO: reduzir o modelo amplo um termo por vez, registrando cada
#   passo.
# O QUE O CODIGO FAZ: em cada volta do laco, roda drop1(test = "Chisq"),
#   ignora os termos fixos, procura o maior p; se ele for >= 0,05, retira
#   esse termo com update() e imprime o passo; senao, para.
# POR QUE ISSO E NECESSARIO: o enunciado pede o caminho documentado. O
#   laco imprime termo retirado, p e AIC antes e depois, para que a
#   decisao possa ser refeita e criticada.
# COMO INTERPRETAR: a ordem de saida vai do termo mais dispensavel ao
#   menos dispensavel; os p crescem ao longo do caminho porque, a cada
#   retirada, os termos restantes sao reavaliados.
ALFA_SELECAO  <- 0.05
TERMOS_FIXOS  <- "cultivar"

modelo_atual    <- modelo_amplo
caminho_backward <- data.frame()
passo <- 0
repeat {
    d <- drop1(modelo_atual, test = "Chisq")
    d <- d[-1, ]                                      # tira a linha <none>
    d <- d[!rownames(d) %in% TERMOS_FIXOS, ]
    if (nrow(d) == 0 || max(d[["Pr(>Chi)"]]) < ALFA_SELECAO) break
    remover <- rownames(d)[which.max(d[["Pr(>Chi)"]])]
    passo   <- passo + 1
    aic_antes <- AIC(modelo_atual)
    modelo_atual <- update(modelo_atual, as.formula(paste(". ~ . -", remover)))
    caminho_backward <- rbind(caminho_backward, data.frame(
        passo = passo, termo_retirado = remover,
        p_LRT = max(d[["Pr(>Chi)"]]), AIC_antes = aic_antes, AIC_depois = AIC(modelo_atual)))
    cat(sprintf("   passo %d: retira %-22s  p = %.4f   AIC %.2f -> %.2f\n",
                passo, remover, max(d[["Pr(>Chi)"]]), aic_antes, AIC(modelo_atual)))
}
modelo_final <- modelo_atual
verificar(modelo_final, "modelo_final")
termos_final <- attr(terms(modelo_final), "term.labels")

cat("\n   termos retirados:", nrow(caminho_backward), "| termos mantidos:", length(termos_final), "\n")
cat("   modelo final (backward):", deparse(formula(modelo_final)), "\n")
cat(sprintf("   AIC = %.2f | BIC = %.2f | deviance = %.2f com %d gl\n",
            AIC(modelo_final), BIC(modelo_final), deviance(modelo_final), df.residual(modelo_final)))
cat("\n   LRT de cada termo no modelo final:\n")
print(drop1(modelo_final, test = "Chisq"))

# ---------------------------------------------------------------------
# Segunda opiniao: todos os subconjuntos por AICc e por BIC
# ---------------------------------------------------------------------
# Em linguagem simples: em vez de retirar um termo de cada vez, ajustar
# todos os modelos possiveis com as candidatas e ordena-los por uma nota
# que premia o ajuste e penaliza o numero de parametros.
#
# Em termos estatisticos (slides, secao 7.5):
#   AIC = -2 log-verossimilhanca + 2k
#   BIC = -2 log-verossimilhanca + k log(n)
# onde k e o numero de parametros e n o numero de observacoes. O BIC
# penaliza mais (log(340) e maior que 2), entao tende a escolher modelos
# menores. O AICc e o AIC com uma correcao para amostras pequenas; com
# n = 340 a correcao e minima (Burnham e Anderson, 2002).
#
# Soma de pesos de Akaike (sw): cada modelo recebe um peso proporcional a
# exp(-dAICc/2); a soma dos pesos dos modelos que contem um termo mede
# com que frequencia esse termo aparece entre os bons modelos (0 a 1).
#
# OBJETIVO: conferir o backward com um criterio independente.
# O QUE O CODIGO FAZ: MuMIn::dredge() ajusta todos os subconjuntos das
#   candidatas com cultivar fixa (2^12 = 4096 modelos), ordena por AICc
#   e, separadamente, por BIC; sw() calcula a soma de pesos.
# COMO INTERPRETAR: se o modelo do backward esta no topo (dAICc < 2), os
#   dois criterios concordam. Se AICc e BIC discordam entre si, isso e
#   esperado (penalidades diferentes) e deve ser discutido, nao escondido.
cat("\najustando todos os subconjuntos (dredge) - leva alguns segundos...\n")
tempo_dredge <- system.time({
    dredge_aicc <- dredge(modelo_amplo, rank = "AICc", fixed = TERMOS_FIXOS, trace = FALSE)
    dredge_bic  <- dredge(modelo_amplo, rank = "BIC",  fixed = TERMOS_FIXOS, trace = FALSE)
})
cat(sprintf("   %d modelos ajustados em %.0f s\n", nrow(dredge_aicc), tempo_dredge[["elapsed"]]))

cat("\n   os 8 melhores por AICc:\n")
print(head(dredge_aicc, 8))
cat("\n   os 5 melhores por BIC:\n")
print(head(dredge_bic, 5))

tab_sw <- data.frame(termo = names(sw(dredge_aicc)), soma_pesos = as.vector(sw(dredge_aicc)))
tab_sw <- tab_sw[order(-tab_sw$soma_pesos), ]
cat("\n   soma de pesos de Akaike (AICc) por termo:\n")
print(data.frame(termo = tab_sw$termo, soma_pesos = round(tab_sw$soma_pesos, 3)), row.names = FALSE)

# ---------------------------------------------------------------------
# Os tres criterios concordam?
# ---------------------------------------------------------------------
# OBJETIVO: comparar o modelo do backward com o melhor por AICc e o melhor
#   por BIC.
# O QUE O CODIGO FAZ: recupera os dois melhores modelos do dredge, lista
#   os termos de cada um, e calcula o AICc do modelo do backward em
#   relacao ao melhor.
# COMO INTERPRETAR: termos presentes no melhor AICc mas ausentes no
#   backward sao candidatos a "voltar" - e passam pela regra 3.
melhor_aicc <- get.models(dredge_aicc, 1)[[1]]
melhor_bic  <- get.models(dredge_bic, 1)[[1]]
termos_melhor_aicc <- attr(terms(melhor_aicc), "term.labels")
termos_melhor_bic  <- attr(terms(melhor_bic),  "term.labels")
daicc_final <- AICc(modelo_final) - min(dredge_aicc$AICc)
dbic_final  <- BIC(modelo_final)  - min(dredge_bic$BIC)

cat("\n   backward    :", paste(termos_final, collapse = " + "), "\n")
cat("   melhor AICc :", paste(termos_melhor_aicc, collapse = " + "), "\n")
cat("   melhor BIC  :", paste(termos_melhor_bic,  collapse = " + "), "\n")
cat(sprintf("   modelo do backward: dAICc = %.2f em relacao ao melhor AICc; dBIC = %.2f em relacao ao melhor BIC\n",
            daicc_final, dbic_final))

so_no_aicc <- setdiff(termos_melhor_aicc, termos_final)
so_no_bw   <- setdiff(termos_final, termos_melhor_bic)
if (length(so_no_aicc) > 0) {
    cat("\n   termos no melhor AICc que o backward retirou:", paste(so_no_aicc, collapse = ", "), "\n")
    for (t in so_no_aicc) {
        m_t <- update(modelo_final, as.formula(paste(". ~ . +", t)))
        lrt_t <- anova(modelo_final, m_t, test = "Chisq")
        cat(sprintf("      %-22s LRT p = %.4f  (regra 3: so volta com p < %.2f e sentido agronomico)\n",
                    t, lrt_t[["Pr(>Chi)"]][2], ALFA_SELECAO))
    }
}
if (length(so_no_bw) > 0)
    cat("\n   termos no backward que o BIC dispensa:", paste(so_no_bw, collapse = ", "),
        "\n   (o BIC penaliza mais; mantidos pelo criterio principal, com LRT significativo)\n")

# ---------------------------------------------------------------------
# A simplificacao custou ajuste? LRT final x amplo
# ---------------------------------------------------------------------
# OBJETIVO: testar, de uma vez, se os termos retirados faziam falta em
#   conjunto.
# O QUE O CODIGO FAZ: anova(modelo_final, modelo_amplo, test = "Chisq") -
#   queda da deviance entre os dois, com gl = numero de parametros
#   retirados.
# COMO INTERPRETAR: p grande = os termos retirados nao explicavam nada
#   alem do ruido; p pequeno = simplificou demais.
lrt_final_amplo <- anova(modelo_final, modelo_amplo, test = "Chisq")
cat("\n   LRT: modelo final contra modelo amplo\n")
print(lrt_final_amplo)
p_final_amplo <- lrt_final_amplo[["Pr(>Chi)"]][2]

# ---------------------------------------------------------------------
# Irrelevantes por construcao
# ---------------------------------------------------------------------
# O enunciado avisa que ha variaveis geradas sem relacao com a resposta
# e pede que sejam identificadas. Regra 5: p > 0,5 no LRT do modelo amplo
# E soma de pesos < 0,30. Os dois criterios juntos evitam chamar de
# "irrelevante" um termo que so ficou fraco por acaso em um deles.
# Observacao: "irrelevante por construcao" e uma hipotese sobre como os
# dados foram gerados; o que se observa e apenas ausencia de evidencia
# de efeito nesta amostra.
tab_irrelevantes <- merge(tab_wald_lrt[, c("termo", "p_LRT")], tab_sw, by = "termo")
tab_irrelevantes$irrelevante <- tab_irrelevantes$p_LRT > 0.5 & tab_irrelevantes$soma_pesos < 0.30
tab_irrelevantes <- tab_irrelevantes[order(tab_irrelevantes$soma_pesos), ]
cat("\n   criterio 'irrelevante por construcao' (p_LRT > 0,5 e soma de pesos < 0,30):\n")
print(data.frame(termo       = tab_irrelevantes$termo,
                 p_LRT_amplo = round(tab_irrelevantes$p_LRT, 3),
                 soma_pesos  = round(tab_irrelevantes$soma_pesos, 3),
                 irrelevante = ifelse(tab_irrelevantes$irrelevante, "sim", "")),
      row.names = FALSE)
irrelevantes <- tab_irrelevantes$termo[tab_irrelevantes$irrelevante]

# ---------------------------------------------------------------------
# Verificacao de sensibilidade: a broca faria diferenca?
# ---------------------------------------------------------------------
# A contagem de brocas foi deixada fora das candidatas (cabecalho). Para
# mostrar que isso nao escondeu um efeito, ela e testada aqui como taxa
# de captura por unidade de esforco (armadilhas x dias), em escala
# logaritmica, com 0,5 somado a contagem para evitar log(0). Se o LRT
# nao rejeitar, a exclusao nao custou nada.
dados_cafe$log_taxa_broca <- log((dados_cafe$n_brocas_capturadas + 0.5) /
                                 (dados_cafe$n_armadilhas * dados_cafe$dias_exposicao))
modelo_com_broca <- update(modelo_final, . ~ . + log_taxa_broca)
lrt_broca <- anova(modelo_final, modelo_com_broca, test = "Chisq")
p_broca   <- lrt_broca[["Pr(>Chi)"]][2]
cat(sprintf("\n   sensibilidade: log da taxa de captura de broca no modelo final -> LRT p = %.3f\n", p_broca))

# ---------------------------------------------------------------------
# Figura C4: soma de pesos por termo
# ---------------------------------------------------------------------
# Pergunta que responde: com que frequencia cada termo aparece entre os
# melhores modelos?
# Como ler: barra = soma de pesos de Akaike (0 a 1); verde = termo no
# modelo final; linha tracejada = limite 0,30 da regra 5. Termos abaixo
# da linha e com p > 0,5 no amplo sao os "irrelevantes por construcao".
# Limitacao: a soma de pesos depende do conjunto de candidatas e do
# criterio (AICc); nao e um teste.
tab_sw$no_final <- ifelse(tab_sw$termo %in% termos_final, "no modelo final", "retirado")
g_sw <- ggplot(tab_sw, aes(x = reorder(termo, soma_pesos), y = soma_pesos, fill = no_final)) +
    geom_col(colour = "black", width = .7) +
    geom_hline(yintercept = 0.30, linetype = "dashed", colour = VERMELHO) +
    coord_flip() +
    scale_fill_manual(values = c(`no modelo final` = VERDE, retirado = "grey80"), name = NULL) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
    labs(x = NULL, y = "soma de pesos de Akaike (AICc)",
         title = "Com que frequencia cada termo aparece entre os melhores modelos",
         subtitle = sprintf("%d subconjuntos com cultivar fixa; linha tracejada: limite 0,30 do criterio de irrelevancia",
                            nrow(dredge_aicc)),
         caption = "cultivar tem soma 1 por ser fixa, nao por evidencia") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
salvar(g_sw, "C4_importancia_termos.png", 9, 6)

# ---------------------------------------------------------------------
# Parecer do Step 4
# ---------------------------------------------------------------------
# CONCLUSAO LITERAL:
# O que ficou, o que saiu, o que nunca importou.
#
# CONCLUSAO ESTATISTICA:
# Criterio, caminho, concordancia entre criterios, LRT final x amplo.
#
# CRITICA DO PROCESSO:
# O que o metodo nao garante.
#
# O QUE FICA PARA OS PROXIMOS STEPS:
# Diagnosticar o modelo final antes de interpreta-lo.
parecer("CONCLUSAO LITERAL. Das ", length(termos_amplo), " caracteristicas candidatas, ",
        length(termos_final), " ficam no modelo final: ",
        paste(termos_final, collapse = ", "), ". As demais sairam uma a uma, nesta ordem: ",
        paste(caminho_backward$termo_retirado, collapse = " > "), ". ",
        if (length(irrelevantes) > 0)
            paste0("Tres criterios apontam ", paste(irrelevantes, collapse = ", "),
                   " como candidatas a 'irrelevantes por construcao': sairam cedo, tem p > 0,5 ",
                   "no modelo amplo e soma de pesos abaixo de 0,30. ")
        else "Nenhum termo atende ao criterio de 'irrelevante por construcao'. ",
        "A contagem de brocas, testada como taxa de captura, nao acrescenta nada ao ",
        "modelo final (p = ", sprintf("%.2f", p_broca), "). ",
        "CONCLUSAO ESTATISTICA. Selecao por eliminacao backward com LRT a ",
        sprintf("%.0f%%", 100 * ALFA_SELECAO), ", mantendo cultivar por desenho, em ",
        nrow(caminho_backward), " passos. O modelo final fica a ",
        sprintf("%.2f", daicc_final), " unidades de AICc do melhor entre os ",
        nrow(dredge_aicc), " subconjuntos",
        if (length(so_no_aicc) > 0)
            paste0(" (o melhor acrescenta ", paste(so_no_aicc, collapse = ", "),
                   ", sem LRT significativo: nao volta pela regra 3)") else "",
        ". O BIC, que penaliza mais, prefere ", length(termos_melhor_bic), " termos (",
        paste(termos_melhor_bic, collapse = ", "), ")",
        if (length(so_no_bw) > 0)
            paste0("; ", paste(so_no_bw, collapse = ", "),
                   " sao mantidos pelo criterio principal, com LRT significativo e sentido agronomico")
        else "", ". O LRT do modelo final contra o amplo tem p = ",
        sprintf("%.2f", p_final_amplo),
        if (p_final_amplo >= 0.05) ": a simplificacao nao custou ajuste. "
        else ": a simplificacao custou ajuste e deve ser revista. ",
        "CRITICA DO PROCESSO. A selecao backward examina um caminho entre muitos; o ",
        "dredge confirma que esse caminho leva a um modelo equivalente ao melhor por ",
        "AICc, mas a discordancia com o BIC mostra que a fronteira entre 'fica' e ",
        "'sai' e uma escolha de criterio, nao um fato dos dados. Manter cultivar por ",
        "desenho e uma decisao substantiva, nao estatistica. Os p-valores dos termos ",
        "mantidos foram obtidos apos a selecao e sao, por isso, otimistas. E 'irrelevante ",
        "por construcao' e uma hipotese sobre a origem dos dados; o que se observa e ",
        "ausencia de evidencia nesta amostra. ",
        "O QUE FICA PARA OS PROXIMOS STEPS: o modelo final so sera interpretado depois ",
        "do diagnostico (Step 5): residuos, influencia, ligacao e escala das numericas.")


# =====================================================================
# ---- STEP 5 ---- DIAGNOSTICO DO MODELO FINAL
# =====================================================================
secao("STEP 5 - DIAGNOSTICO: RESIDUOS")

# Pergunta do step: o modelo escolhido e digno de confianca? Os erros
# que ele comete tem o padrao que a teoria preve, algum talhao sozinho
# comanda o resultado, a ligacao esta certa, e as caracteristicas
# numericas agem mesmo "em linha reta" no logit?
#
# Por que diagnosticar antes de interpretar: se o diagnostico mostrar
# que uma numerica precisa de outra escala, ou que um talhao muda os
# sinais, o modelo final muda - e tudo o que tivesse sido interpretado
# antes teria de ser refeito. E a ordem do ciclo de trabalho dos slides
# (secao 12.1): selecionar, diagnosticar, so entao interpretar.
#
# A semente e fixada de novo aqui para que este bloco seja reprodutivel
# mesmo se for executado isoladamente: residuos quantilicos e envelope
# simulado envolvem sorteio.
set.seed(20260920)

numericas_final <- intersect(termos_final, numericas)
fatores_final   <- intersect(termos_final, fatores)
prob_final      <- fitted(modelo_final)
eta_final       <- predict(modelo_final, type = "link")
n_obs           <- nrow(dados_cafe)
n_par_final     <- length(coef(modelo_final))

# ---------------------------------------------------------------------
# Residuos quantilicos aleatorizados
# ---------------------------------------------------------------------
# Por que os residuos comuns nao servem aqui (slides, secao 8.1): como a
# resposta so vale 0 ou 1, a diferenca "observado - previsto" forma duas
# faixas curvas, uma para os 0s e outra para os 1s, e nao uma nuvem. Um
# analista pode ler essas faixas como falta de ajuste; elas sao apenas a
# discretude da resposta. Isso vale para os residuos ordinarios, de
# Pearson e de deviance (apostila, secao 3.4.2).
#
# Residuo quantilico aleatorizado (Dunn e Smyth, 1996; slides, secao
# 8.1.1), em linguagem simples: trocar cada observacao pela "posicao"
# que ela ocupa na distribuicao que o modelo previu para ela, sorteando
# um ponto dentro do degrau que a resposta 0/1 cria, e converter essa
# posicao para a escala normal. Se o modelo estiver correto, o resultado
# e exatamente uma amostra da normal padrao N(0, 1) - sem faixas, sem
# degraus. Por depender de sorteio, muda se a semente mudar.
#
# OBJETIVO: obter residuos legiveis para uma resposta 0/1 e ver se eles
#   tem o padrao esperado.
# O QUE O CODIGO FAZ: statmod::qresiduals() calcula os residuos
#   quantilicos; o grafico os coloca contra o preditor linear com uma
#   curva suavizada; shapiro.test() testa a normalidade.
# COMO INTERPRETAR: nuvem sem tendencia em torno de zero, com a maior
#   parte dos pontos entre -2 e 2, e o padrao esperado. Uma curva
#   suavizada que se afasta de zero sugere problema de ligacao ou de
#   escala (verificados adiante). O teste de Shapiro-Wilk com n = 340 e
#   sensivel a desvios pequenos; o grafico pesa mais do que o p-valor.
res_quantilico <- qresiduals(modelo_final)
shapiro_rq     <- shapiro.test(res_quantilico)

cat(sprintf("residuos quantilicos aleatorizados: media = %.3f | desvio-padrao = %.3f\n",
            mean(res_quantilico), sd(res_quantilico)))
cat(sprintf("   fora de [-2, 2]: %d de %d (esperado sob N(0,1): cerca de %.0f)\n",
            sum(abs(res_quantilico) > 2), n_obs, 0.0455 * n_obs))
cat(sprintf("   Shapiro-Wilk: W = %.4f, p = %.3f (com n = %d, ler junto com o grafico)\n",
            shapiro_rq$statistic, shapiro_rq$p.value, n_obs))

# Sensibilidade ao sorteio. Como o residuo quantilico e aleatorizado, o
# desvio-padrao e a contagem fora de [-2, 2] mudam a cada semente. Para
# nao ler um sorteio extremo como sinal de problema (nem um sorteio
# benigno como prova de ajuste), os residuos sao sorteados de novo 50
# vezes e resume-se a faixa. A semente do script e restaurada em seguida.
N_SORTEIOS <- 50
sorteios_rq <- t(sapply(seq_len(N_SORTEIOS), function(s) {
    set.seed(s)
    r <- qresiduals(modelo_final)
    c(sd = sd(r), fora = sum(abs(r) > 2))
}))
set.seed(20260920)
cat(sprintf("   em %d sorteios: desvio-padrao de %.2f a %.2f (media %.2f); fora de [-2, 2] de %d a %d (media %.1f)\n",
            N_SORTEIOS, min(sorteios_rq[, "sd"]), max(sorteios_rq[, "sd"]), mean(sorteios_rq[, "sd"]),
            min(sorteios_rq[, "fora"]), max(sorteios_rq[, "fora"]), mean(sorteios_rq[, "fora"])))
cat("   (o sorteio da semente do script e um entre esses; o envelope abaixo e a verificacao formal)\n")

# Figura C5: residuos quantilicos contra o preditor linear, e Q-Q
# Pergunta que responde: os erros do modelo tem padrao?
# Como ler (esquerda): eixo x = preditor linear ajustado (logit); eixo y
# = residuo quantilico; linhas em -2 e 2; a curva verde e um suavizador
# (loess) com faixa de incerteza. Curva plana perto de zero = bom.
# Como ler (direita): quantis dos residuos contra os de uma normal
# padrao; pontos sobre a reta = compativel com N(0, 1).
# Limitacao: depende da semente; o Q-Q sem envelope nao separa desvio
# real de variacao amostral - por isso o envelope vem a seguir.
df_rq <- data.frame(eta = eta_final, rq = res_quantilico)
g_rq_eta <- ggplot(df_rq, aes(x = eta, y = rq)) +
    geom_hline(yintercept = 0, colour = CINZA) +
    geom_hline(yintercept = c(-2, 2), colour = VERMELHO, linetype = "dashed") +
    geom_point(shape = 21, fill = "white", colour = "black", size = 2, alpha = .6) +
    geom_smooth(method = "loess", formula = y ~ x, colour = VERDE, fill = VERDE, alpha = .15, linewidth = .9) +
    labs(x = "preditor linear ajustado (logit)", y = "residuo quantilico aleatorizado",
         title = "Residuos quantilicos x preditor linear",
         subtitle = "esperado: nuvem sem tendencia, curva suavizada perto de zero") +
    theme_bw(base_size = 11)
g_rq_qq <- ggplot(df_rq, aes(sample = rq)) +
    stat_qq_line(colour = VERMELHO, linetype = "dashed") +
    stat_qq(shape = 21, fill = "white", colour = "black", size = 2, alpha = .6) +
    labs(x = "quantis teoricos N(0, 1)", y = "quantis dos residuos",
         title = "Q-Q dos residuos quantilicos",
         subtitle = sprintf("Shapiro-Wilk p = %.3f", shapiro_rq$p.value)) +
    theme_bw(base_size = 11)
salvar(g_rq_eta | g_rq_qq, "C5_residuos_quantilicos.png", 12, 5.5)

# ---------------------------------------------------------------------
# Envelope simulado (meio-normal)
# ---------------------------------------------------------------------
# Em linguagem simples (slides, secao 8.2): um grafico quantil-quantil
# sozinho nao diz o que e desvio real e o que e variacao amostral. O
# envelope resolve isso: simulam-se muitas bases de dados a partir do
# proprio modelo ajustado, calculam-se os residuos de cada uma, e
# desenha-se a faixa que contem 95% deles. Pontos dentro da faixa sao
# compativeis com o modelo; pontos fora, em grupo e de forma
# sistematica, indicam familia ou ligacao erradas. Os slides o chamam de
# "o diagnostico mais informativo disponivel para MLG, e o mais
# subutilizado".
#
# OBJETIVO: verificar se a familia binomial com ligacao logit e
#   compativel com os residuos observados.
# O QUE O CODIGO FAZ: hnp() simula 99 bases sob o modelo final, calcula
#   os residuos de deviance de cada uma, ordena os valores absolutos e
#   monta a faixa de 95%. Os residuos observados sao desenhados por cima.
# COMO INTERPRETAR: conta-se quantos pontos ficam fora da faixa; alguns
#   pontos isolados sao esperados (cerca de 5%); um grupo nas caudas nao.
#
# [APROFUNDAMENTO] O envelope simula do modelo ajustado. Isso pressupoe
# que o modelo esta correto para gerar a faixa; ele detecta
# incompatibilidade entre residuos e familia/ligacao, mas nao detecta
# uma covariavel omitida que nao altere a forma dos residuos.
envelope <- hnp(modelo_final, resid.type = "deviance", sim = 99, conf = 0.95,
                how.many.out = TRUE, print.on = FALSE, plot.sim = FALSE)
cat(sprintf("\nenvelope simulado (99 simulacoes, faixa de 95%%): %d de %d pontos fora (%.1f%%)\n",
            envelope$out, envelope$total, 100 * envelope$out / envelope$total))

# Figura C5: envelope
# Pergunta que responde: familia e ligacao estao certas?
# Como ler: eixo x = quantis meio-normais; eixo y = residuos absolutos
# ordenados; faixa cinza = envelope de 95%; linha tracejada = mediana
# simulada; pontos = residuos observados.
df_envelope <- data.frame(x = envelope$x, mediana = envelope$median,
                          inferior = envelope$lower, superior = envelope$upper,
                          residuo = envelope$residuals)
g_envelope <- ggplot(df_envelope, aes(x = x)) +
    geom_ribbon(aes(ymin = inferior, ymax = superior), fill = "grey80", alpha = .6) +
    geom_line(aes(y = mediana), colour = VERDE, linetype = "dashed") +
    geom_point(aes(y = residuo), shape = 21, fill = "white", colour = "black", size = 2, alpha = .6) +
    labs(x = "quantis meio-normais", y = "residuos de deviance absolutos (ordenados)",
         title = "Envelope simulado meio-normal do modelo final",
         subtitle = sprintf("%d de %d pontos fora da faixa de 95%% (99 simulacoes sob o modelo ajustado)",
                            envelope$out, envelope$total),
         caption = "pontos fora em grupo, sobretudo nas caudas, indicariam familia ou ligacao inadequadas") +
    theme_bw(base_size = 12)
salvar(g_envelope, "C5_envelope_simulado.png", 8, 6)

# ---------------------------------------------------------------------
# Dispersao
# ---------------------------------------------------------------------
# Como no Step 3: a razao de Pearson sobre os graus de liberdade e
# reportada por exigencia do requisito (c), com a ressalva de que, em
# dados 0/1 nao agrupados, a superdispersao nao e identificavel (slides,
# secao 9.1.4). A verificacao util e o envelope acima.
pearson_gl_final <- sum(residuals(modelo_final, type = "pearson")^2) / df.residual(modelo_final)
cat(sprintf("\ndispersao: estatistica de Pearson / gl (modelo final) = %.3f\n", pearson_gl_final))
cat("   (nao interpretavel como dispersao em dados binarios nao agrupados; ver envelope)\n")


# =====================================================================
# ---- STEP 5 ---- DIAGNOSTICO: AJUSTE GLOBAL (CALIBRACAO)
# =====================================================================
secao("STEP 5 - DIAGNOSTICO: CALIBRACAO")

# Em linguagem simples: dividir os talhoes em dez grupos pela
# probabilidade prevista e comparar, em cada grupo, quantos atingiram o
# padrao de fato com quantos o modelo esperava. Se o modelo diz "30%" a
# um grupo e cerca de 30% dele atingiu o padrao, o modelo esta
# calibrado.
#
# Em termos estatisticos: e o teste de Hosmer-Lemeshow, que os slides
# (secao 7.4) recomendam para dados binarios no lugar da deviance
# residual. A estatistica soma, sobre os g grupos, (O - E)^2 / (E(1 - p))
# com O = observados, E = esperados e p = probabilidade media do grupo;
# sob calibracao adequada ela segue aproximadamente qui-quadrado com
# g - 2 graus de liberdade (Hosmer, Lemeshow e Sturdivant, 2013).
#
# OBJETIVO: avaliar o ajuste global do modelo final.
# O QUE O CODIGO FAZ: corta as probabilidades ajustadas em decis, calcula
#   observados, esperados e a estatistica; desenha observados contra
#   esperados por decil.
# COMO INTERPRETAR: p grande = sem evidencia de descalibracao; pontos
#   perto da diagonal no grafico = o modelo "vale o que diz".
# Limitacao: o teste depende de como os grupos sao formados; e um
# indicador, nao uma prova de ajuste.
G_CALIBRACAO <- 10
decil <- cut(prob_final, quantile(prob_final, seq(0, 1, length.out = G_CALIBRACAO + 1)),
             include.lowest = TRUE)
tab_calibracao <- data.frame(
    grupo      = seq_len(nlevels(decil)),
    faixa_p    = levels(decil),
    n          = as.vector(table(decil)),
    p_media    = as.vector(tapply(prob_final, decil, mean)),
    observados = as.vector(tapply(resposta_exportacao, decil, sum)))
tab_calibracao$esperados <- tab_calibracao$n * tab_calibracao$p_media
tab_calibracao$prop_obs  <- tab_calibracao$observados / tab_calibracao$n
ic_cal <- mapply(function(x, n) prop.test(x, n, correct = FALSE)$conf.int,
                 tab_calibracao$observados, tab_calibracao$n)
tab_calibracao$ic_inf <- ic_cal[1, ]
tab_calibracao$ic_sup <- ic_cal[2, ]

estat_HL <- with(tab_calibracao, sum((observados - esperados)^2 / (esperados * (1 - p_media))))
gl_HL    <- nrow(tab_calibracao) - 2
p_HL     <- pchisq(estat_HL, gl_HL, lower.tail = FALSE)

print(data.frame(grupo      = tab_calibracao$grupo,
                 n          = tab_calibracao$n,
                 p_media    = round(tab_calibracao$p_media, 3),
                 esperados  = round(tab_calibracao$esperados, 1),
                 observados = tab_calibracao$observados,
                 prop_obs   = round(tab_calibracao$prop_obs, 3)),
      row.names = FALSE)
cat(sprintf("\n   Hosmer-Lemeshow: X2 = %.2f com %d gl, p = %.3f\n", estat_HL, gl_HL, p_HL))

# Figura C5: calibracao
# Pergunta que responde: a probabilidade prevista vale o que diz?
# Como ler: eixo x = probabilidade media prevista no decil; eixo y =
# proporcao observada (com IC de Wilson); a diagonal e a calibracao
# perfeita. Pontos cujo IC cruza a diagonal sao compativeis com ela.
g_calibracao <- ggplot(tab_calibracao, aes(x = p_media, y = prop_obs)) +
    geom_abline(slope = 1, intercept = 0, colour = VERMELHO, linetype = "dashed") +
    geom_errorbar(aes(ymin = ic_inf, ymax = ic_sup), width = 0, colour = CINZA) +
    geom_line(colour = VERDE, linewidth = .5) +
    geom_point(aes(size = n), colour = VERDE) +
    scale_size_continuous(range = c(2, 4), guide = "none") +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "probabilidade media prevista no decil", y = "proporcao observada no decil",
         title = "Calibracao do modelo final por decil de probabilidade prevista",
         subtitle = sprintf("Hosmer-Lemeshow: X2 = %.2f, gl = %d, p = %.3f", estat_HL, gl_HL, p_HL),
         caption = "diagonal = calibracao perfeita; barras = IC de Wilson 95%") +
    theme_bw(base_size = 12)
salvar(g_calibracao, "C5_calibracao.png", 7, 7)


# =====================================================================
# ---- STEP 5 ---- DIAGNOSTICO: TALHOES INFLUENTES
# =====================================================================
secao("STEP 5 - DIAGNOSTICO: INFLUENCIA")

# Duas ideias, uma de cada vez (slides, secao 8.3; apostila, 3.3.3):
#
# Alavancagem: quanto um talhao e "incomum" nas suas caracteristicas
# (por exemplo, o unico Bourbon a 1300 m). Talhoes incomuns tem mais
# poder de puxar o ajuste. Em um MLG, a alavancagem depende tambem da
# probabilidade ajustada, nao so do delineamento.
#
# Distancia de Cook: quanto as estimativas do modelo mudariam se aquele
# talhao fosse retirado. E a medida de influencia propriamente dita.
#
# Conduta correta (slides): identificar o ponto influente pelo seu
# identificador, reajustar o modelo sem ele, comparar as estimativas e
# relatar as duas versoes. Identificar nao autoriza remover.
#
# OBJETIVO: saber se algum talhao, sozinho, comanda o resultado.
# O QUE O CODIGO FAZ: calcula alavancagem e distancia de Cook; lista os
#   cinco talhoes mais influentes com suas caracteristicas; reajusta o
#   modelo sem os talhoes acima do limite pratico 4/n e compara os
#   coeficientes.
# COMO INTERPRETAR: a regra 4/n e apenas uma regua para olhar os
#   maiores; o que importa e se o valor absoluto de Cook e grande e se
#   os coeficientes mudam de sinal ou de conclusao sem esses talhoes.
alavancagem <- hatvalues(modelo_final)
dist_cook   <- cooks.distance(modelo_final)
limite_cook <- 4 / n_obs
limite_hat  <- 2 * n_par_final / n_obs

tab_influencia <- data.frame(id_talhao = dados_cafe$id_talhao,
                             dados_cafe[, termos_final],
                             prob     = round(prob_final, 3),
                             resposta = resposta_exportacao,
                             cook     = dist_cook,
                             alavanc  = alavancagem)
tab_influencia <- tab_influencia[order(-tab_influencia$cook), ]
acima_cook <- which(dist_cook > limite_cook)
acima_hat  <- which(alavancagem > limite_hat)

cat(sprintf("distancia de Cook: maxima = %.4f | acima de 4/n = %.4f: %d talhoes\n",
            max(dist_cook), limite_cook, length(acima_cook)))
cat(sprintf("alavancagem: maxima = %.3f | acima de 2p/n = %.3f: %d talhoes\n",
            max(alavancagem), limite_hat, length(acima_hat)))
cat("\n   os 5 talhoes mais influentes (distancia de Cook):\n")
tab_top5 <- head(tab_influencia, 5)
tab_top5$cook    <- round(tab_top5$cook, 4)
tab_top5$alavanc <- round(tab_top5$alavanc, 3)
print(tab_top5, row.names = FALSE)

# Reajuste sem os talhoes acima de 4/n e comparacao dos coeficientes
modelo_sem_influentes <- update(modelo_final, data = dados_cafe[-acima_cook, ])
tab_comparacao <- data.frame(
    coeficiente = names(coef(modelo_final)),
    com_todos   = coef(modelo_final),
    sem_influentes = coef(modelo_sem_influentes))
tab_comparacao$mudanca_pct <- 100 * (tab_comparacao$sem_influentes - tab_comparacao$com_todos) /
                              abs(tab_comparacao$com_todos)
tab_comparacao$mesmo_sinal <- sign(tab_comparacao$com_todos) == sign(tab_comparacao$sem_influentes)
cat(sprintf("\n   reajuste sem os %d talhoes acima de 4/n (%s):\n",
            length(acima_cook), paste(dados_cafe$id_talhao[acima_cook], collapse = ", ")))
print(data.frame(coeficiente    = tab_comparacao$coeficiente,
                 com_todos      = round(tab_comparacao$com_todos, 4),
                 sem_influentes = round(tab_comparacao$sem_influentes, 4),
                 mudanca_pct    = round(tab_comparacao$mudanca_pct, 1),
                 mesmo_sinal    = ifelse(tab_comparacao$mesmo_sinal, "sim", "NAO")),
      row.names = FALSE)
trocas_sinal <- sum(!tab_comparacao$mesmo_sinal[-1])   # ignora o intercepto
maior_mudanca <- max(abs(tab_comparacao$mudanca_pct[-1]))
cat(sprintf("   coeficientes que trocam de sinal: %d | maior mudanca relativa: %.0f%%\n",
            trocas_sinal, maior_mudanca))

# Figura C5: influencia
# Pergunta que responde: algum talhao comanda o resultado?
# Como ler: eixo x = posicao do talhao na base; eixo y = distancia de
# Cook; linha tracejada = regua 4/n; os cinco maiores estao rotulados
# pelo id_talhao. Olhar o valor absoluto do eixo y, nao so quem passa
# da linha.
df_cook <- data.frame(indice = seq_len(n_obs), cook = dist_cook, id = dados_cafe$id_talhao)
df_cook$rotulo <- ifelse(df_cook$id %in% tab_top5$id_talhao, as.character(df_cook$id), "")
# rotulos alternam o lado (esquerda/direita) pela ordem de influencia, para nao se sobrepor
df_cook$lado <- ifelse(match(df_cook$id, tab_top5$id_talhao) %% 2 == 1, 1.15, -0.15)
g_influencia <- ggplot(df_cook, aes(x = indice, y = cook)) +
    geom_hline(yintercept = limite_cook, colour = VERMELHO, linetype = "dashed") +
    geom_segment(aes(xend = indice, yend = 0), colour = CINZA, linewidth = .3) +
    geom_point(colour = VERDE, size = 1.8) +
    geom_text(aes(label = rotulo, hjust = lado), vjust = -.4, size = 3, colour = VERMELHO) +
    scale_y_continuous(expand = expansion(mult = c(0, .15))) +
    labs(x = "talhao (ordem na base)", y = "distancia de Cook",
         title = "Distancia de Cook por talhao",
         subtitle = sprintf("regua 4/n = %.4f (tracejada); maxima = %.4f; %d talhoes acima da regua",
                            limite_cook, max(dist_cook), length(acima_cook)),
         caption = "identificar nao autoriza remover: o reajuste sem esses talhoes e reportado no console") +
    theme_bw(base_size = 12)
salvar(g_influencia, "C5_influencia.png", 10, 5.5)


# =====================================================================
# ---- STEP 5 ---- DIAGNOSTICO: LIGACAO E ESCALA DAS COVARIAVEIS
# =====================================================================
secao("STEP 5 - DIAGNOSTICO: LIGACAO E ESCALA")

# ---------------------------------------------------------------------
# A funcao de ligacao esta adequada?
# ---------------------------------------------------------------------
# Em linguagem simples: se a curva que liga o preditor linear a
# probabilidade estivesse errada, o proprio preditor linear elevado ao
# quadrado "explicaria" parte do que sobrou.
#
# Em termos estatisticos (apostila, secao 3.5): adiciona-se eta^2 (o
# preditor linear ajustado ao quadrado) como covariavel extra e
# examina-se a mudanca na deviance pelo LRT. Se a queda for grande, ha
# evidencia de que a funcao de ligacao e insatisfatoria.
#
# OBJETIVO: testar a ligacao logit com os dados.
# O QUE O CODIGO FAZ: guarda o preditor linear ajustado na base, ajusta
#   o modelo final mais o termo eta^2 e compara por anova(test = "Chisq").
# COMO INTERPRETAR: p grande = sem evidencia contra a ligacao logit.
dados_cafe$eta_final <- eta_final
modelo_eta2 <- update(modelo_final, . ~ . + I(eta_final^2))
lrt_eta2    <- anova(modelo_final, modelo_eta2, test = "Chisq")
p_eta2      <- lrt_eta2[["Pr(>Chi)"]][2]
cat(sprintf("verificacao da ligacao (eta^2 como covariavel extra): LRT = %.3f com %d gl, p = %.3f\n",
            lrt_eta2$Deviance[2], lrt_eta2$Df[2], p_eta2))

# ---------------------------------------------------------------------
# A escala de cada numerica esta adequada?
# ---------------------------------------------------------------------
# Em linguagem simples: a altitude entra no modelo "como esta", somando
# a mesma quantidade ao logit por metro. Mas poderia ser que o efeito
# fosse maior em altitudes baixas e menor nas altas (ou o contrario),
# pedindo log, quadrado ou outra forma. Tres ferramentas, da informal a
# mais formal:
#
# (i) Residuos parciais (apostila, secao 3.7): para cada covariavel, um
#     residuo que "devolve" o efeito dela, u = z - eta + beta * x, onde
#     z e a resposta de trabalho do ajuste. O grafico de u contra x deve
#     ser aproximadamente reto se a escala esta boa; a forma da curva
#     sugere a transformacao, se houver.
#
# (ii) Variavel construida (apostila, secao 3.7): adiciona-se ao modelo o
#     termo x * log(x) e testa-se o seu coeficiente por LRT. Se nao for
#     significativo, a escala linear e adequada. Esse e o teste de
#     Box-Tidwell (1962), que corresponde a variavel construida da
#     familia de potencias.
#
# (iii) Suavizador (modelo aditivo generalizado, GAM; Wood, 2017): troca-se
#     o termo linear beta * x por uma curva flexivel s(x) e le-se o
#     numero efetivo de graus de liberdade (edf). Um edf proximo de 1
#     significa que a curva e, na pratica, uma reta.
#
# OBJETIVO: verificar a escala das numericas do modelo final.
# O QUE O CODIGO FAZ: residuals(type = "partial") para (i); update() com
#   I(x * log(x)) e anova() para (ii); mgcv::gam() com s(x) para (iii).
# COMO INTERPRETAR: (i) reto, (ii) p grande e (iii) edf perto de 1
#   apontam juntos para a escala linear. Se discordarem, o grafico (i)
#   e a curva do GAM dizem qual transformacao considerar.
res_parciais <- residuals(modelo_final, type = "partial")
df_parciais <- bind_rows(lapply(numericas_final, function(v)
    data.frame(variavel = v, x = dados_cafe[[v]], u = res_parciais[, v])))
df_parciais$variavel <- factor(df_parciais$variavel, levels = numericas_final)

tab_escala <- bind_rows(lapply(numericas_final, function(v) {
    m_bt <- update(modelo_final, as.formula(sprintf(". ~ . + I(%s * log(%s))", v, v)))
    lrt  <- anova(modelo_final, m_bt, test = "Chisq")
    data.frame(variavel = v, LRT_x_logx = lrt$Deviance[2], p_x_logx = lrt[["Pr(>Chi)"]][2])
}))

formula_gam <- as.formula(paste("padrao_exportacao ~",
                                paste(c(fatores_final, sprintf("s(%s)", numericas_final)), collapse = " + ")))
modelo_gam  <- mgcv::gam(formula_gam, family = binomial, data = dados_cafe)
edf_gam     <- summary(modelo_gam)$s.table[, "edf"]
tab_escala$edf_gam <- as.vector(edf_gam)

cat("\nescala das numericas do modelo final:\n")
print(data.frame(variavel   = tab_escala$variavel,
                 LRT_x_logx = round(tab_escala$LRT_x_logx, 3),
                 p_x_logx   = round(tab_escala$p_x_logx, 3),
                 edf_gam    = round(tab_escala$edf_gam, 2)),
      row.names = FALSE)
cat("   p_x_logx: teste da variavel construida x*log(x) (Box-Tidwell); edf_gam: graus de liberdade do suavizador\n")

# Figura C5: residuos parciais
# Pergunta que responde: a escala de cada numerica esta certa?
# Como ler: eixo x = covariavel; eixo y = residuo parcial; a reta
# vermelha e a relacao linear que o modelo assume; a curva verde e um
# suavizador. Curva colada a reta = escala adequada; curva que se afasta
# de forma sistematica = considerar transformacao.
# Os pontos formam duas faixas (uma dos talhoes com 0, outra dos com 1):
# e a discretude da resposta, a mesma dos residuos comuns, e nao falta
# de ajuste. Nao se leem os pontos; le-se a curva suavizada contra a reta.
# Limitacao: distorcoes podem aparecer se a escala de outra covariavel
# estiver errada; por isso se olham todas.
g_parciais <- ggplot(df_parciais, aes(x = x, y = u)) +
    geom_point(shape = 21, fill = "white", colour = "black", size = 1.8, alpha = .5) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = VERMELHO, linewidth = .8) +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE, colour = VERDE, fill = VERDE,
                alpha = .15, linewidth = .9) +
    facet_wrap(~ variavel, scales = "free", ncol = length(numericas_final)) +
    labs(x = NULL, y = "residuo parcial",
         title = "Residuos parciais das numericas do modelo final",
         subtitle = "reta vermelha: escala linear assumida; curva verde: suavizador com faixa de incerteza",
         caption = "metodo da apostila (secao 3.7); complementado pelo teste de Box-Tidwell e pelo edf do GAM") +
    theme_bw(base_size = 11)
salvar(g_parciais, "C5_residuos_parciais.png", 11, 5)

# ---------------------------------------------------------------------
# Decisao: o modelo final se mantem?
# ---------------------------------------------------------------------
# Cada verificacao vira um sinal. Os limiares sao praticos: envelope com
# mais de 5% dos pontos fora; Hosmer-Lemeshow, eta^2 ou x*log(x) com
# p < 0,05; edf do GAM acima de 2; troca de sinal de algum coeficiente
# sem os talhoes influentes. Um sinal nao muda o modelo automaticamente:
# ele exige que a revisao seja discutida e registrada.
sinais <- c(
    envelope    = envelope$out / envelope$total > 0.05,
    calibracao  = p_HL < 0.05,
    influencia  = trocas_sinal > 0,
    ligacao     = p_eta2 < 0.05,
    escala_BT   = any(tab_escala$p_x_logx < 0.05),
    escala_GAM  = any(tab_escala$edf_gam > 2))
cat("\nsinais de alerta do diagnostico:\n")
for (s in names(sinais))
    cat(sprintf("   %-12s %s\n", s, if (sinais[[s]]) "[ATENCAO]" else "[ok]"))
modelo_confirmado <- !any(sinais)
cat(if (modelo_confirmado) "\n   modelo final CONFIRMADO: segue para a interpretacao.\n"
    else "\n   modelo final com sinais de alerta: revisar antes de interpretar.\n")

# ---------------------------------------------------------------------
# Parecer do Step 5
# ---------------------------------------------------------------------
# RESIDUOS / DISPERSAO / CALIBRACAO / INFLUENCIA / LIGACAO / ESCALA:
# um resultado por bloco, com o numero que o sustenta.
#
# CONCLUSAO:
# Modelo confirmado ou revisado.
#
# O QUE FICA PARA OS PROXIMOS STEPS:
# Interpretar (Steps 6 a 8) e decidir (Step 9).
parecer("RESIDUOS. Os residuos quantilicos aleatorizados do modelo final tem media ",
        sprintf("%.2f", mean(res_quantilico)), " e desvio-padrao ",
        sprintf("%.2f", sd(res_quantilico)), "; ",
        sum(abs(res_quantilico) > 2), " de ", n_obs, " ficam fora de [-2, 2] e o ",
        "Shapiro-Wilk da p = ", sprintf("%.3f", shapiro_rq$p.value),
        "; como o residuo e aleatorizado, esses resumos variam com o sorteio (em ",
        N_SORTEIOS, " sorteios, desvio-padrao medio de ",
        sprintf("%.2f", mean(sorteios_rq[, "sd"])), " e ",
        sprintf("%.1f", mean(sorteios_rq[, "fora"])), " pontos fora, em media). ",
        "O envelope simulado meio-normal (99 simulacoes) deixa ", envelope$out,
        " de ", envelope$total, " pontos fora da faixa de 95%",
        if (envelope$out / envelope$total <= 0.05) ", sem afastamento sistematico: familia e ligacao compativeis com os dados. "
        else ", acima do esperado: familia ou ligacao merecem revisao. ",
        "DISPERSAO. A razao de Pearson sobre os graus de liberdade vale ",
        sprintf("%.2f", pearson_gl_final), ", reportada por exigencia do enunciado; em ",
        "dados 0/1 nao agrupados ela nao identifica superdispersao, e a verificacao ",
        "que vale e a do envelope. ",
        "CALIBRACAO. Em ", G_CALIBRACAO, " decis de probabilidade prevista, o teste de ",
        "Hosmer-Lemeshow da X2 = ", sprintf("%.2f", estat_HL), " com ", gl_HL,
        " gl (p = ", sprintf("%.3f", p_HL), ")",
        if (p_HL >= 0.05) ": sem evidencia de descalibracao. " else ": ha evidencia de descalibracao. ",
        "INFLUENCIA. A maior distancia de Cook e ", sprintf("%.4f", max(dist_cook)),
        ", pequena em valor absoluto; ", length(acima_cook), " talhoes passam da regua 4/n (",
        paste(head(dados_cafe$id_talhao[acima_cook], 5), collapse = ", "),
        if (length(acima_cook) > 5) ", ..." else "", "). Reajustando sem eles, ",
        if (trocas_sinal == 0) "nenhum coeficiente troca de sinal" else
            paste0(trocas_sinal, " coeficiente(s) troca(m) de sinal"),
        " e a maior mudanca relativa e de ", sprintf("%.0f%%", maior_mudanca),
        "; as duas versoes sao reportadas e o modelo com todos os talhoes e mantido. ",
        "LIGACAO. O termo eta^2 adicionado ao modelo tem LRT com p = ",
        sprintf("%.3f", p_eta2),
        if (p_eta2 >= 0.05) ": sem evidencia contra a ligacao logit. " else ": a ligacao logit merece revisao. ",
        "ESCALA. Para ", paste(numericas_final, collapse = " e "),
        ", o teste da variavel construida x*log(x) da p = ",
        paste(sprintf("%.3f", tab_escala$p_x_logx), collapse = " e "),
        ", e o suavizador do GAM usa edf = ",
        paste(sprintf("%.2f", tab_escala$edf_gam), collapse = " e "),
        if (!sinais[["escala_BT"]] && !sinais[["escala_GAM"]])
            ": os residuos parciais, o teste e o suavizador apontam juntos para a escala linear no logit. "
        else ": ha sinal de que alguma escala precisa ser revista; ver residuos parciais. ",
        "CONCLUSAO. ",
        if (modelo_confirmado)
            "Nenhum dos seis diagnosticos acusa problema: o modelo final e uma descricao compativel com os dados e pode ser interpretado. "
        else "Ha sinais de alerta; a revisao do modelo final precisa ser discutida e registrada antes da interpretacao. ",
        "O QUE FICA PARA OS PROXIMOS STEPS: a comparacao com e sem a severidade da ",
        "ferrugem (Step 6), as razoes de chances com intervalos (Step 7), a capacidade ",
        "preditiva (Step 8) e a recomendacao (Step 9).")


# =====================================================================
# ---- STEP 6 ---- EFEITO TOTAL E EFEITO DIRETO: A FERRUGEM NO CAMINHO?
# =====================================================================
secao("STEP 6 - EFEITO TOTAL x EFEITO DIRETO (TAREFA 2)")

# Pergunta do step: a ferrugem - que tambem e um desfecho deste projeto -
# esta no caminho entre a cultivar e a exportacao? Se sim, qual dos dois
# modelos responde a pergunta de quem ainda vai escolher a cultivar?
#
# A ordem dos acontecimentos, em linguagem simples: a cultivar e
# escolhida no plantio; a ferrugem aparece depois, ao longo da safra, e
# o Desafio B mostrou que sua severidade depende da cultivar; a
# classificacao do lote vem por ultimo. A severidade da ferrugem esta,
# portanto, NO MEIO do caminho entre a cultivar e a exportacao. Uma
# variavel nessa posicao e chamada de mediadora (ou desfecho
# intermediario). Aqui isso e uma hipotese a investigar, nao um fato.
#
# Efeito total e efeito direto (enunciado, Tarefa 2):
# - efeito TOTAL da cultivar: tudo o que trocar de cultivar faz a chance
#   de exportar, INCLUINDO o que faz por meio da ferrugem. E o que estima
#   o modelo SEM a severidade.
# - efeito DIRETO da cultivar: o que a cultivar faz "por outras vias",
#   mantendo a ferrugem fixa. E o que estima o modelo COM a severidade -
#   sob hipoteses, discutidas abaixo.
#
# Diagrama causal (DAG) adotado como hipotese:
#
#   irrigacao, densidade, umidade (Desafio B)
#                    |
#                    v
#   cultivar ---> severidade_ferrugem ---> padrao_exportacao
#       |                                        ^
#       +----------------------------------------+
#   manejo, altitude, idade da lavoura -----------+
#
#   (fator nao medido?) --> severidade  e  --> exportacao   [hipotese]
#
# As setas dizem o que se supoe, nao o que se provou. A ultima linha
# marca o risco: se existir um fator nao medido que afete a ferrugem E a
# qualidade do lote, "segurar a ferrugem fixa" mistura talhoes diferentes
# e o efeito direto fica viesado (vies de colisor; Pearl, 2009;
# VanderWeele, 2015). Esse fator nao pode ser verificado com estes dados.
#
# [APROFUNDAMENTO] Nao colapsabilidade da razao de chances.
# Na escala do logit, adicionar uma covariavel forte muda os coeficientes
# das outras MESMO que nao haja mediacao nem confusao (Greenland, Robins
# e Pearl, 1999). Logo, "o coeficiente da cultivar mudou" nao prova, por
# si so, que a ferrugem e mediadora. A leitura precisa do diagrama e da
# relacao cultivar -> severidade demonstrada no Desafio B.
#
# [APROFUNDAMENTO] "Table 2 fallacy" (Westreich e Greenland, 2013): no
# modelo com a mediadora, o coeficiente da cultivar e o da severidade
# respondem a perguntas diferentes; nao se le a tabela inteira com a
# mesma lente.

# Verificacao da mediadora antes de usa-la
severidade <- dados_cafe$severidade_ferrugem
if (any(is.na(severidade)) || any(severidade <= 0 | severidade >= 1))
    stop("severidade_ferrugem precisa estar em (0, 1) sem NA.", call. = FALSE)
cat(sprintf("severidade_ferrugem: %.3f a %.3f, sem NA  [ok]\n", min(severidade), max(severidade)))

# ---------------------------------------------------------------------
# Figura C6: o diagrama causal
# ---------------------------------------------------------------------
# Pergunta que responde: qual e a hipotese causal adotada?
# Como ler: setas cheias = caminhos supostos; setas tracejadas = o fator
# nao medido que, se existir, vicia o efeito direto. E uma hipotese,
# nao um resultado.
nos <- data.frame(
    nome = c("cultivar", "severidade da\nferrugem", "padrao de\nexportacao",
             "manejo, altitude,\nidade da lavoura", "irrigacao, densidade,\numidade (Desafio B)",
             "fator nao\nmedido?"),
    x = c(1, 3, 5, 1, 3, 5.2), y = c(2, 3.3, 2, 0.6, 4.7, 4.4),
    tipo = c("var", "var", "var", "var", "var", "hip"))
setas <- data.frame(
    x    = c(1.55, 1.55, 3.7, 1.75, 3.0, 4.75, 5.2),
    y    = c(2.25, 2.0,  3.1, 0.75, 4.35, 4.25, 4.05),
    xend = c(2.35, 4.35, 4.5, 4.35, 3.0,  3.65, 5.0),
    yend = c(3.05, 2.0,  2.25, 1.8, 3.65, 3.5,  2.35),
    tipo = c("solid", "solid", "solid", "solid", "solid", "dashed", "dashed"))
g_dag <- ggplot() +
    geom_segment(data = setas, aes(x = x, y = y, xend = xend, yend = yend, linetype = tipo),
                 arrow = arrow(length = unit(.22, "cm"), type = "closed"), colour = "black") +
    geom_label(data = nos, aes(x = x, y = y, label = nome, fill = tipo), size = 3.6,
               label.padding = unit(.35, "lines"), colour = "black") +
    scale_fill_manual(values = c(var = "white", hip = "#F3E1DF"), guide = "none") +
    scale_linetype_identity() +
    coord_cartesian(xlim = c(0.2, 6), ylim = c(0, 5.4)) +
    labs(title = "Hipotese causal: a severidade da ferrugem esta no caminho entre cultivar e exportacao",
         subtitle = "setas cheias: caminhos supostos; tracejadas: fator nao medido que viciaria o efeito direto",
         caption = "diagrama de hipoteses; nenhuma seta e provada por estes dados") +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), plot.margin = margin(8, 8, 8, 8))
salvar(g_dag, "C6_dag.png", 10, 6)

# ---------------------------------------------------------------------
# Os dois modelos
# ---------------------------------------------------------------------
# OBJETIVO: estimar o efeito total (sem a severidade) e o efeito direto
#   (com a severidade) da cultivar.
# O QUE O CODIGO FAZ: modelo_total e o modelo final do Step 4;
#   modelo_direto acrescenta severidade_ferrugem; verificar() confere
#   convergencia; um envelope rapido confere que o modelo direto tambem
#   e compativel com os dados; anova() testa se a severidade acrescenta
#   ajuste.
# COMO INTERPRETAR: se a severidade tem coeficiente negativo e
#   significativo, mais ferrugem = menor chance de exportar, ajustado
#   pelas demais caracteristicas. O sinal do coeficiente sera usado para
#   conferir a coerencia da decomposicao a seguir.
modelo_total  <- modelo_final
modelo_direto <- update(modelo_final, . ~ . + severidade_ferrugem)
verificar(modelo_direto, "modelo_direto")
envelope_direto <- hnp(modelo_direto, resid.type = "deviance", sim = 99, conf = 0.95,
                       how.many.out = TRUE, print.on = FALSE, plot.sim = FALSE)
cat(sprintf("   envelope do modelo direto: %d de %d pontos fora\n", envelope_direto$out, envelope_direto$total))
lrt_sev <- anova(modelo_total, modelo_direto, test = "Chisq")
beta_sev <- coef(modelo_direto)[["severidade_ferrugem"]]
cat(sprintf("   severidade_ferrugem no modelo direto: beta = %.3f | LRT = %.2f com 1 gl, p = %s\n",
            beta_sev, lrt_sev$Deviance[2], format.pval(lrt_sev[["Pr(>Chi)"]][2], digits = 3)))

# ---------------------------------------------------------------------
# Comparacao dos coeficientes de cultivar
# ---------------------------------------------------------------------
# OBJETIVO: ver como o coeficiente de cada cultivar muda ao incluir a
#   severidade, e se a mudanca e coerente com a hipotese de mediacao.
# O QUE O CODIGO FAZ: extrai os coeficientes e razoes de chances de
#   cultivar nos dois modelos; calcula a diferenca (total - direto), que
#   sob a hipotese do diagrama e a parte do efeito que passa pela
#   ferrugem; calcula a severidade media de cada cultivar e a diferenca
#   em relacao a Bourbon.
# COMO INTERPRETAR: coerencia esperada sob mediacao - uma cultivar com
#   MENOS ferrugem que Bourbon (diferenca negativa), combinada com um
#   coeficiente NEGATIVO da severidade, deve ter parte indireta POSITIVA
#   (a menor ferrugem ajuda a exportar). Nesse caso, o efeito total fica
#   MENOS negativo que o direto. A coluna "coerente" confere isso.
#   Coerencia nao e prova: a nao colapsabilidade tambem move o
#   coeficiente na mesma direcao.
or_total  <- tabela_or(modelo_total)
or_direto <- tabela_or(modelo_direto)
linhas_cultivar <- grep("^cultivar", or_total$termo)
niveis_cultivar <- sub("^cultivar", "", or_total$termo[linhas_cultivar])
referencia_cultivar <- levels(dados_cafe$cultivar)[1]

sev_media <- tapply(severidade, dados_cafe$cultivar, mean)
tab_total_direto <- data.frame(
    cultivar     = niveis_cultivar,
    beta_total   = or_total$beta[linhas_cultivar],
    beta_direto  = or_direto$beta[match(or_total$termo[linhas_cultivar], or_direto$termo)],
    OR_total     = or_total$OR[linhas_cultivar],
    OR_direto    = or_direto$OR[match(or_total$termo[linhas_cultivar], or_direto$termo)],
    sev_media    = as.vector(sev_media[niveis_cultivar]),
    sev_dif_ref  = as.vector(sev_media[niveis_cultivar]) - sev_media[[referencia_cultivar]])
tab_total_direto$indireto  <- tab_total_direto$beta_total - tab_total_direto$beta_direto
tab_total_direto$coerente  <- sign(tab_total_direto$indireto) == sign(beta_sev * tab_total_direto$sev_dif_ref)

cat(sprintf("\ncoeficientes de cultivar (referencia: %s; severidade media da referencia: %.3f):\n",
            referencia_cultivar, sev_media[[referencia_cultivar]]))
print(data.frame(cultivar    = tab_total_direto$cultivar,
                 beta_total  = round(tab_total_direto$beta_total, 3),
                 beta_direto = round(tab_total_direto$beta_direto, 3),
                 indireto    = round(tab_total_direto$indireto, 3),
                 OR_total    = round(tab_total_direto$OR_total, 3),
                 OR_direto   = round(tab_total_direto$OR_direto, 3),
                 sev_media   = round(tab_total_direto$sev_media, 3),
                 sev_vs_ref  = round(tab_total_direto$sev_dif_ref, 3),
                 coerente    = ifelse(tab_total_direto$coerente, "sim", "NAO")),
      row.names = FALSE)
cat("   indireto = beta_total - beta_direto (parte que passaria pela ferrugem, sob a hipotese do diagrama)\n")

# A mesma comparacao no modelo amplo: o fenomeno nao e artefato da selecao
modelo_amplo_direto <- update(modelo_amplo, . ~ . + severidade_ferrugem)
b_amplo_total  <- coef(modelo_amplo)[or_total$termo[linhas_cultivar]]
b_amplo_direto <- coef(modelo_amplo_direto)[or_total$termo[linhas_cultivar]]
cat("\n   no modelo amplo (13 candidatas), a mesma comparacao:\n")
print(data.frame(cultivar = niveis_cultivar,
                 beta_total_amplo  = round(unname(b_amplo_total), 3),
                 beta_direto_amplo = round(unname(b_amplo_direto), 3),
                 indireto_amplo    = round(unname(b_amplo_total - b_amplo_direto), 3)),
      row.names = FALSE)

# Razao de chances da propria severidade, por 0,10 de severidade.
# A severidade varia de 0,001 a 0,61; "uma unidade" (0 -> 1) nao existe
# na base. Um incremento de 0,10 (10 pontos percentuais de area foliar
# lesionada) e uma variacao concreta: exp(0,10 * beta).
ic_sev <- suppressMessages(confint(modelo_direto))["severidade_ferrugem", ]
or_sev_010 <- exp(0.10 * c(beta_sev, ic_sev))
cat(sprintf("\n   razao de chances por +0,10 de severidade: %.3f [IC95%%: %.3f; %.3f]\n",
            or_sev_010[1], or_sev_010[2], or_sev_010[3]))

# ---------------------------------------------------------------------
# Figura C6: razoes de chances de cultivar nos dois modelos
# ---------------------------------------------------------------------
# Pergunta que responde: como a razao de chances de cada cultivar muda
# ao incluir a mediadora?
# Como ler: eixo x em escala log; linha em 1 = igual a Bourbon; para
# cada cultivar, o ponto do modelo total (sem severidade) e o do direto
# (com severidade), com IC de perfil. O deslocamento entre os dois e o
# efeito de incluir a mediadora.
# Limitacao: deslocamento nao prova mediacao (nao colapsabilidade).
df_or_cult <- rbind(
    cbind(or_total[linhas_cultivar, ],  modelo = "total (sem severidade)"),
    cbind(or_direto[match(or_total$termo[linhas_cultivar], or_direto$termo), ], modelo = "direto (com severidade)"))
df_or_cult$cultivar <- sub("^cultivar", "", df_or_cult$termo)
df_or_cult$modelo   <- factor(df_or_cult$modelo, levels = c("total (sem severidade)", "direto (com severidade)"))
g_total_direto <- ggplot(df_or_cult, aes(x = OR, y = cultivar, colour = modelo)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = CINZA) +
    geom_errorbar(aes(xmin = IC_inf, xmax = IC_sup), width = .25,
                  position = position_dodge(width = .6)) +
    geom_point(size = 3, position = position_dodge(width = .6)) +
    scale_x_log10() +
    scale_colour_manual(values = c(`total (sem severidade)` = VERDE, `direto (com severidade)` = VERMELHO), name = NULL) +
    labs(x = sprintf("razao de chances em relacao a %s (escala log)", referencia_cultivar), y = NULL,
         title = "Cultivar: efeito total e efeito direto na chance de exportar",
         subtitle = "o efeito direto mantem a severidade da ferrugem fixa; o total nao",
         caption = "IC de verossimilhanca perfilada (95%); o deslocamento entre os pontos nao prova mediacao") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
salvar(g_total_direto, "C6_total_vs_direto.png", 9, 5)

# ---------------------------------------------------------------------
# Parecer do Step 6
# ---------------------------------------------------------------------
# CONCLUSAO LITERAL:
# O que muda no coeficiente da cultivar e o que a ferrugem faz.
#
# CONCLUSAO ESTATISTICA:
# Total x direto, coerencia com a hipotese, e o que nao esta provado.
#
# RESPOSTA A TAREFA 2:
# Qual modelo serve a quem escolhe a cultivar antes do plantio.
#
# LIMITACOES:
# Hipoteses do diagrama; nao colapsabilidade; fator nao medido.
cultivar_menor_sev <- names(which.min(sev_media))
cultivar_menor_or  <- tab_total_direto$cultivar[which.min(tab_total_direto$OR_total)]
n_coerentes <- sum(tab_total_direto$coerente)

parecer("CONCLUSAO LITERAL. Mais ferrugem, menor chance de exportar: cada 0,10 a mais ",
        "de severidade multiplica a chance por ", sprintf("%.2f", or_sev_010[1]),
        " (IC95% ", sprintf("%.2f a %.2f", or_sev_010[2], or_sev_010[3]),
        "), ajustado pelas demais caracteristicas. Ao incluir a severidade no modelo, ",
        "os coeficientes de cultivar mudam: ",
        paste(sprintf("%s de %.2f para %.2f", tab_total_direto$cultivar,
                      tab_total_direto$beta_total, tab_total_direto$beta_direto), collapse = "; "),
        " (escala do logit, em relacao a ", referencia_cultivar, "). A cultivar com menor ",
        "severidade media e ", cultivar_menor_sev, " (", sprintf("%.3f", sev_media[[cultivar_menor_sev]]),
        " contra ", sprintf("%.3f", sev_media[[referencia_cultivar]]), " de ", referencia_cultivar,
        "), e a com menor razao de chances total de exportar e ", cultivar_menor_or, ". ",
        "CONCLUSAO ESTATISTICA. O modelo sem a severidade estima o efeito total da ",
        "cultivar; o modelo com a severidade estima o efeito direto, sob as hipoteses ",
        "do diagrama. Em ", n_coerentes, " das ", nrow(tab_total_direto),
        " cultivares, o sinal da diferenca total - direto e coerente com a mediacao ",
        "(cultivares com menos ferrugem que ", referencia_cultivar, " tem parte indireta ",
        "positiva, o que torna o efeito total menos negativo que o direto). O mesmo ",
        "padrao aparece no modelo amplo, logo nao e artefato da selecao. A coerencia, ",
        "porem, nao prova mediacao: a razao de chances nao e colapsavel, e parte do ",
        "deslocamento ocorreria mesmo sem mediacao. A evidencia de mediacao apoia-se ",
        "no diagrama e na relacao cultivar -> severidade demonstrada no Desafio B. ",
        "RESPOSTA A TAREFA 2. O produtor que escolhe a cultivar antes do plantio ",
        "recebe a ferrugem junto com a cultivar: ele nao consegue 'segurar a ferrugem ",
        "fixa'. A pergunta dele e respondida pelo efeito TOTAL, isto e, pelo modelo SEM ",
        "a severidade. O efeito direto responde a outra pergunta - 'e se a ferrugem ",
        "fosse controlada ao mesmo nivel em todas as cultivares?' -, util para quem ",
        "planeja o programa de controle da doenca, nao para quem escolhe o que plantar. ",
        "LIMITACOES. O efeito direto so e interpretavel se nao houver fator nao medido ",
        "que afete a ferrugem e a qualidade do lote ao mesmo tempo; com estes dados ",
        "isso e uma hipotese, nao uma verificacao. ",
        "O QUE FICA PARA OS PROXIMOS STEPS: as magnitudes com intervalos (Step 7) e a ",
        "reconciliacao com o Desafio B na recomendacao (Step 9).")


# =====================================================================
# ---- STEP 7 ---- RAZOES DE CHANCES COM INTERVALOS DE CONFIANCA
# =====================================================================
secao("STEP 7 - RAZOES DE CHANCES (TAREFA 4)")

# Pergunta do step: em quanto cada caracteristica multiplica a chance de
# um lote ser exportavel - e com que incerteza?
#
# Chance e probabilidade, de novo (Step 0): probabilidade e "de cada 100
# lotes, quantos exportam"; chance e "quantos exportam para cada um que
# nao exporta", p / (1 - p).
#
# Razao de chances (slides, secao 9.1.2), em linguagem simples: por
# quanto a chance e multiplicada ao trocar de nivel (fatores) ou ao
# subir uma unidade (numericas), mantendo o resto igual. Formula:
# OR = exp(beta). E constante em toda a escala: o mesmo multiplicador
# vale para um talhao que parte de chance 0,5 ou de chance 2.
#
# Aviso dos slides: RAZAO DE CHANCES NAO E RISCO RELATIVO. OR = 2 nao
# significa "o dobro da probabilidade". A aproximacao OR ~ RR so vale
# para desfechos raros (p < 0,10); aqui o desfecho e frequente (cerca de
# metade dos lotes), e a diferenca e grande. O bloco final deste step
# mostra isso com numeros.
#
# Intervalo de confianca (apostila, secao 2.11): faixa de valores
# compativeis com os dados. Se contem 1, os dados nao excluem "sem
# efeito". Dois jeitos de calcula-lo: Wald (beta +- 1,96 EP, simetrico
# no logit) e perfil de verossimilhanca (confint() no glm; assimetrico
# na escala da razao de chances e mais confiavel em amostras moderadas).
#
# OBJETIVO: interpretar o modelo final na escala das razoes de chances.
# O QUE O CODIGO FAZ: tabela_or() extrai beta, erro-padrao, OR e IC de
#   perfil; um segundo bloco reescala as numericas para unidades
#   agronomicas (100 m de altitude; 5 anos de lavoura) por exp(c * beta);
#   um terceiro compara Wald e perfil no termo mais incerto.
# COMO INTERPRETAR: OR > 1 = chance maior que a referencia (ou que
#   subir uma unidade aumenta a chance); OR < 1 = menor. IC que cruza 1 =
#   sem evidencia clara. A magnitude em probabilidade depende do ponto
#   de partida (Step 9).
tab_or_final <- tabela_or(modelo_final)
tab_or_final$p_Wald <- summary(modelo_final)$coefficients[, "Pr(>|z|)"]

cat(sprintf("razoes de chances do modelo final (referencias: %s; %s):\n",
            paste(sapply(fatores_final, function(f) sprintf("%s = %s", f, levels(dados_cafe[[f]])[1])), collapse = "; "),
            "IC de verossimilhanca perfilada, 95%"))
print(data.frame(termo  = tab_or_final$termo,
                 beta   = round(tab_or_final$beta, 4),
                 OR     = round(tab_or_final$OR, 3),
                 IC_inf = round(tab_or_final$IC_inf, 3),
                 IC_sup = round(tab_or_final$IC_sup, 3),
                 p_Wald = round(tab_or_final$p_Wald, 4),
                 cruza_1 = ifelse(tab_or_final$IC_inf < 1 & tab_or_final$IC_sup > 1, "sim", "")),
      row.names = FALSE)

# ---------------------------------------------------------------------
# Reescalar as numericas para unidades com sentido agronomico
# ---------------------------------------------------------------------
# OR por 1 metro de altitude e um numero perto de 1 e dificil de ler.
# OR por 100 m = exp(100 * beta) - e NAO 100 * exp(beta). O mesmo vale
# para o IC: exp(100 * limite do IC de beta). Escalas adotadas:
ESCALAS <- c(altitude_m = 100, idade_lavoura_anos = 5)
tab_or_reescalada <- bind_rows(lapply(numericas_final, function(v) {
    c_esc <- if (v %in% names(ESCALAS)) ESCALAS[[v]] else 1
    linha <- tab_or_final[tab_or_final$termo == v, ]
    ic_b  <- suppressMessages(confint(modelo_final))[v, ]
    data.frame(termo = v, unidade = c_esc,
               OR_por_unidade = exp(c_esc * linha$beta),
               IC_inf = exp(c_esc * ic_b[1]), IC_sup = exp(c_esc * ic_b[2]))
}))
cat("\nnumericas reescaladas (OR por c unidades = exp(c * beta)):\n")
print(data.frame(termo   = tab_or_reescalada$termo,
                 por     = tab_or_reescalada$unidade,
                 OR      = round(tab_or_reescalada$OR_por_unidade, 3),
                 IC_inf  = round(tab_or_reescalada$IC_inf, 3),
                 IC_sup  = round(tab_or_reescalada$IC_sup, 3)),
      row.names = FALSE)

# ---------------------------------------------------------------------
# Wald e perfil: a diferenca no termo mais incerto
# ---------------------------------------------------------------------
# [APROFUNDAMENTO] O IC de Wald supoe que a log-verossimilhanca e uma
# parabola em torno da estimativa; o de perfil percorre a
# log-verossimilhanca de verdade. Quando o erro-padrao e grande, a
# parabola e uma aproximacao pior e os dois intervalos se afastam. Na
# escala da razao de chances, o de Wald e simetrico em torno de beta
# (assimetrico em torno de OR) e o de perfil nao tem essa restricao.
termo_mais_incerto <- tab_or_final$termo[-1][which.max(tab_or_final$ep[-1])]
linha_inc <- tab_or_final[tab_or_final$termo == termo_mais_incerto, ]
ic_wald   <- exp(linha_inc$beta + c(-1, 1) * qnorm(0.975) * linha_inc$ep)
cat(sprintf("\ntermo com maior erro-padrao: %s (EP = %.3f)\n", termo_mais_incerto, linha_inc$ep))
cat(sprintf("   IC95%% de Wald   : [%.3f; %.3f]\n", ic_wald[1], ic_wald[2]))
cat(sprintf("   IC95%% de perfil : [%.3f; %.3f]\n", linha_inc$IC_inf, linha_inc$IC_sup))

# ---------------------------------------------------------------------
# Figura C7: forest plot das razoes de chances
# ---------------------------------------------------------------------
# Pergunta que responde: quanto cada caracteristica multiplica a chance?
# Como ler: eixo x em escala log (distancias iguais = mesmo fator
# multiplicativo); linha em 1 = sem efeito; ponto = OR; barra = IC de
# perfil 95%. Fatores sao lidos em relacao a categoria de referencia,
# escrita no rotulo. Numericas aparecem na unidade reescalada.
# Limitacao: OR nao e risco relativo; a leitura em probabilidade vem
# no Step 9.
rotulo_termo <- function(t) {
    for (f in fatores_final) if (startsWith(t, f))
        return(sprintf("%s: %s (ref.: %s)", f, sub(f, "", t), levels(dados_cafe[[f]])[1]))
    if (t %in% tab_or_reescalada$termo)
        return(sprintf("%s (por %g)", t, tab_or_reescalada$unidade[tab_or_reescalada$termo == t]))
    t
}
df_forest <- tab_or_final[tab_or_final$termo != "(Intercept)", c("termo", "OR", "IC_inf", "IC_sup")]
for (v in tab_or_reescalada$termo) {
    i <- df_forest$termo == v
    df_forest$OR[i]     <- tab_or_reescalada$OR_por_unidade[tab_or_reescalada$termo == v]
    df_forest$IC_inf[i] <- tab_or_reescalada$IC_inf[tab_or_reescalada$termo == v]
    df_forest$IC_sup[i] <- tab_or_reescalada$IC_sup[tab_or_reescalada$termo == v]
}
df_forest$rotulo <- sapply(df_forest$termo, rotulo_termo)
df_forest$rotulo <- factor(df_forest$rotulo, levels = rev(df_forest$rotulo))
df_forest$evidencia <- ifelse(df_forest$IC_inf > 1 | df_forest$IC_sup < 1, "IC exclui 1", "IC cruza 1")
g_forest <- ggplot(df_forest, aes(x = OR, y = rotulo, colour = evidencia)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = CINZA) +
    geom_errorbar(aes(xmin = IC_inf, xmax = IC_sup), width = .25) +
    geom_point(size = 3) +
    geom_text(aes(label = sprintf("%.2f [%.2f; %.2f]", OR, IC_inf, IC_sup)),
              vjust = -.9, size = 3, colour = "black") +
    scale_x_log10() +
    scale_colour_manual(values = c(`IC exclui 1` = VERDE, `IC cruza 1` = CINZA), name = NULL) +
    labs(x = "razao de chances (escala log)", y = NULL,
         title = "Razoes de chances do modelo final, com IC de perfil 95%",
         subtitle = "linha tracejada: OR = 1 (sem efeito); numericas em unidades reescaladas",
         caption = "razao de chances nao e risco relativo: a leitura em probabilidade depende do ponto de partida") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
salvar(g_forest, "C7_forest_or.png", 10, 6)

# ---------------------------------------------------------------------
# Razao de chances nao e risco relativo: um exemplo com numeros
# ---------------------------------------------------------------------
# OBJETIVO: mostrar quanto a mesma razao de chances vale em pontos
#   percentuais, partindo de probabilidades diferentes.
# O QUE O CODIGO FAZ: toma a razao de chances do manejo com maior OR e
#   aplica-a a tres probabilidades de partida: a de um talhao de
#   referencia previsto pelo modelo, e duas hipoteticas (0,10 e 0,50).
#   chance1 = OR * chance0; p1 = chance1 / (1 + chance1).
# COMO INTERPRETAR: a coluna "risco relativo" (p1/p0) muda com o ponto de
#   partida, enquanto a OR e a mesma; a diferenca em pontos percentuais
#   e maior perto de p = 0,5.
termo_manejo_max <- df_forest$termo[startsWith(df_forest$termo, "manejo")][
    which.max(df_forest$OR[startsWith(df_forest$termo, "manejo")])]
or_exemplo <- tab_or_final$OR[tab_or_final$termo == termo_manejo_max]
talhao_base <- dados_cafe[1, termos_final]
for (v in numericas_final) talhao_base[[v]] <- median(dados_cafe[[v]])
for (f in fatores_final) talhao_base[[f]] <- factor(levels(dados_cafe[[f]])[1], levels = levels(dados_cafe[[f]]))
p_base <- predict(modelo_final, newdata = talhao_base, type = "response")
p0 <- c(referencia = unname(p_base), hipotetico_010 = 0.10, hipotetico_080 = 0.80)
chance0 <- p0 / (1 - p0)
p1 <- (or_exemplo * chance0) / (1 + or_exemplo * chance0)
tab_or_rr <- data.frame(ponto_de_partida = names(p0), p0 = round(p0, 3), p1 = round(p1, 3),
                        diferenca_pp = round(100 * (p1 - p0), 1),
                        risco_relativo = round(p1 / p0, 2), razao_de_chances = round(or_exemplo, 2))
cat(sprintf("\na mesma razao de chances (%s, OR = %.2f) em probabilidades diferentes:\n",
            sub("manejo", "manejo ", termo_manejo_max), or_exemplo))
print(tab_or_rr, row.names = FALSE)

# ---------------------------------------------------------------------
# Parecer do Step 7
# ---------------------------------------------------------------------
# CONCLUSAO LITERAL: cada caracteristica em "multiplica a chance por".
# CONCLUSAO ESTATISTICA: OR com IC de perfil; onde o IC cruza 1.
# RESPOSTA A TAREFA 4: interpretacao em razoes de chances com IC.
# O QUE FICA: a leitura em probabilidade (Step 9).
frase_or <- function(i) {
    r <- df_forest[i, ]
    sprintf("%s: OR %.2f [%.2f; %.2f]", as.character(r$rotulo), r$OR, r$IC_inf, r$IC_sup)
}
com_evidencia <- which(df_forest$evidencia == "IC exclui 1")
sem_evidencia <- which(df_forest$evidencia == "IC cruza 1")

parecer("CONCLUSAO LITERAL. Mantidas as demais caracteristicas, ",
        paste(sapply(seq_len(nrow(df_forest)), frase_or), collapse = "; "),
        ". Lendo em linguagem comum: uma razao de chances de 2 significa que, para cada ",
        "lote que nao exporta, o grupo tem o dobro de lotes exportando em relacao a ",
        "referencia - nao 'o dobro da probabilidade'. O exemplo numerico mostra que a ",
        "mesma OR de ", sprintf("%.2f", or_exemplo), " vale ",
        sprintf("%.1f", tab_or_rr$diferenca_pp[1]), " pontos percentuais partindo de ",
        sprintf("%.2f", tab_or_rr$p0[1]), ", mas ", sprintf("%.1f", tab_or_rr$diferenca_pp[2]),
        " partindo de 0,10. ",
        "CONCLUSAO ESTATISTICA. Os intervalos sao de verossimilhanca perfilada; no termo ",
        "mais incerto (", termo_mais_incerto, ") ele vai de ",
        sprintf("%.2f a %.2f", linha_inc$IC_inf, linha_inc$IC_sup), " contra ",
        sprintf("%.2f a %.2f", ic_wald[1], ic_wald[2]), " pelo metodo de Wald. ",
        "Os termos cujo IC exclui 1 sao: ",
        if (length(com_evidencia) > 0) paste(as.character(df_forest$rotulo[com_evidencia]), collapse = "; ") else "nenhum",
        ". Os termos cujo IC cruza 1 sao: ",
        if (length(sem_evidencia) > 0) paste(as.character(df_forest$rotulo[sem_evidencia]), collapse = "; ") else "nenhum",
        " - para eles os dados nao excluem ausencia de efeito, mas o intervalo mostra ",
        "a faixa de magnitudes ainda plausivel, o que nao e o mesmo que 'nao ha efeito'. ",
        "RESPOSTA A TAREFA 4. Os efeitos estao interpretados em razoes de chances com ",
        "intervalos de confianca, em relacao as categorias de referencia (",
        paste(sapply(fatores_final, function(f) levels(dados_cafe[[f]])[1]), collapse = ", "),
        ") e em unidades agronomicas para as numericas. ",
        "O QUE FICA PARA OS PROXIMOS STEPS: a capacidade preditiva (Step 8) e a ",
        "traducao das razoes de chances em probabilidades para cenarios concretos (Step 9).")


# =====================================================================
# ---- STEP 8 ---- CAPACIDADE PREDITIVA E PONTO DE CORTE
# =====================================================================
secao("STEP 8 - CAPACIDADE PREDITIVA (TAREFA 3)")

# Pergunta do step: se usassemos o modelo para apontar, antes da
# classificacao oficial, quais talhoes tendem a produzir lote
# exportavel, quantas vezes acertariamos - e onde colocar a regua,
# sabendo que errar para um lado custa mais do que para o outro?
#
# De probabilidade a decisao, em linguagem simples: o modelo da uma
# probabilidade a cada talhao; para agir e preciso uma regua (ponto de
# corte). Acima dela, o talhao e tratado como "exportavel". A regua
# 0,5 e apenas o padrao ingenuo - ela so e otima quando os dois tipos de
# erro custam o mesmo.
#
# Matriz de confusao: tabela 2 x 2 de acertos e erros.
#   verdadeiro positivo (VP): previsto exportavel, e exportou;
#   falso positivo (FP)     : previsto exportavel, e nao exportou;
#   falso negativo (FN)     : previsto nao exportavel, e exportou;
#   verdadeiro negativo (VN): previsto nao exportavel, e nao exportou.
# Sensibilidade = VP / (VP + FN): dos que exportaram, quantos o modelo
# apontou. Especificidade = VN / (VN + FP): dos que nao exportaram,
# quantos ele deixou de fora. Valor preditivo positivo = VP / (VP + FP):
# dos apontados, quantos de fato exportaram.
# (Hosmer, Lemeshow e Sturdivant, 2013; Fawcett, 2006.)
#
# OBJETIVO: medir o acerto do modelo com a regua ingenua, como ponto de
#   partida.
# O QUE O CODIGO FAZ: a funcao metricas() classifica com um corte e
#   devolve a matriz e as medidas.
# COMO INTERPRETAR: a acuracia sozinha engana quando as classes sao
#   desiguais; sensibilidade e especificidade dizem em que direcao o
#   modelo erra. E aqui que o possivel desbalanceamento das classes
#   (Step 0) e retomado: as contagens de cada classe sao impressas ao
#   lado das medidas, junto com a acuracia de quem sempre apostasse na
#   classe maior - a referencia minima que qualquer modelo deve superar.
prob_predita <- fitted(modelo_final)

metricas <- function(prob, obs, corte) {
    pred <- as.integer(prob >= corte)
    VP <- sum(pred == 1 & obs == 1); FP <- sum(pred == 1 & obs == 0)
    FN <- sum(pred == 0 & obs == 1); VN <- sum(pred == 0 & obs == 0)
    c(corte = corte, VP = VP, FP = FP, FN = FN, VN = VN,
      sensibilidade = VP / (VP + FN), especificidade = VN / (VN + FP),
      acuracia = (VP + VN) / length(obs),
      VPP = if (VP + FP > 0) VP / (VP + FP) else NA,
      VPN = if (VN + FN > 0) VN / (VN + FN) else NA)
}

m_05 <- metricas(prob_predita, resposta_exportacao, 0.5)
cat(sprintf("regua ingenua (corte 0,5): VP = %d | FP = %d | FN = %d | VN = %d\n",
            m_05[["VP"]], m_05[["FP"]], m_05[["FN"]], m_05[["VN"]]))
cat(sprintf("   sensibilidade = %.3f | especificidade = %.3f | acuracia = %.3f | VPP = %.3f | VPN = %.3f\n",
            m_05[["sensibilidade"]], m_05[["especificidade"]], m_05[["acuracia"]], m_05[["VPP"]], m_05[["VPN"]]))
cat(sprintf("   classes: %d com padrao (%.1f%%) e %d sem; acuracia de sempre prever a classe maior: %.3f\n",
            frequencia_exportacao[["1"]], 100 * proporcao_exportacao, frequencia_exportacao[["0"]],
            max(proporcao_exportacao, 1 - proporcao_exportacao)))

# Figura C8: matriz de confusao no corte ingenuo
# Pergunta que responde: onde o modelo acerta e onde erra, com a regua 0,5?
# Como ler: linhas = previsto; colunas = observado; diagonal = acertos.
df_matriz <- data.frame(previsto  = factor(c("exportavel", "exportavel", "nao exportavel", "nao exportavel"),
                                           levels = c("nao exportavel", "exportavel")),
                        observado = factor(c("atingiu (1)", "nao atingiu (0)", "atingiu (1)", "nao atingiu (0)"),
                                           levels = c("nao atingiu (0)", "atingiu (1)")),
                        n = c(m_05[["VP"]], m_05[["FP"]], m_05[["FN"]], m_05[["VN"]]),
                        tipo = c("VP", "FP", "FN", "VN"))
df_matriz$acerto <- df_matriz$tipo %in% c("VP", "VN")
g_matriz <- ggplot(df_matriz, aes(x = observado, y = previsto, fill = acerto)) +
    geom_tile(colour = "black") +
    geom_text(aes(label = sprintf("%s\n%d", tipo, n)), size = 5) +
    scale_fill_manual(values = c(`TRUE` = "#CFE3DE", `FALSE` = "#F3E1DF"), guide = "none") +
    labs(x = "observado", y = "previsto pelo modelo (corte 0,5)",
         title = "Matriz de confusao com a regua ingenua (0,5)",
         subtitle = sprintf("sensibilidade %.2f | especificidade %.2f | acuracia %.2f",
                            m_05[["sensibilidade"]], m_05[["especificidade"]], m_05[["acuracia"]]),
         caption = "0,5 so e a regua otima quando falso positivo e falso negativo custam o mesmo") +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank())
salvar(g_matriz, "C8_matriz_confusao.png", 7, 6)

# ---------------------------------------------------------------------
# Curva ROC e AUC
# ---------------------------------------------------------------------
# Em linguagem simples: em vez de escolher uma regua, experimentar todas.
# Para cada corte possivel, anota-se a sensibilidade e a especificidade;
# a curva ROC desenha uma contra a outra. A area sob a curva (AUC) e a
# probabilidade de o modelo dar nota maior a um talhao que exportou do
# que a um que nao exportou: 0,5 e uma moeda; 1 e separacao perfeita.
# O intervalo de confianca da AUC e o de DeLong (pROC::ci.auc).
#
# OBJETIVO: medir a capacidade de ordenar do modelo, independente da regua.
# O QUE O CODIGO FAZ: pROC::roc() varre os cortes; auc() e ci.auc()
#   calculam a area e o IC; coords() devolve sensibilidade e
#   especificidade por corte, e o corte de Youden (maximo de
#   sensibilidade + especificidade - 1), criterio que ignora custos.
# COMO INTERPRETAR: AUC perto de 0,5 = ordena mal; perto de 1 = ordena
#   bem. AUC nao diz nada sobre calibracao (Step 5): sao coisas
#   diferentes.
curva_roc  <- pROC::roc(response = resposta_exportacao, predictor = prob_predita,
                        levels = c(0, 1), direction = "<", quiet = TRUE)
auc_amostra <- as.numeric(pROC::auc(curva_roc))
ic_auc      <- as.numeric(pROC::ci.auc(curva_roc))
corte_youden <- pROC::coords(curva_roc, "best", best.method = "youden",
                             ret = c("threshold", "sensitivity", "specificity"))
cat(sprintf("\nAUC na propria amostra: %.3f [IC95%% DeLong: %.3f; %.3f]\n", auc_amostra, ic_auc[1], ic_auc[3]))
cat(sprintf("   corte de Youden: %.3f (sensibilidade %.3f, especificidade %.3f) - criterio sem custos\n",
            corte_youden$threshold, corte_youden$sensitivity, corte_youden$specificity))

# ---------------------------------------------------------------------
# Otimismo e validacao cruzada
# ---------------------------------------------------------------------
# Em linguagem simples: medir o acerto nos mesmos talhoes que ajustaram
# o modelo superestima, porque o modelo "ja viu" as respostas. Na
# validacao cruzada em 10 partes, o modelo e ajustado em 9 partes e
# avaliado na decima, em rodizio; cada talhao recebe uma probabilidade
# prevista por um modelo que nao o usou. A AUC dessas probabilidades e
# uma estimativa honesta (James et al., An Introduction to Statistical
# Learning, cap. 5). Depende do sorteio das partes: semente fixa.
K_FOLDS <- 10
set.seed(20260920)
fold <- sample(rep(seq_len(K_FOLDS), length.out = n_obs))
prob_cv <- numeric(n_obs)
for (k in seq_len(K_FOLDS)) {
    m_k <- glm(formula(modelo_final), family = binomial(link = "logit"), data = dados_cafe[fold != k, ])
    prob_cv[fold == k] <- predict(m_k, newdata = dados_cafe[fold == k, ], type = "response")
}
auc_cv <- as.numeric(pROC::auc(pROC::roc(resposta_exportacao, prob_cv, levels = c(0, 1), direction = "<", quiet = TRUE)))
cat(sprintf("   AUC por validacao cruzada (%d partes): %.3f | otimismo (amostra - CV): %.3f\n",
            K_FOLDS, auc_cv, auc_amostra - auc_cv))

# Referencia: AUC na amostra dos outros modelos ja ajustados
auc_de <- function(m) as.numeric(pROC::auc(pROC::roc(resposta_exportacao, fitted(m), levels = c(0, 1), direction = "<", quiet = TRUE)))
cat(sprintf("   referencia (na amostra): modelo amplo %.3f | modelo direto (com severidade) %.3f | modelo nulo %.3f\n",
            auc_de(modelo_amplo), auc_de(modelo_direto), 0.5))

# Figura C8: curva ROC
# Pergunta que responde: quao bem o modelo ordena os talhoes?
# Como ler: eixo x = 1 - especificidade (falsos positivos); eixo y =
# sensibilidade; a diagonal e a moeda; quanto mais a curva se afasta
# dela para o canto superior esquerdo, melhor. Os pontos marcam a regua
# 0,5 e o corte de Youden.
df_roc <- pROC::coords(curva_roc, "all", ret = c("threshold", "sensitivity", "specificity"))
pontos_roc <- rbind(
    data.frame(regua = "corte 0,5", sensibilidade = m_05[["sensibilidade"]], especificidade = m_05[["especificidade"]]),
    data.frame(regua = "Youden", sensibilidade = corte_youden$sensitivity, especificidade = corte_youden$specificity))
g_roc <- ggplot(df_roc, aes(x = 1 - specificity, y = sensitivity)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = CINZA) +
    geom_path(colour = VERDE, linewidth = 1) +
    geom_point(data = pontos_roc, aes(x = 1 - especificidade, y = sensibilidade, shape = regua),
               size = 3.5, colour = VERMELHO) +
    coord_equal() +
    labs(x = "1 - especificidade (fracao de falsos positivos)", y = "sensibilidade",
         shape = NULL,
         title = "Curva ROC do modelo final",
         subtitle = sprintf("AUC = %.3f [%.3f; %.3f] na amostra; %.3f por validacao cruzada",
                            auc_amostra, ic_auc[1], ic_auc[3], auc_cv),
         caption = "diagonal = classificacao ao acaso; AUC mede ordenacao, nao calibracao") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
salvar(g_roc, "C8_roc.png", 7, 7.5)

# ---------------------------------------------------------------------
# O ponto de corte como decisao de custo
# ---------------------------------------------------------------------
# O que e cada erro, para a cooperativa (definicao adotada neste trabalho):
#   FALSO POSITIVO: tratar como exportavel um lote que nao atinge o
#     padrao - o lote entra na cadeia de exportacao (contrato, logistica,
#     certificacao) e e rejeitado ou rebaixado; custo de retrabalho,
#     frete e reputacao junto ao comprador.
#   FALSO NEGATIVO: tratar como nao exportavel um lote que atingiria o
#     padrao - o lote vai ao mercado interno e o premio de exportacao e
#     perdido.
# Qual custa mais e uma decisao da cooperativa, nao do modelo. A razao
# de custos abaixo e uma SUPOSICAO declarada para este trabalho: a
# rejeicao de um lote ja comprometido com a exportacao (FP) e tomada
# como duas vezes mais cara do que o premio perdido (FN). A tabela mostra
# o que muda se a suposicao mudar.
#
# Regra de decisao com custos (dedutivel igualando o custo esperado das
# duas decisoes): classificar como exportavel quando
#     p > C_FP / (C_FP + C_FN).
# Com custos iguais, o corte e 0,5; se FP custa o dobro, 2/3; se FN custa
# o dobro, 1/3. O corte sobe quando o erro "para cima" custa mais.
CUSTO_FP <- 2
CUSTO_FN <- 1
corte_adotado <- CUSTO_FP / (CUSTO_FP + CUSTO_FN)

razoes_custo <- c(`1:1` = 1, `2:1` = 2, `3:1` = 3, `1:2` = 0.5)     # C_FP / C_FN
tab_cortes <- bind_rows(lapply(names(razoes_custo), function(nome) {
    r <- razoes_custo[[nome]]
    corte <- r / (1 + r)
    m <- metricas(prob_predita, resposta_exportacao, corte)
    data.frame(razao_FP_FN = nome, corte = corte,
               VP = m[["VP"]], FP = m[["FP"]], FN = m[["FN"]], VN = m[["VN"]],
               sensibilidade = m[["sensibilidade"]], especificidade = m[["especificidade"]],
               custo_por_lote = (m[["FP"]] * r + m[["FN"]] * 1) / n_obs)
}))
m_youden <- metricas(prob_predita, resposta_exportacao, corte_youden$threshold)
tab_cortes <- rbind(tab_cortes,
    data.frame(razao_FP_FN = "Youden (sem custos)", corte = corte_youden$threshold,
               VP = m_youden[["VP"]], FP = m_youden[["FP"]], FN = m_youden[["FN"]], VN = m_youden[["VN"]],
               sensibilidade = m_youden[["sensibilidade"]], especificidade = m_youden[["especificidade"]],
               custo_por_lote = NA))
# Regras triviais, para comparacao: "nunca apontar" (so falsos negativos:
# custo = n1 * C_FN / n) e "apontar todos" (so falsos positivos: custo =
# n0 * C_FP / n). Um modelo so tem valor de decisao se custar menos do
# que a melhor regra trivial para a mesma razao de custos.
tab_cortes$custo_nunca_apontar <- ifelse(is.na(tab_cortes$custo_por_lote), NA,
                                         frequencia_exportacao[["1"]] * 1 / n_obs)
tab_cortes$custo_apontar_todos <- ifelse(is.na(tab_cortes$custo_por_lote), NA,
                                         frequencia_exportacao[["0"]] * razoes_custo[tab_cortes$razao_FP_FN] / n_obs)
tab_cortes$ganho_vs_trivial <- pmin(tab_cortes$custo_nunca_apontar, tab_cortes$custo_apontar_todos, na.rm = FALSE) -
                               tab_cortes$custo_por_lote

cat("\ncortes por razao de custos (custos por lote em unidades de C_FN; ganho = melhor regra trivial - modelo):\n")
print(data.frame(razao_FP_FN    = tab_cortes$razao_FP_FN,
                 corte          = round(tab_cortes$corte, 3),
                 VP = tab_cortes$VP, FP = tab_cortes$FP, FN = tab_cortes$FN, VN = tab_cortes$VN,
                 sensibilidade  = round(tab_cortes$sensibilidade, 3),
                 especificidade = round(tab_cortes$especificidade, 3),
                 custo_modelo   = round(tab_cortes$custo_por_lote, 3),
                 custo_nunca    = round(tab_cortes$custo_nunca_apontar, 3),
                 custo_todos    = round(tab_cortes$custo_apontar_todos, 3),
                 ganho          = round(tab_cortes$ganho_vs_trivial, 3)),
      row.names = FALSE)

# Custo esperado ao longo de todos os cortes, por razao de custos
grade_cortes <- seq(0.02, 0.98, by = 0.01)
df_custo <- bind_rows(lapply(names(razoes_custo), function(nome) {
    r <- razoes_custo[[nome]]
    custo <- sapply(grade_cortes, function(cte) {
        m <- metricas(prob_predita, resposta_exportacao, cte)
        (m[["FP"]] * r + m[["FN"]]) / n_obs
    })
    data.frame(razao_FP_FN = nome, corte = grade_cortes, custo = custo)
}))
df_custo$razao_FP_FN <- factor(df_custo$razao_FP_FN, levels = names(razoes_custo))
df_min <- df_custo %>% group_by(razao_FP_FN) %>% slice_min(custo, n = 1, with_ties = FALSE) %>% ungroup()

# Figura C8: custo esperado por corte
# Pergunta que responde: onde colocar a regua?
# Como ler: uma curva por razao de custos (C_FP : C_FN); eixo y = custo
# medio por lote, em unidades do custo de um falso negativo; o ponto
# marca o minimo empirico de cada curva; a linha vertical e o corte
# teorico C_FP / (C_FP + C_FN) da razao adotada.
# Limitacao: os custos sao suposicoes declaradas; o minimo empirico e
# calculado nos mesmos dados que ajustaram o modelo.
g_custo <- ggplot(df_custo, aes(x = corte, y = custo, colour = razao_FP_FN)) +
    geom_vline(xintercept = corte_adotado, linetype = "dashed", colour = CINZA) +
    geom_line(linewidth = .9) +
    geom_point(data = df_min, size = 3) +
    scale_colour_manual(values = c(`1:1` = CINZA, `2:1` = VERDE, `3:1` = "#8A6D3B", `1:2` = VERMELHO),
                        name = "custo FP : custo FN") +
    labs(x = "ponto de corte (probabilidade)", y = "custo esperado por lote (em custos de FN)",
         title = "Custo esperado de classificacao ao longo do ponto de corte",
         subtitle = sprintf("linha tracejada: corte teorico %.3f para a razao adotada %d:%d; pontos: minimo empirico de cada curva",
                            corte_adotado, CUSTO_FP, CUSTO_FN),
         caption = "custos sao suposicoes declaradas; quem decide a razao e a cooperativa") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
salvar(g_custo, "C8_custo_por_corte.png", 9, 6)

# Matriz no corte adotado, e o minimo empirico da curva de custo para a
# razao adotada (pode diferir do corte teorico se o modelo nao estiver
# perfeitamente calibrado ou por variacao amostral).
m_adotado <- metricas(prob_predita, resposta_exportacao, corte_adotado)
cat(sprintf("\ncorte adotado (%d:%d): %.3f -> VP = %d | FP = %d | FN = %d | VN = %d\n",
            CUSTO_FP, CUSTO_FN, corte_adotado, m_adotado[["VP"]], m_adotado[["FP"]], m_adotado[["FN"]], m_adotado[["VN"]]))
cat(sprintf("   sensibilidade = %.3f | especificidade = %.3f | VPP = %.3f | VPN = %.3f\n",
            m_adotado[["sensibilidade"]], m_adotado[["especificidade"]], m_adotado[["VPP"]], m_adotado[["VPN"]]))
linha_adotada <- tab_cortes[tab_cortes$razao_FP_FN == sprintf("%d:%d", CUSTO_FP, CUSTO_FN), ]
min_adotado   <- df_min[df_min$razao_FP_FN == sprintf("%d:%d", CUSTO_FP, CUSTO_FN), ]
cat(sprintf("   custo por lote: modelo %.3f | nunca apontar %.3f | apontar todos %.3f | ganho do modelo %.3f\n",
            linha_adotada$custo_por_lote, linha_adotada$custo_nunca_apontar,
            linha_adotada$custo_apontar_todos, linha_adotada$ganho_vs_trivial))
cat(sprintf("   minimo empirico da curva de custo para %d:%d: corte %.2f, custo %.3f\n",
            CUSTO_FP, CUSTO_FN, min_adotado$corte, min_adotado$custo))

# ---------------------------------------------------------------------
# Parecer do Step 8
# ---------------------------------------------------------------------
# CONCLUSAO LITERAL: quantas vezes acertariamos; o que muda com a regua.
# CONCLUSAO ESTATISTICA: AUC, otimismo, regra de corte por custo.
# RESPOSTA A TAREFA 3: matriz, ROC, AUC e ponto de corte discutido.
# O QUE FICA: onde esta a variacao nao explicada.
parecer("CONCLUSAO LITERAL. Com a regua ingenua de 0,5, o modelo aponta ",
        sprintf("%.0f%%", 100 * m_05[["sensibilidade"]]), " dos lotes que atingiram o padrao ",
        "e deixa de fora ", sprintf("%.0f%%", 100 * m_05[["especificidade"]]),
        " dos que nao atingiram (acuracia ", sprintf("%.0f%%", 100 * m_05[["acuracia"]]),
        ", contra ", sprintf("%.0f%%", 100 * max(proporcao_exportacao, 1 - proporcao_exportacao)),
        " de quem sempre apostasse na classe maior). Com a razao de custos adotada (",
        CUSTO_FP, ":", CUSTO_FN, ", rejeicao na exportacao mais cara que premio perdido), a ",
        "regua sobe para ", sprintf("%.2f", corte_adotado), ": a sensibilidade cai para ",
        sprintf("%.0f%%", 100 * m_adotado[["sensibilidade"]]), " e a especificidade sobe para ",
        sprintf("%.0f%%", 100 * m_adotado[["especificidade"]]),
        " - o modelo fica mais exigente para apontar um lote como exportavel. ",
        "Comparado com as regras triviais, o modelo nesse corte custa ",
        sprintf("%.3f", linha_adotada$custo_por_lote), " por lote contra ",
        sprintf("%.3f", linha_adotada$custo_nunca_apontar), " de 'nunca apontar' e ",
        sprintf("%.3f", linha_adotada$custo_apontar_todos), " de 'apontar todos': ",
        if (linha_adotada$ganho_vs_trivial > 0.02)
            "ha ganho de decisao em usa-lo. "
        else if (linha_adotada$ganho_vs_trivial > 0)
            "o ganho sobre a melhor regra trivial e pequeno - sob essa razao de custos, o modelo acrescenta pouco a decisao. "
        else "ele nao supera a melhor regra trivial - sob essa razao de custos, o modelo nao tem valor de decisao. ",
        "CONCLUSAO ESTATISTICA. A AUC e ", sprintf("%.3f", auc_amostra), " [",
        sprintf("%.3f; %.3f", ic_auc[1], ic_auc[3]), "] na propria amostra e ",
        sprintf("%.3f", auc_cv), " por validacao cruzada em ", K_FOLDS,
        " partes: o modelo ordena melhor que o acaso, mas com capacidade ",
        if (auc_cv < 0.7) "modesta" else if (auc_cv < 0.8) "moderada" else "boa",
        ", e o otimismo da avaliacao na amostra e de ",
        sprintf("%.3f", auc_amostra - auc_cv), ". O ponto de corte nao e uma propriedade ",
        "do modelo: e a decisao p > C_FP/(C_FP + C_FN), que vale 0,5 so com custos ",
        "iguais; o corte de Youden (", sprintf("%.2f", corte_youden$threshold),
        ") e reportado como referencia sem custos. ",
        "RESPOSTA A TAREFA 3. Matriz de confusao, curva ROC e AUC estao calculadas; a ",
        "escolha do corte foi discutida em funcao do custo assimetrico, com a razao de ",
        "custos declarada como suposicao e a tabela mostrando o efeito de razoes ",
        "alternativas. ",
        "O QUE FICA PARA OS PROXIMOS STEPS: a capacidade preditiva modesta indica que ",
        "boa parte do que decide o padrao de exportacao nao esta nas caracteristicas do ",
        "talhao medidas nesta base (colheita, pos-colheita, beneficiamento). A ",
        "recomendacao do Step 9 precisa dizer isso.")


# =====================================================================
# ---- STEP 9 ---- CENARIOS, CONCLUSAO E RECOMENDACAO
# =====================================================================
secao("STEP 9 - CENARIOS E RECOMENDACAO (REQUISITO D)")

# Pergunta do step: para um produtor concreto, qual e a probabilidade de
# exportar sob diferentes escolhas - e o que a cooperativa deve
# recomendar para a proxima safra?
#
# Da razao de chances a probabilidade (slides, secao 12, item 6): os
# coeficientes vivem na escala do preditor linear; interpretar exige
# voltar a escala da resposta. Para um cenario, soma-se o preditor
# linear (intercepto + coeficientes das caracteristicas escolhidas) e
# converte-se com plogis(). Como a ligacao e nao linear, a mesma razao
# de chances produz diferencas em pontos percentuais que DEPENDEM do
# ponto de partida (Step 7).
#
# Intervalo de confianca da predicao: calculado na escala do preditor
# linear (eta +- 1,96 EP) e so depois convertido com plogis(). Um
# intervalo feito como p +- 1,96 EP poderia sair de [0, 1].
#
# "Mantendo as demais constantes": as numericas ficam na mediana da base
# e os fatores no nivel indicado em cada cenario. A escolha e declarada
# aqui, e nao escondida.
#
# OBJETIVO: traduzir o modelo em probabilidades para cenarios concretos.
# O QUE O CODIGO FAZ: monta a grade cultivar x manejo com altitude e
#   idade na mediana; predict(type = "link", se.fit = TRUE); plogis()
#   para a probabilidade e o IC; diferencas em pontos percentuais.
# COMO INTERPRETAR: a probabilidade de cada cenario e a leitura para o
#   produtor; os IC largos sao parte da resposta (n = 340).
mediana_altitude <- median(dados_cafe$altitude_m)
mediana_idade    <- median(dados_cafe$idade_lavoura_anos)

grade_cenarios <- expand.grid(cultivar = levels(dados_cafe$cultivar),
                              manejo   = levels(dados_cafe$manejo),
                              altitude_m = mediana_altitude,
                              idade_lavoura_anos = mediana_idade)
prever <- function(novos) {
    pr <- predict(modelo_final, newdata = novos, type = "link", se.fit = TRUE)
    novos$p      <- plogis(pr$fit)
    novos$ic_inf <- plogis(pr$fit - qnorm(0.975) * pr$se.fit)
    novos$ic_sup <- plogis(pr$fit + qnorm(0.975) * pr$se.fit)
    novos
}
pred_cenarios <- prever(grade_cenarios)
p_ref <- pred_cenarios$p[pred_cenarios$cultivar == levels(dados_cafe$cultivar)[1] &
                         pred_cenarios$manejo   == levels(dados_cafe$manejo)[1]]
pred_cenarios$dif_pp_vs_ref <- 100 * (pred_cenarios$p - p_ref)

cat(sprintf("talhao tipico: altitude %.0f m e lavoura de %.1f anos (medianas da base)\n",
            mediana_altitude, mediana_idade))
cat("probabilidade de atingir o padrao, por cultivar x manejo:\n")
print(data.frame(cultivar = pred_cenarios$cultivar, manejo = pred_cenarios$manejo,
                 p = round(pred_cenarios$p, 3),
                 IC_inf = round(pred_cenarios$ic_inf, 3), IC_sup = round(pred_cenarios$ic_sup, 3),
                 dif_pp_vs_ref = round(pred_cenarios$dif_pp_vs_ref, 1)),
      row.names = FALSE)
cat(sprintf("   referencia (dif_pp_vs_ref = 0): %s, %s\n",
            levels(dados_cafe$cultivar)[1], levels(dados_cafe$manejo)[1]))

# O mesmo manejo, pontos de partida diferentes: ganho em pp do manejo
# com maior probabilidade em relacao ao convencional, dentro de cada
# cultivar. Mostra que a mesma razao de chances vale pontos percentuais
# diferentes conforme a cultivar.
manejo_ref <- levels(dados_cafe$manejo)[1]
ganho_manejo <- pred_cenarios %>%
    group_by(cultivar) %>%
    summarise(p_convencional = p[manejo == manejo_ref],
              melhor_manejo  = as.character(manejo[which.max(p)]),
              p_melhor       = max(p),
              ganho_pp       = 100 * (p_melhor - p_convencional), .groups = "drop")
cat("\nganho em pontos percentuais do melhor manejo sobre o convencional, por cultivar:\n")
print(data.frame(cultivar = ganho_manejo$cultivar, p_convencional = round(ganho_manejo$p_convencional, 3),
                 melhor_manejo = ganho_manejo$melhor_manejo, p_melhor = round(ganho_manejo$p_melhor, 3),
                 ganho_pp = round(ganho_manejo$ganho_pp, 1)),
      row.names = FALSE)

# Altitude e idade: variacoes dentro do suporte dos dados (quartis), no
# cenario de referencia (cultivar e manejo de referencia).
q_alt   <- quantile(dados_cafe$altitude_m, c(.25, .75))
q_idade <- quantile(dados_cafe$idade_lavoura_anos, c(.25, .75))
cen_num <- rbind(
    data.frame(cenario = sprintf("altitude %.0f m (Q1)", q_alt[1]), altitude_m = q_alt[1], idade_lavoura_anos = mediana_idade),
    data.frame(cenario = sprintf("altitude %.0f m (Q3)", q_alt[2]), altitude_m = q_alt[2], idade_lavoura_anos = mediana_idade),
    data.frame(cenario = sprintf("lavoura de %.1f anos (Q1)", q_idade[1]), altitude_m = mediana_altitude, idade_lavoura_anos = q_idade[1]),
    data.frame(cenario = sprintf("lavoura de %.1f anos (Q3)", q_idade[2]), altitude_m = mediana_altitude, idade_lavoura_anos = q_idade[2]))
cen_num$cultivar <- factor(levels(dados_cafe$cultivar)[1], levels = levels(dados_cafe$cultivar))
cen_num$manejo   <- factor(manejo_ref, levels = levels(dados_cafe$manejo))
pred_num <- prever(cen_num)
cat(sprintf("\nvariacoes de altitude e idade no cenario de referencia (%s, %s):\n",
            levels(dados_cafe$cultivar)[1], manejo_ref))
print(data.frame(cenario = pred_num$cenario, p = round(pred_num$p, 3),
                 IC_inf = round(pred_num$ic_inf, 3), IC_sup = round(pred_num$ic_sup, 3)),
      row.names = FALSE)
dif_alt_pp   <- 100 * (pred_num$p[2] - pred_num$p[1])
dif_idade_pp <- 100 * (pred_num$p[4] - pred_num$p[3])
cat(sprintf("   de Q1 a Q3 de altitude: %+.1f pontos percentuais | de Q1 a Q3 de idade: %+.1f pontos percentuais\n",
            dif_alt_pp, dif_idade_pp))

# ---------------------------------------------------------------------
# Figura C9: cenarios cultivar x manejo
# ---------------------------------------------------------------------
# Pergunta que responde: qual a probabilidade de exportar em cada escolha?
# Como ler: um painel por manejo; ponto = probabilidade prevista para o
# talhao tipico; barra = IC 95% (calculado no logit e convertido); linha
# tracejada = proporcao geral da base. Comparar alturas e a sobreposicao
# dos intervalos.
# Limitacao: altitude e idade fixas na mediana; outros valores mudam as
# probabilidades (mas nao as razoes de chances).
g_cenarios <- ggplot(pred_cenarios, aes(x = cultivar, y = p)) +
    geom_hline(yintercept = proporcao_exportacao, linetype = "dashed", colour = CINZA) +
    geom_errorbar(aes(ymin = ic_inf, ymax = ic_sup), width = .2) +
    geom_point(size = 3.5, colour = VERDE) +
    geom_text(aes(label = sprintf("%.2f", p)), hjust = -.45, size = 3.2) +
    facet_wrap(~ manejo) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
    labs(x = NULL, y = "probabilidade de atingir o padrao de exportacao",
         title = "Probabilidade de exportar por cultivar e manejo, num talhao tipico",
         subtitle = sprintf("altitude %.0f m e lavoura de %.1f anos (medianas); barras: IC 95%%; tracejada: proporcao geral",
                            mediana_altitude, mediana_idade),
         caption = "IC calculado na escala do logit e convertido; os pontos percentuais de cada troca dependem do ponto de partida") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 15, hjust = 1))
salvar(g_cenarios, "C9_cenarios.png", 11, 6)

# ---------------------------------------------------------------------
# Figura C9: curvas ao longo da altitude e da idade
# ---------------------------------------------------------------------
# Pergunta que responde: como altitude e idade movem a probabilidade?
# Como ler: uma curva por cultivar, manejo de referencia; faixa = IC 95%.
# A inclinacao muda ao longo do eixo: e a nao linearidade da ligacao na
# escala da probabilidade. Curvas so dentro da faixa observada dos dados.
grade_alt <- expand.grid(cultivar = levels(dados_cafe$cultivar), manejo = manejo_ref,
                         altitude_m = seq(min(dados_cafe$altitude_m), max(dados_cafe$altitude_m), length.out = 80),
                         idade_lavoura_anos = mediana_idade)
grade_idade <- expand.grid(cultivar = levels(dados_cafe$cultivar), manejo = manejo_ref,
                           altitude_m = mediana_altitude,
                           idade_lavoura_anos = seq(min(dados_cafe$idade_lavoura_anos), max(dados_cafe$idade_lavoura_anos), length.out = 80))
grade_alt$manejo   <- factor(grade_alt$manejo, levels = levels(dados_cafe$manejo))
grade_idade$manejo <- factor(grade_idade$manejo, levels = levels(dados_cafe$manejo))
pred_alt   <- prever(grade_alt)
pred_idade <- prever(grade_idade)
cores_cultivar <- c(VERDE, "#4A7FA5", VERMELHO, "#8A6D3B")
names(cores_cultivar) <- levels(dados_cafe$cultivar)
g_alt <- ggplot(pred_alt, aes(x = altitude_m, y = p, colour = cultivar, fill = cultivar)) +
    geom_ribbon(aes(ymin = ic_inf, ymax = ic_sup), alpha = .12, colour = NA) +
    geom_line(linewidth = 1) +
    scale_colour_manual(values = cores_cultivar) + scale_fill_manual(values = cores_cultivar) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "altitude (m)", y = "probabilidade de exportar",
         title = "Ao longo da altitude", subtitle = sprintf("manejo %s; lavoura de %.1f anos", manejo_ref, mediana_idade)) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
g_idade <- ggplot(pred_idade, aes(x = idade_lavoura_anos, y = p, colour = cultivar, fill = cultivar)) +
    geom_ribbon(aes(ymin = ic_inf, ymax = ic_sup), alpha = .12, colour = NA) +
    geom_line(linewidth = 1) +
    scale_colour_manual(values = cores_cultivar) + scale_fill_manual(values = cores_cultivar) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "idade da lavoura (anos)", y = NULL,
         title = "Ao longo da idade da lavoura", subtitle = sprintf("manejo %s; altitude %.0f m", manejo_ref, mediana_altitude)) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
salvar((g_alt | g_idade) + plot_layout(guides = "collect") & theme(legend.position = "bottom"),
       "C9_altitude_idade.png", 12, 6)

# ---------------------------------------------------------------------
# Reconciliacao com o Desafio B: resistencia a ferrugem x exportacao
# ---------------------------------------------------------------------
# A cultivar mais resistente a ferrugem (menor severidade media, Step 6)
# e a que mais exporta? A tabela coloca as duas ordenacoes lado a lado,
# no cenario tipico com manejo de referencia.
tab_reconcilia <- pred_cenarios %>%
    filter(manejo == manejo_ref) %>%
    transmute(cultivar = as.character(cultivar), p_exportar = p) %>%
    mutate(severidade_media = as.vector(sev_media[cultivar])) %>%
    arrange(desc(p_exportar))
cat("\ncultivares: probabilidade de exportar (talhao tipico, manejo de referencia) x severidade media de ferrugem:\n")
print(data.frame(cultivar = tab_reconcilia$cultivar, p_exportar = round(tab_reconcilia$p_exportar, 3),
                 severidade_media = round(tab_reconcilia$severidade_media, 3)), row.names = FALSE)
cultivar_mais_exporta <- tab_reconcilia$cultivar[1]
cultivar_mais_resist  <- tab_reconcilia$cultivar[which.min(tab_reconcilia$severidade_media)]

# ---------------------------------------------------------------------
# Parecer final do Desafio C
# ---------------------------------------------------------------------
# RECOMENDACAO: o que fazer, para quem, com que efeito em probabilidade.
# MANEJAVEL x NAO MANEJAVEL: o que se escolhe e o que so orienta onde agir.
# RECONCILIACAO COM A e B: a cultivar resistente e a que exporta?
# LIMITES: o que o modelo nao consegue dizer.
melhor_cenario <- pred_cenarios[which.max(pred_cenarios$p), ]
pior_cenario   <- pred_cenarios[which.min(pred_cenarios$p), ]

parecer("RECOMENDACAO. Num talhao tipico (altitude ", sprintf("%.0f", mediana_altitude),
        " m, lavoura de ", sprintf("%.1f", mediana_idade), " anos), a probabilidade de ",
        "atingir o padrao de exportacao vai de ", sprintf("%.2f", pior_cenario$p), " (",
        pior_cenario$cultivar, ", ", pior_cenario$manejo, ") a ",
        sprintf("%.2f", melhor_cenario$p), " (", melhor_cenario$cultivar, ", ",
        melhor_cenario$manejo, "), com intervalos largos (n = ", n_obs, "). Entre as ",
        "escolhas do produtor: (1) o manejo - trocar o convencional pelo ",
        ganho_manejo$melhor_manejo[1], " acrescenta de ",
        sprintf("%.0f", min(ganho_manejo$ganho_pp)), " a ",
        sprintf("%.0f", max(ganho_manejo$ganho_pp)), " pontos percentuais",
        if (diff(range(ganho_manejo$ganho_pp)) < 3)
            ", quase o mesmo em todas as cultivares porque todas partem de probabilidades proximas do meio da escala, onde a mesma razao de chances rende mais pontos percentuais"
        else ", conforme a cultivar, porque a mesma razao de chances rende mais pontos percentuais perto de 50% do que longe dele",
        "; (2) a cultivar - ", cultivar_mais_exporta,
        " tem a maior probabilidade de exportar no cenario tipico, e o efeito total ",
        "(Step 6) e o que vale para quem escolhe antes do plantio; (3) a idade da ",
        "lavoura - entre o primeiro e o terceiro quartil de idade a probabilidade muda ",
        sprintf("%+.0f", dif_idade_pp), " pontos percentuais, o que sustenta programar a ",
        "renovacao das lavouras mais velhas. ",
        "MANEJAVEL E NAO MANEJAVEL. A altitude nao se escolhe, mas orienta onde ",
        "concentrar o esforco: entre o primeiro e o terceiro quartil de altitude a ",
        "probabilidade muda ", sprintf("%+.0f", dif_alt_pp), " pontos percentuais; talhoes ",
        "altos sao os candidatos naturais aos contratos de exportacao, e talhoes baixos ",
        "os que mais precisam do manejo para compensar. ",
        "RECONCILIACAO COM OS DESAFIOS A E B. A cultivar mais resistente a ferrugem e ",
        cultivar_mais_resist, " (severidade media ",
        sprintf("%.3f", tab_reconcilia$severidade_media[tab_reconcilia$cultivar == cultivar_mais_resist]),
        "), e a que mais exporta e ", cultivar_mais_exporta, " (severidade media ",
        sprintf("%.3f", tab_reconcilia$severidade_media[tab_reconcilia$cultivar == cultivar_mais_exporta]),
        "). ",
        if (cultivar_mais_resist != cultivar_mais_exporta)
            paste0("Nao sao a mesma: a resistencia a ferrugem, sozinha, nao garante o padrao ",
                   "de exportacao, e a cooperativa precisa decidir o objetivo - reduzir a doenca ",
                   "ou maximizar a exportacao - ou combinar cultivares por talhao conforme a ",
                   "pressao de ferrugem local (Desafio B) e a altitude. ")
        else "Sao a mesma cultivar, o que simplifica a recomendacao. ",
        "A contagem de brocas (Desafio A) nao acrescentou informacao ao modelo (Step 4). ",
        "LIMITES. O modelo ordena os talhoes melhor que o acaso, mas com capacidade ",
        "modesta (AUC ", sprintf("%.2f", auc_cv), " por validacao cruzada): boa parte do ",
        "que decide o padrao de exportacao nao esta nas caracteristicas de talhao desta ",
        "base - colheita, pos-colheita e beneficiamento nao foram medidos. As ",
        "probabilidades acima valem dentro da faixa observada de altitude e idade e para ",
        "as combinacoes de cultivar e manejo presentes na base; os intervalos de ",
        "confianca devem acompanhar qualquer numero repassado ao produtor. E as ",
        "associacoes estimadas sao observacionais: o modelo nao substitui um ensaio em ",
        "que cultivar e manejo sejam atribuidos aos talhoes.")

# ---------------------------------------------------------------------
# Reprodutibilidade: versoes usadas nesta execucao
# ---------------------------------------------------------------------
# Para o relatorio, os pareceres de todos os steps podem ser extraidos
# da saida deste script com:
#   grep -A40 ">> PARECER" saida.txt
cat("\nfiguras salvas em:", FIG, "\n")
cat("\nversoes desta execucao:\n")
print(sessionInfo())


# =======================================================================================
#                  Creative Commons License 4.0
#                       (CC BY-NC-SA 4.0)
#
#  This is a human-readable summary of (and not a substitute for) the
#  license (https://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
#
#  TODO (licenca): o nome acima (BY-NC-SA) e o endereco (by-nc-nd) apontam
#  para licencas diferentes, e o arquivo LICENSE na raiz do repositorio e
#  CC0 1.0. A mesma inconsistencia existe nos scripts A, B e da AED.
#  Rodape preservado ate o grupo decidir qual licenca adotar.
# =======================================================================================
