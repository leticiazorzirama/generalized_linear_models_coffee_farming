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
