#!/usr/bin/env bash
# =======================================================================================
# BOOTSTRAP DO AMBIENTE - Linux e macOS
#
#   ./relatorio/setup.sh
#
# Resolve o caso em que a maquina nao tem NADA: acha o R, e se nao existir,
# descobre o gerenciador de pacotes e diz o comando exato. Depois entrega o
# trabalho pesado para o setup.R, que roda igual nos tres sistemas.
#
# Nao instala nada sem perguntar.
# =======================================================================================
set -uo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VERDE='\033[0;32m'; VERM='\033[0;31m'; AMAR='\033[0;33m'; NEUTRO='\033[0m'
ok()    { printf "${VERDE}  [ok]   ${NEUTRO}%s\n" "$1"; }
aviso() { printf "${AMAR}  [!]    ${NEUTRO}%s\n" "$1"; }
erro()  { printf "${VERM}  [FALHA]${NEUTRO}%s\n" "$1"; }

echo
echo "======================================================================"
echo "BOOTSTRAP DO AMBIENTE - RELATORIO PBL CAFEICULTURA"
echo "======================================================================"

# --- sistema ---------------------------------------------------------------------------
case "$(uname -s)" in
    Linux*)  SO="Linux" ;;
    Darwin*) SO="macOS" ;;
    *)       SO="$(uname -s)" ;;
esac
echo "  sistema : $SO ($(uname -m))"

# --- gerenciador de pacotes -------------------------------------------------------------
GERENCIADOR=""; INSTALA_R=""
if   command -v brew    >/dev/null 2>&1; then GERENCIADOR="brew";    INSTALA_R="brew install r"
elif command -v apt-get >/dev/null 2>&1; then GERENCIADOR="apt";     INSTALA_R="sudo apt-get update && sudo apt-get install -y r-base"
elif command -v dnf     >/dev/null 2>&1; then GERENCIADOR="dnf";     INSTALA_R="sudo dnf install -y R"
elif command -v pacman  >/dev/null 2>&1; then GERENCIADOR="pacman";  INSTALA_R="sudo pacman -S --noconfirm r"
elif command -v zypper  >/dev/null 2>&1; then GERENCIADOR="zypper";  INSTALA_R="sudo zypper install -y R-base"
elif command -v apk     >/dev/null 2>&1; then GERENCIADOR="apk";     INSTALA_R="sudo apk add R R-dev"
fi
[ -n "$GERENCIADOR" ] && echo "  gerenciador: $GERENCIADOR" || echo "  gerenciador: nao identificado"
echo

# --- procurar o Rscript ------------------------------------------------------------------
echo "[1/3] Procurando o R"
echo "----------------------------------------------------------------------"
RSCRIPT=""
if command -v Rscript >/dev/null 2>&1; then
    RSCRIPT="$(command -v Rscript)"
else
    # instalacoes fora do PATH
    for c in /usr/local/bin/Rscript /opt/homebrew/bin/Rscript /usr/bin/Rscript \
             /Library/Frameworks/R.framework/Resources/bin/Rscript \
             /opt/R/*/bin/Rscript; do
        [ -x "$c" ] && { RSCRIPT="$c"; break; }
    done
fi

if [ -z "$RSCRIPT" ]; then
    erro "R nao encontrado."
    echo
    if [ -n "$INSTALA_R" ]; then
        echo "  Instale com:"
        echo "      $INSTALA_R"
        echo
        read -r -p "  Quer que eu rode esse comando agora? [s/N] " resp
        case "$resp" in
            [sS]|[sS][iI][mM]|[yY])
                echo
                eval "$INSTALA_R" || { erro "a instalacao falhou"; exit 1; }
                RSCRIPT="$(command -v Rscript || true)"
                ;;
            *) echo; echo "  Instale o R e rode este script de novo."; exit 1 ;;
        esac
    else
        echo "  Nao identifiquei seu gerenciador de pacotes."
        echo "  Baixe o R em https://cran.r-project.org e rode este script de novo."
        exit 1
    fi
fi

[ -z "$RSCRIPT" ] && { erro "R instalado mas Rscript nao ficou no PATH. Reabra o terminal."; exit 1; }
ok "Rscript: $RSCRIPT"
ok "$("$RSCRIPT" -e 'cat(R.version.string)' 2>/dev/null)"

# --- dependencias de sistema para compilar pacotes R --------------------------------------
echo
echo "[2/3] Dependencias de compilacao"
echo "----------------------------------------------------------------------"
# No Linux, varios pacotes R vem como fonte e precisam de cabecalhos de sistema.
# No macOS, os binarios do CRAN dispensam isso.
if [ "$SO" = "Linux" ] && [ -n "$GERENCIADOR" ]; then
    FALTANDO=""
    for h in /usr/include/curl/curl.h /usr/include/openssl/ssl.h; do
        [ -f "$h" ] || FALTANDO="sim"
    done
    if [ -n "$FALTANDO" ]; then
        aviso "faltam cabecalhos de desenvolvimento (curl/openssl)."
        echo "         Varios pacotes R sao compilados da fonte no Linux e precisam deles."
        case "$GERENCIADOR" in
            apt)    echo "         sudo apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev" ;;
            dnf)    echo "         sudo dnf install -y libcurl-devel openssl-devel libxml2-devel" ;;
            pacman) echo "         sudo pacman -S curl openssl libxml2" ;;
            zypper) echo "         sudo zypper install -y libcurl-devel libopenssl-devel libxml2-devel" ;;
        esac
        echo "         Se o passo 3 falhar instalando pacote, rode o comando acima."
    else
        ok "cabecalhos presentes"
    fi
else
    ok "nada a verificar neste sistema"
fi

# --- entregar ao setup.R -------------------------------------------------------------------
echo
echo "[3/3] Preparando R e LaTeX"
echo "----------------------------------------------------------------------"
echo "  (daqui em diante o setup.R assume - ele roda igual nos tres sistemas)"
exec "$RSCRIPT" "$DIR/setup.R"
