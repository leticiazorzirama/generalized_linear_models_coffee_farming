#!/usr/bin/env bash
# Atalho Linux / macOS. Toda a logica esta em build.R, que roda igual nos tres sistemas.
#   ./relatorio/build.sh
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if ! command -v Rscript >/dev/null 2>&1; then
    echo "Rscript nao encontrado no PATH."
    echo "  Ubuntu/Debian : sudo apt install r-base"
    echo "  macOS         : brew install r     (ou https://cran.r-project.org)"
    exit 1
fi

exec Rscript "$DIR/build.R" "$@"
