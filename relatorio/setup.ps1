# =======================================================================================
# BOOTSTRAP DO AMBIENTE - Windows
#
#   .\relatorio\setup.ps1
#
# Se o PowerShell recusar por politica de execucao, rode antes (uma vez so):
#   Set-ExecutionPolicy -Scope CurrentUser -ExecutionPolicy RemoteSigned
#
# Resolve o caso em que a maquina nao tem NADA: acha o R mesmo fora do PATH
# (instalador do Windows nao o adiciona), e oferece instalar via winget.
# Depois entrega o trabalho pesado ao setup.R, que roda igual nos tres sistemas.
#
# Nao instala nada sem perguntar.
# =======================================================================================
$ErrorActionPreference = 'Stop'
$dir = $PSScriptRoot

function Ok    ($m) { Write-Host "  [ok]   $m"   -ForegroundColor Green }
function Aviso ($m) { Write-Host "  [!]    $m"   -ForegroundColor Yellow }
function Falha ($m) { Write-Host "  [FALHA]$m"   -ForegroundColor Red }

Write-Host ""
Write-Host "======================================================================"
Write-Host "BOOTSTRAP DO AMBIENTE - RELATORIO PBL CAFEICULTURA"
Write-Host "======================================================================"
Write-Host "  sistema : Windows $([Environment]::OSVersion.Version) ($env:PROCESSOR_ARCHITECTURE)"
Write-Host "  shell   : PowerShell $($PSVersionTable.PSVersion)"
Write-Host ""

# =======================================================================================
Write-Host "[1/3] Procurando o R"
Write-Host "----------------------------------------------------------------------"

$rscript = $null

# a) no PATH
$cmd = Get-Command Rscript -ErrorAction SilentlyContinue
if ($cmd) { $rscript = $cmd.Source }

# b) instalacao padrao - o instalador do R nao mexe no PATH
if (-not $rscript) {
    $cand = Get-ChildItem 'C:\Program Files\R\*\bin\Rscript.exe',
                          'C:\Program Files (x86)\R\*\bin\Rscript.exe' `
                          -ErrorAction SilentlyContinue |
            Sort-Object FullName -Descending | Select-Object -First 1
    if ($cand) { $rscript = $cand.FullName }
}

# c) registro
if (-not $rscript) {
    foreach ($k in 'HKLM:\SOFTWARE\R-core\R', 'HKCU:\SOFTWARE\R-core\R') {
        try {
            $p = (Get-ItemProperty $k -ErrorAction Stop).InstallPath
            if ($p -and (Test-Path (Join-Path $p 'bin\Rscript.exe'))) {
                $rscript = Join-Path $p 'bin\Rscript.exe'; break
            }
        } catch {}
    }
}

# d) nao achou: oferecer instalacao
if (-not $rscript) {
    Falha "R nao encontrado."
    Write-Host ""
    $winget = Get-Command winget -ErrorAction SilentlyContinue
    if ($winget) {
        Write-Host "  Posso instalar via winget:"
        Write-Host "      winget install --id RProject.R -e"
        Write-Host ""
        $r = Read-Host "  Instalar agora? [s/N]"
        if ($r -match '^[sSyY]') {
            Write-Host ""
            & winget install --id RProject.R -e --accept-source-agreements --accept-package-agreements
            $cand = Get-ChildItem 'C:\Program Files\R\*\bin\Rscript.exe' -ErrorAction SilentlyContinue |
                    Sort-Object FullName -Descending | Select-Object -First 1
            if ($cand) { $rscript = $cand.FullName }
        }
    } else {
        Write-Host "  winget nao disponivel."
    }

    if (-not $rscript) {
        Write-Host ""
        Write-Host "  Baixe o R em https://cran.r-project.org/bin/windows/base/"
        Write-Host "  Instale e rode este script de novo."
        exit 1
    }
}

Ok "Rscript: $rscript"
Ok (& $rscript -e 'cat(R.version.string)' 2>$null)

# =======================================================================================
Write-Host ""
Write-Host "[2/3] Ambiente do Windows"
Write-Host "----------------------------------------------------------------------"

# Rtools nao e necessario: no Windows o CRAN entrega binarios pre-compilados.
if (Get-ChildItem 'C:\rtools*' -ErrorAction SilentlyContinue) {
    Ok "Rtools presente (nao e obrigatorio - o CRAN entrega binarios)"
} else {
    Ok "sem Rtools - tudo bem, o CRAN entrega binarios para Windows"
}

# Caminho com espaco quebra algumas ferramentas do TeX Live.
if ($dir -match ' ') {
    Aviso "o caminho do projeto tem espaco:"
    Write-Host "         $dir"
    Write-Host "         Costuma funcionar, mas se o LaTeX reclamar de arquivo nao"
    Write-Host "         encontrado, mova o repositorio para um caminho sem espacos."
} else {
    Ok "caminho sem espacos"
}

# OneDrive sincronizando .aux/.pdf durante a compilacao causa erro de arquivo travado.
if ($dir -match 'OneDrive|Dropbox|Google Drive') {
    Aviso "o projeto esta numa pasta sincronizada na nuvem."
    Write-Host "         A sincronizacao pode travar arquivos no meio da compilacao."
    Write-Host "         Se der 'permission denied', pause a sincronizacao e tente de novo."
}

# =======================================================================================
Write-Host ""
Write-Host "[3/3] Preparando R e LaTeX"
Write-Host "----------------------------------------------------------------------"
Write-Host "  (daqui em diante o setup.R assume - ele roda igual nos tres sistemas)"

& $rscript (Join-Path $dir 'setup.R')
exit $LASTEXITCODE
