# Atalho Windows. Toda a logica esta em build.R, que roda igual nos tres sistemas.
#   .\relatorio\build.ps1
$ErrorActionPreference = 'Stop'
$dir = $PSScriptRoot

$rscript = Get-Command Rscript -ErrorAction SilentlyContinue
if (-not $rscript) {
    # Instalacao padrao do R no Windows nao entra no PATH sozinha
    $cand = Get-ChildItem 'C:\Program Files\R\*\bin\Rscript.exe' -ErrorAction SilentlyContinue |
            Sort-Object FullName -Descending | Select-Object -First 1
    if ($cand) {
        $rscript = $cand.FullName
    } else {
        Write-Host 'Rscript nao encontrado.' -ForegroundColor Red
        Write-Host 'Instale o R em https://cran.r-project.org/bin/windows/base/'
        exit 1
    }
} else {
    $rscript = $rscript.Source
}

& $rscript (Join-Path $dir 'build.R') @args
