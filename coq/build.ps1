# Сборка Coq-файлов без make (Windows).
# Требуется coqc из состава Coq и установленный mathcomp.

$ErrorActionPreference = "Stop"

$root = Split-Path -Parent $MyInvocation.MyCommand.Path
Push-Location $root
try {
    $files = @(
        "Definitions.v",
        "Ddet_succ.v",
        "Ddet_eq.v",
        "Ps_step3.v",
        "Qs_step3.v",
        "S_to_J_contraction.v",
        "All.v"
    )

    foreach ($f in $files) {
        Write-Host "==> coqc -Q . GolombDickmanCoq $f"
        & coqc -Q . GolombDickmanCoq $f
        if ($LASTEXITCODE -ne 0) {
            throw "coqc failed for $f"
        }
    }

    Write-Host "OK: all Coq files compiled."
}
finally {
    Pop-Location
}
