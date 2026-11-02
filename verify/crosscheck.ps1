# Cross-check the certified result (pub\lambda_30000.txt) against other files:
#  * the copy of the independent Python computation shipped in verify\ (always
#    present in this package);
#  * the older files of the same project kept next to pub\ in
#    improved\golomb_dickman (optional: reported as MISSING if absent).
#
# "first" is the index inside the string "0."+digits, so 30002 means "identical
# on all 30000 digits". Comparison uses the first line starting with "0.".
#
# Run:  powershell -ExecutionPolicy Bypass -File verify\crosscheck.ps1
# (ASCII only on purpose: Windows PowerShell 5.1 reads .ps1 without BOM as ANSI.)

param([string]$Pub = $PSScriptRoot)

$ErrorActionPreference = 'Stop'
$prog = Split-Path $Pub -Parent                     # program directory
$new  = @(Get-Content -LiteralPath (Join-Path $Pub 'lambda_30000.txt') |
          Where-Object { $_ -match '^0\.' })[0]

# label, path, description
$cands = @(
    @('verify/lambda_30000_fix_python.txt', (Join-Path $Pub 'verify\lambda_30000_fix_python.txt'),
      'python golomb_dickman.py, 23.10 (run30kfix)'),
    @('30000/out/lambda_30000_fixed.txt', (Join-Path $prog '30000\out\lambda_30000_fixed.txt'),
      'main_fixed.py, 15.10'),
    @('30000/out/lambda_30000.txt', (Join-Path $prog '30000\out\lambda_30000.txt'),
      'main.py, 13.10 (withdrawn)'),
    @('lambda_30000_cpp.txt', (Join-Path $prog 'lambda_30000_cpp.txt'),
      'earlier C++, 25.10'),
    @('lambda_30000_clai.txt', (Join-Path $prog 'lambda_30000_clai.txt'),
      'earlier python')
)

'{0,-34} {1,-44} {2,6} {3,6}  {4}' -f 'file', 'what', 'len', 'first', 'ball.radius (header)'
'-' * 122
foreach ($c in $cands) {
    $p = $c[1]
    if (-not (Test-Path -LiteralPath $p)) { '{0,-34} MISSING' -f $c[0]; continue }
    $hdr  = Get-Content -LiteralPath $p -TotalCount 1
    $ball = 'none'
    if ($hdr -match 'ball\.radius=10\^(-?\d+)') { $ball = '10^' + $Matches[1] }
    if ($hdr -match 'NOT CONFIRMED|\*\*\*')     { $ball += '  ***' }
    $o = @(Get-Content -LiteralPath $p | Where-Object { $_ -match '^0\.' })[0]
    $n = [Math]::Min($new.Length, $o.Length)
    $i = 0
    while ($i -lt $n -and $new[$i] -eq $o[$i]) { $i++ }
    '{0,-34} {1,-44} {2,6} {3,6}  {4}' -f $c[0], $c[2], $o.Length, $i, $ball
}
