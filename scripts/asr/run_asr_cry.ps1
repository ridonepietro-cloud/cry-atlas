# Ancestors (PAML codeml + GRASP) of the full-toxin trees of the Cry pack -> site/data/cry/asr_<key>.json + asr_index.json.
# The tree columns are those occupied by >= 20 % of the entries (20_trees.py), so ASR_MIN_OCC = 0.2.   Usage: powershell -File run_asr_cry.ps1 struct_A [seq_A]
param([string[]]$keys = @("struct_A"))
$AS = Split-Path -Parent $MyInvocation.MyCommand.Path; Set-Location $AS
New-Item -ItemType Directory -Force -Path (Join-Path $AS "out") | Out-Null
$env:ASR_MIN_OCC = '0.2'; $env:ASR_PAML_TIMEOUT = '43200'; $env:ASR_GRASP_THREADS = '4'; $env:ASR_GRASP_MEM = '6g'; $env:PYTHONIOENCODING = 'utf-8'
foreach ($k in $keys) {
    $log = Join-Path $AS "out\asr_cry_$k.log"
    "$(Get-Date -Format 'HH:mm:ss') start $k" | Add-Content $log
    & py asr_export.py cry $k *>> $log
    "$(Get-Date -Format 'HH:mm:ss') end $k (exit $LASTEXITCODE)" | Add-Content $log
}
