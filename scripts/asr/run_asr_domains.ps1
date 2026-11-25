# Ancestors (PAML codeml + GRASP PSP in interleaved blocks) of the six domain trees of the Cry pack, all at once (codeml is single-threaded).
# ASR_MIN_OCC = 0.1 = the occupancy of the columns of the domain trees (20_trees.py); output site/data/cry/asr_<key>.json + asr_index.json (merged).
$AS = Split-Path -Parent $MyInvocation.MyCommand.Path; Set-Location $AS
New-Item -ItemType Directory -Force -Path (Join-Path $AS "out") | Out-Null
$env:ASR_MIN_OCC = '0.1'; $env:ASR_PAML_TIMEOUT = '43200'; $env:ASR_GRASP_THREADS = '2'; $env:ASR_GRASP_MEM = '2g'; $env:ASR_GRASP_CHUNKS = '2'; $env:PYTHONIOENCODING = 'utf-8'
$procs = @()
foreach ($k in "struct_DI", "struct_DII", "struct_DIII", "seq_DI", "seq_DII", "seq_DIII") {
    $procs += Start-Process -FilePath "py" -ArgumentList @("asr_export.py", "cry", $k, "--indel", "PSP") -WorkingDirectory $AS -NoNewWindow -PassThru -RedirectStandardOutput (Join-Path $AS "out\asr_cry_$k.out.txt") -RedirectStandardError (Join-Path $AS "out\asr_cry_$k.err.txt")
}
$procs | Wait-Process
"all six finished $(Get-Date -Format 'HH:mm:ss')" | Out-File -Append -Encoding utf8 (Join-Path $AS "out\asr_cry_domains.log")
