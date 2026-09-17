param(
    [string]$AtlasTable,
    [int]$Port = 8765
)
$ErrorActionPreference = 'Stop'
$repoRoot = Split-Path -Parent $PSScriptRoot
$pythonPath = Join-Path $repoRoot '.venv\Scripts\python.exe'
if (-not (Test-Path -LiteralPath $pythonPath)) {
    throw 'Create the isolated environment first: python -m venv .venv; .\.venv\Scripts\python.exe -m pip install -e . -r requirements-model.txt'
}
if ($AtlasTable) {
    $AtlasTable = (Resolve-Path -LiteralPath $AtlasTable).Path
}
Push-Location -LiteralPath $repoRoot
try {
    if ($AtlasTable) {
        & $pythonPath 'bridge\serve_table.py' $AtlasTable --port $Port
    } else {
        & $pythonPath -m quadcond.service --port $Port
    }
    if ($LASTEXITCODE -ne 0) { throw "Workbench exited with code $LASTEXITCODE" }
} finally {
    Pop-Location
}
