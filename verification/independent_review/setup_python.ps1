param([string]$BasePython='python')
$ErrorActionPreference='Stop'
$verificationDir=Split-Path -Parent $PSScriptRoot
$venvDir=Join-Path $verificationDir '.venv'
$pythonExe=Join-Path $venvDir 'Scripts\python.exe'
if(-not (Test-Path -LiteralPath $pythonExe)) {
 & $BasePython -m venv $venvDir
 if($LASTEXITCODE -ne 0) { throw 'Could not create the dedicated Python environment.' }
}
$lockFile=Join-Path $verificationDir 'requirements.lock.txt'
& $pythonExe -m pip install -r $lockFile
if($LASTEXITCODE -ne 0) { throw 'Could not install the locked dependencies.' }
Write-Host ('Python environment ready: ' + $pythonExe)
