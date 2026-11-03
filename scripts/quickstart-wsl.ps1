<#
.SYNOPSIS
  Runs the qchem-vqe quickstart on Windows by running it inside WSL (Ubuntu).

.DESCRIPTION
  qchem-vqe uses POSIX APIs and Linux packages, so on Windows it is built and run
  inside the Windows Subsystem for Linux. This script:
    1. checks that WSL and the Ubuntu distribution are installed,
    2. converts this repository's path to its WSL path,
    3. runs the same one-line quickstart as the README inside Ubuntu.
  sudo may ask for your Ubuntu password in the console.

.PARAMETER Distro
  WSL distribution to use. Default: Ubuntu.

.EXAMPLE
  powershell -ExecutionPolicy Bypass -File scripts\quickstart-wsl.ps1
#>
[CmdletBinding()]
param(
    [string]$Distro = 'Ubuntu'
)

# Not 'Stop': Windows PowerShell 5.1 turns any stderr line of a native command
# (wsl.exe warnings included) into a terminating error. Exit codes are checked below.
$ErrorActionPreference = 'Continue'

function Stop-With([string]$Message) {
    Write-Host "error: $Message" -ForegroundColor Red
    exit 1
}

# 1. WSL present?
if (-not (Get-Command wsl.exe -ErrorAction SilentlyContinue)) {
    Stop-With "WSL is not installed. In an elevated PowerShell run 'wsl --install -d $Distro', reboot, then run this script again."
}

# 2. Distribution installed? (wsl --list prints UTF-16 with NUL characters, so strip them)
$names = @(wsl.exe --list --quiet 2>$null | ForEach-Object { ($_ -replace "`0", '').Trim() } | Where-Object { $_ })
if ($names -notcontains $Distro) {
    Stop-With "WSL distribution '$Distro' is not installed. In an elevated PowerShell run 'wsl --install -d $Distro', reboot, then run this script again."
}

# 3. Repository path as seen from WSL (the script lives in <repo>\scripts)
$repoWin = Split-Path -Parent $PSScriptRoot
# --exec runs the program directly instead of through the Linux shell, which would
# strip backslashes and re-split arguments. Forward slashes are a second safeguard.
$repoWinFwd = $repoWin -replace '\\', '/'
$repoWsl = (wsl.exe -d $Distro --exec wslpath -a -u $repoWinFwd 2>$null)
if ($repoWsl) { $repoWsl = "$repoWsl".Trim() }
if ($LASTEXITCODE -ne 0 -or -not $repoWsl) {
    Stop-With "could not convert '$repoWin' to a WSL path."
}
if ($repoWsl -like '/mnt/*') {
    Write-Warning "The repository is on the Windows drive ($repoWsl). Builds are slower there; cloning into the Linux home of $Distro is faster."
}

Write-Host "Running the quickstart in $Distro at $repoWsl" -ForegroundColor Cyan

# 4. The same one-line quickstart as the README (Ubuntu/Debian packages, libcint, build, tests)
$quickstart = 'sudo apt-get update && sudo apt-get install -y build-essential curl unzip zstd libopenblas-dev liblapacke-dev libnlopt-dev libcjson-dev libcurl4-openssl-dev && make fetch-libcint && make && make test'
wsl.exe -d $Distro --cd $repoWsl --exec bash -lc $quickstart
if ($LASTEXITCODE -ne 0) {
    Stop-With "the quickstart failed inside $Distro (exit code $LASTEXITCODE). See the output above."
}

Write-Host ""
Write-Host "Done. Run the program inside $Distro with:" -ForegroundColor Green
Write-Host "  wsl -d $Distro --cd $repoWsl -- ./build/qchem-vqe rhf"
