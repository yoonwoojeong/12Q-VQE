# PowerShell One-Command GitHub & Qiskit Ecosystem Submission Script
# Run from PowerShell:
#   powershell -ExecutionPolicy Bypass -File C:\Users\sshs4\12Q-VQE\scripts\submit_ecosystem.ps1

$ErrorActionPreference = "Stop"
Set-Location -Path (Split-Path -Parent $PSScriptRoot)

Write-Host "==> 1. Pushing latest main and tag v0.2.0 to origin..." -ForegroundColor Cyan
git push origin main
git push origin v0.2.0 --force

Write-Host "==> 2. Retrieving GitHub token from Git Credential Manager..." -ForegroundColor Cyan
$credInput = "protocol=https`nhost=github.com`n`n"
$credOutput = $credInput | & "C:\Program Files\Git\mingw64\bin\git-credential-manager.exe" get
$token = ($credOutput | Where-Object { $_ -match "^password=" }) -replace "^password=", ""
if (-not $token) {
    throw "Could not retrieve GitHub token from Git Credential Manager."
}

$headers = @{
    "Authorization" = "Bearer $token"
    "Accept"        = "application/vnd.github+json"
    "X-GitHub-Api-Version" = "2022-11-28"
    "User-Agent"    = "qchem-vqe-submitter"
}

Write-Host "==> 3. Setting GitHub repository description and topics (required by Qiskit Ecosystem CI [G14]/[G15])..." -ForegroundColor Cyan
$repoPatch = @{
    description = "C17 ab initio quantum chemistry, adjoint L-BFGS/ADAPT-VQE, and Sample-Based Quantum Diagonalization (SQD) engine with IBM Quantum support"
} | ConvertTo-Json
Invoke-RestMethod -Uri "https://api.github.com/repos/yoonwoojeong/qchem-vqe" -Method Patch -Headers $headers -Body $repoPatch -ContentType "application/json" | Out-Null

$topicsPut = @{
    names = @("qiskit", "quantum-chemistry", "vqe", "adapt-vqe", "sqd", "openqasm3", "ibm-quantum", "c17")
} | ConvertTo-Json
Invoke-RestMethod -Uri "https://api.github.com/repos/yoonwoojeong/qchem-vqe/topics" -Method Put -Headers $headers -Body $topicsPut -ContentType "application/json" | Out-Null
Write-Host "    Repository description and 'qiskit' topic configured." -ForegroundColor Green

Write-Host "==> 4. Opening Qiskit Ecosystem submission issue on Qiskit/ecosystem..." -ForegroundColor Cyan
$issueBody = @"
### Project name

qchem-vqe

### Description

``````text
C17 quantum chemistry engine (RHF, FCIDUMP, adjoint VQE, SQD) with OpenQASM 3 and IBM Quantum Estimator V2 support.
``````

### Contact email

154740938+yoonwoojeong@users.noreply.github.com

### Category

Application

### Labels

chemistry, Hamiltonian simulation, openqasm, circuit optimization, HPC

### Interface/API

Python, C, Command-line interface (CLI)

### Stability and support expecations

production-ready

### GitHub repository

https://github.com/yoonwoojeong/qchem-vqe

### Home page

_No response_

### Documentation

_No response_

### Packages

``````text
https://github.com/yoonwoojeong/qchem-vqe/blob/main/pyproject.toml
``````

### Reference paper

https://github.com/yoonwoojeong/qchem-vqe/blob/main/paper/paper.md

### Code of Conduct

- [x] I agree to follow the [Qiskit Code of Conduct](https://qisk.it/coc)
"@

$issuePayload = @{
    title  = "[Submission]: qchem-vqe"
    body   = $issueBody
} | ConvertTo-Json -Depth 5

$issueResp = Invoke-RestMethod -Uri "https://api.github.com/repos/Qiskit/ecosystem/issues" -Method Post -Headers $headers -Body $issuePayload -ContentType "application/json; charset=utf-8"
Write-Host "    Submitted Qiskit Ecosystem issue: $($issueResp.html_url)" -ForegroundColor Green
