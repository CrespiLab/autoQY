[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$ProjectRoot,
    [Parameter(Mandatory = $true)]
    [string]$EnvironmentName,
    [Parameter(Mandatory = $true)]
    [string]$CondaCommand,
    [Parameter(Mandatory = $true)]
    [string]$ShortcutDirectory
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

function Read-UninstallConfirmation {
    param([string]$Message)
    while ($true) {
        $answer = (Read-Host "$Message [y/N]").Trim()
        if (-not $answer -or $answer -match "^(?i:n|no)$") { return $false }
        if ($answer -match "^(?i:y|yes)$") { return $true }
        Write-Host "Please answer yes or no."
    }
}

try {
    $projectFullPath = [System.IO.Path]::GetFullPath($ProjectRoot).TrimEnd("\")
    $driveRoot = [System.IO.Path]::GetPathRoot($projectFullPath).TrimEnd("\")
    if ([string]::Equals(
        $projectFullPath, $driveRoot, [StringComparison]::OrdinalIgnoreCase
    )) {
        throw "Refusing to remove a drive root."
    }
    if (-not (Test-Path -LiteralPath (Join-Path $projectFullPath "pyproject.toml")) -or
        -not (Test-Path -LiteralPath (Join-Path $projectFullPath "autoqy_core"))) {
        throw "The recorded installation folder no longer looks like an AutoQY installation: $projectFullPath"
    }

    Write-Host "AutoQY uninstaller" -ForegroundColor Cyan
    Write-Host "Installation folder: $projectFullPath"
    Write-Host "Conda environment: $EnvironmentName"
    if (-not (Read-UninstallConfirmation "Are you sure you want to remove the AutoQY installation folder?")) {
        Write-Host "Uninstall cancelled; no files or environments were changed."
        Read-Host "Press Enter to close"
        exit 0
    }

    $removeEnvironment = Read-UninstallConfirmation "Also remove the Conda environment '$EnvironmentName'?"
    if ($removeEnvironment) {
        Write-Host "Removing Conda environment '$EnvironmentName'..." -ForegroundColor Cyan
        try {
            & $CondaCommand @("env", "remove", "--name", $EnvironmentName, "--yes")
            if ($LASTEXITCODE -ne 0) {
                throw "Conda exited with code $LASTEXITCODE."
            }
            Write-Host "Conda environment removed." -ForegroundColor Green
        }
        catch {
            Write-Host "The Conda environment could not be removed ($($_.Exception.Message)); continuing with the installation folder." -ForegroundColor Yellow
        }
    }
    else {
        Write-Host "Conda environment '$EnvironmentName' was kept." -ForegroundColor Yellow
    }

    Write-Host "Removing $projectFullPath..." -ForegroundColor Cyan
    Remove-Item -LiteralPath $projectFullPath -Recurse -Force
    Write-Host "AutoQY installation folder removed." -ForegroundColor Green

    $shortcutFullPath = [System.IO.Path]::GetFullPath($ShortcutDirectory).TrimEnd("\")
    $desktopPath = [System.IO.Path]::GetFullPath(
        [Environment]::GetFolderPath("Desktop")
    ).TrimEnd("\")
    $shortcutParent = Split-Path -Parent $shortcutFullPath
    if ((Split-Path -Leaf $shortcutFullPath) -eq "AutoQY" -and
        [string]::Equals($shortcutParent, $desktopPath, [StringComparison]::OrdinalIgnoreCase) -and
        (Test-Path -LiteralPath $shortcutFullPath)) {
        Remove-Item -LiteralPath $shortcutFullPath -Recurse -Force
        Write-Host "Desktop shortcuts removed." -ForegroundColor Green
    }
    else {
        Write-Host "Desktop shortcut folder was not removed because its recorded path was unexpected: $shortcutFullPath" -ForegroundColor Yellow
    }
}
catch {
    Write-Host "Uninstall failed: $($_.Exception.Message)" -ForegroundColor Red
}

Read-Host "Press Enter to close"
