[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$ProjectRoot,
    [Parameter(Mandatory = $true)]
    [string]$EnvironmentName,
    [Parameter(Mandatory = $true)]
    [string]$CondaCommand,
    [Parameter(Mandatory = $true)]
    [string]$ShortcutDirectory,
    [string]$EnvironmentPath = ""
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"
$uninstallExitCode = 0

function Read-UninstallConfirmation {
    param([string]$Message)
    while ($true) {
        $answer = (Read-Host "$Message [y/N]").Trim()
        if (-not $answer -or $answer -match "^(?i:n|no)$") { return $false }
        if ($answer -match "^(?i:y|yes)$") { return $true }
        Write-Host "Please answer yes or no."
    }
}

function Test-UninstallPathWithin {
    param([string]$Candidate, [string]$Directory)
    if (-not $Candidate -or -not $Directory) { return $false }
    $directoryFullPath = [IO.Path]::GetFullPath($Directory).TrimEnd("\", "/")
    try { $candidateFullPath = [IO.Path]::GetFullPath($Candidate) }
    catch { return $false }
    return ($candidateFullPath.StartsWith(
        $directoryFullPath + "\", [StringComparison]::OrdinalIgnoreCase
    ) -or [string]::Equals(
        $candidateFullPath, $directoryFullPath, [StringComparison]::OrdinalIgnoreCase
    ))
}

function Get-UninstallEnvironmentPath {
    param([string]$Command, [string]$Name)
    $environmentJson = & $Command @("env", "list", "--json")
    if ($LASTEXITCODE -ne 0) { throw "Conda could not list its environments." }
    $environmentList = ($environmentJson -join "`n") | ConvertFrom-Json
    $environmentMatches = @($environmentList.envs | Where-Object {
        [string]::Equals((Split-Path -Leaf $_), $Name, [StringComparison]::OrdinalIgnoreCase)
    })
    if ($environmentMatches.Count -gt 1) { throw "More than one Conda environment is named '$Name'." }
    if ($environmentMatches.Count -eq 1) { return [IO.Path]::GetFullPath($environmentMatches[0]).TrimEnd("\") }
    return $null
}

function Get-UninstallBlockingProcesses {
    param([string]$Installation, [string]$Environment, [bool]$RemovingEnvironment)
    # Never kill processes: users must be able to save before closing a GUI.
    # Fail closed if Windows cannot check processes; do not guess that it is safe.
    $processes = Get-CimInstance -ClassName Win32_Process -ErrorAction Stop
    $installationPattern = [regex]::Escape($Installation.TrimEnd("\", "/")) + '(?=[\\/"\s]|$)'
    foreach ($process in $processes) {
        if ($process.ProcessId -eq $PID) { continue }
        $executable = [string]$process.ExecutablePath
        $commandLine = [string]$process.CommandLine
        $name = [string]$process.Name
        $usesInstallation = (Test-UninstallPathWithin $executable $Installation) -or
                            ($commandLine -match $installationPattern)
        $usesEnvironment = (Test-UninstallPathWithin $executable $Environment)
        $isAutoQY = ($name -match '^autoqy-(core|analysis-gui|smoother-gui|power-gui)\.exe$') -or
                     ($commandLine -match '\b(autoqy_core|autoqy-(core|analysis-gui|smoother-gui|power-gui))\b')
        if ($usesInstallation -or
            ($usesEnvironment -and ($RemovingEnvironment -or $isAutoQY)) -or
            ($isAutoQY -and (-not $Environment -or -not $executable))) {
            $process
        }
    }
}

function Wait-UninstallProcessesClosed {
    param([string]$Installation, [string]$Environment, [bool]$RemovingEnvironment)
    while ($true) {
        $blocking = @(Get-UninstallBlockingProcesses $Installation $Environment $RemovingEnvironment)
        if ($blocking.Count -eq 0) { return $true }
        Write-Host "AutoQY or a program using its files is still running." -ForegroundColor Yellow
        Write-Host "Save your work and close all AutoQY windows. Nothing will be deleted while they are open."
        if ($RemovingEnvironment) {
            Write-Host "Also close any terminal or Python program using the AutoQY Conda environment."
        }
        foreach ($process in $blocking) {
            Write-Host "   $($process.Name) (process $($process.ProcessId))"
        }
        $answer = (Read-Host "Press Enter to check again, or type cancel").Trim()
        if ($answer -match '^(?i:c|cancel)$') { return $false }
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
    # New shortcuts supply the exact environment. Older shortcuts still work:
    # resolve their environment before checking running GUI/Python processes.
    if (-not $EnvironmentPath) {
        try {
            $EnvironmentPath = Get-UninstallEnvironmentPath $CondaCommand $EnvironmentName
        }
        catch {
            if ($removeEnvironment) { throw }
            Write-Host "Conda could not be checked; checking running AutoQY programs instead." -ForegroundColor Yellow
        }
    }
    if ($EnvironmentPath) {
        $EnvironmentPath = [IO.Path]::GetFullPath($EnvironmentPath).TrimEnd("\")
        if ([string]::Equals($EnvironmentPath, [IO.Path]::GetPathRoot($EnvironmentPath).TrimEnd("\"),
                            [StringComparison]::OrdinalIgnoreCase)) {
            throw "Refusing to use a drive root as the Conda environment."
        }
    }
    if (-not (Wait-UninstallProcessesClosed $projectFullPath $EnvironmentPath $removeEnvironment)) {
        Write-Host "Uninstall cancelled; no files or environments were changed."
        Read-Host "Press Enter to close"
        exit 0
    }
    if ($removeEnvironment -and $EnvironmentPath -and
        (Test-Path -LiteralPath $EnvironmentPath -PathType Container)) {
        Write-Host "Removing Conda environment '$EnvironmentName'..." -ForegroundColor Cyan
        try {
            & $CondaCommand @("env", "remove", "--prefix", $EnvironmentPath, "--yes")
            if ($LASTEXITCODE -ne 0) {
                throw "Conda exited with code $LASTEXITCODE."
            }
            Write-Host "Conda environment removed." -ForegroundColor Green
        }
        catch {
            throw "The Conda environment could not be removed ($($_.Exception.Message)). The installation folder and shortcuts were kept. Close programs using AutoQY and try uninstall again."
        }
    }
    elseif ($removeEnvironment) {
        Write-Host "Conda environment '$EnvironmentName' is already absent; nothing to remove." -ForegroundColor Yellow
    }
    else {
        Write-Host "Conda environment '$EnvironmentName' was kept." -ForegroundColor Yellow
    }

    # Recheck after Conda finishes: a GUI may have been launched in the meantime.
    if (-not (Wait-UninstallProcessesClosed $projectFullPath $EnvironmentPath $false)) {
        Write-Host "Uninstall stopped. The installation folder and shortcuts were kept."
        Read-Host "Press Enter to close"
        exit 0
    }
    Write-Host "Removing $projectFullPath..." -ForegroundColor Cyan
    if (Test-Path -LiteralPath $projectFullPath) {
        Remove-Item -LiteralPath $projectFullPath -Recurse -Force
    }
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
    $uninstallExitCode = 1
    Write-Host "Uninstall failed: $($_.Exception.Message)" -ForegroundColor Red
}

Read-Host "Press Enter to close"
exit $uninstallExitCode
