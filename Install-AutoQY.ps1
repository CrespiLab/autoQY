[CmdletBinding()]
param(
    [string]$EnvironmentName = "autoqy-core",
    [string]$RepositoryUrl = "https://github.com/CrespiLab/autoQY.git",
    [string]$Branch = "main",
    [switch]$CheckOnly
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"
$InstallerDirectory = Split-Path -Parent $PSCommandPath
$CurrentDirectory = (Get-Location).Path
$InstallTimer = [System.Diagnostics.Stopwatch]::StartNew()

function Get-ElapsedText {
    return $InstallTimer.Elapsed.ToString("hh\:mm\:ss")
}

function Write-Step {
    param([string]$Message)
    Write-Host ""
    Write-Host "[$(Get-ElapsedText)] ==> $Message" -ForegroundColor Cyan
}

function Read-Confirmation {
    param(
        [string]$Message,
        [bool]$DefaultYes = $false
    )
    $suffix = if ($DefaultYes) { "[Y/n]" } else { "[y/N]" }
    while ($true) {
        $answer = (Read-Host "$Message $suffix").Trim()
        if (-not $answer) { return $DefaultYes }
        if ($answer -match "^(?i:y|yes)$") { return $true }
        if ($answer -match "^(?i:n|no)$") { return $false }
        Write-Host "Please answer yes or no."
    }
}

function Select-InstallDirectory {
    param([string]$CurrentPath)

    $currentFullPath = [System.IO.Path]::GetFullPath($CurrentPath)
    Write-Host "Current folder: $currentFullPath"
    if (Read-Confirmation "Use the current folder as the AutoQY installation folder?" $true) {
        return $currentFullPath
    }

    while ($true) {
        $enteredPath = (Read-Host "Enter the full path to the installation folder").Trim().Trim('"')
        if (-not $enteredPath) {
            Write-Host "Please enter a folder path."
            continue
        }

        $enteredPath = [Environment]::ExpandEnvironmentVariables($enteredPath)
        if (-not [System.IO.Path]::IsPathRooted($enteredPath)) {
            $enteredPath = Join-Path $currentFullPath $enteredPath
        }

        try {
            $selectedPath = [System.IO.Path]::GetFullPath($enteredPath)
        }
        catch {
            Write-Host "The folder path is not valid: $enteredPath" -ForegroundColor Yellow
            continue
        }

        if (Test-Path -LiteralPath $selectedPath) {
            if (-not (Get-Item -LiteralPath $selectedPath).PSIsContainer) {
                Write-Host "The selected path is not a folder: $selectedPath" -ForegroundColor Yellow
                continue
            }
            return (Resolve-Path -LiteralPath $selectedPath).Path
        }

        if (Read-Confirmation "The folder does not exist. Create '$selectedPath'?" $true) {
            New-Item -ItemType Directory -Path $selectedPath -Force | Out-Null
            return (Resolve-Path -LiteralPath $selectedPath).Path
        }
    }
}

function Get-CondaCommand {
    $available = Get-Command conda -ErrorAction SilentlyContinue
    if ($available) {
        if ($available.CommandType -eq "Application") { return $available.Source }
        return "conda"
    }
    $candidates = @(
        (Join-Path $env:USERPROFILE "anaconda3\Scripts\conda.exe"),
        (Join-Path $env:USERPROFILE "miniconda3\Scripts\conda.exe"),
        (Join-Path $env:LOCALAPPDATA "anaconda3\Scripts\conda.exe"),
        (Join-Path $env:LOCALAPPDATA "miniconda3\Scripts\conda.exe")
    )
    return $candidates | Where-Object { Test-Path -LiteralPath $_ } | Select-Object -First 1
}

function Invoke-Checked {
    param(
        [string]$Command,
        [string[]]$Arguments,
        [string]$Activity
    )
    $commandTimer = [System.Diagnostics.Stopwatch]::StartNew()
    Write-Host "   $Activity"
    Write-Host "> $Command $($Arguments -join ' ')" -ForegroundColor DarkGray
    & $Command @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "Command failed with exit code $LASTEXITCODE`: $Command"
    }
    Write-Host "   Completed in $($commandTimer.Elapsed.ToString('hh\:mm\:ss'))." -ForegroundColor Green
}

function Get-CondaEnvironmentPath {
    param(
        [string]$CondaCommand,
        [string]$Name
    )
    $jsonText = & $CondaCommand env list --json
    if ($LASTEXITCODE -ne 0) { throw "Conda could not list its environments." }
    $environmentList = (($jsonText -join "`n") | ConvertFrom-Json).envs
    return $environmentList | Where-Object {
        (Split-Path -Leaf ([System.IO.Path]::GetFullPath($_).TrimEnd("\"))) -eq $Name
    } | Select-Object -First 1
}

function Test-AutoQYProject {
    param([string]$Path)
    return (Test-Path -LiteralPath (Join-Path $Path "pyproject.toml")) -and
           (Test-Path -LiteralPath (Join-Path $Path "autoqy_core"))
}

function Get-EnvironmentGit {
    param([string]$EnvironmentPath)
    $candidates = @(
        (Join-Path $EnvironmentPath "Library\bin\git.exe"),
        (Join-Path $EnvironmentPath "Scripts\git.exe"),
        (Join-Path $EnvironmentPath "git.exe")
    )
    return $candidates | Where-Object { Test-Path -LiteralPath $_ } | Select-Object -First 1
}

function Get-CondaHookPath {
    param([string]$CondaCommand)
    $condaBaseText = & $CondaCommand info --base
    if ($LASTEXITCODE -ne 0 -or -not $condaBaseText) {
        throw "Conda could not locate its activation script."
    }
    $condaBase = ($condaBaseText | Select-Object -Last 1).Trim()
    $condaHook = Join-Path $condaBase "shell\condabin\conda-hook.ps1"
    if (-not (Test-Path -LiteralPath $condaHook)) {
        throw "Conda activation hook not found: $condaHook"
    }
    return $condaHook
}

function Activate-AutoQYEnvironment {
    param(
        [string]$CondaCommand,
        [string]$Name,
        [string]$ExpectedPath
    )
    $condaHook = Get-CondaHookPath -CondaCommand $CondaCommand
    . $condaHook
    conda activate $Name
    if (-not $env:CONDA_PREFIX) { throw "Conda could not activate '$Name'." }
    $activatedPath = [System.IO.Path]::GetFullPath($env:CONDA_PREFIX).TrimEnd("\")
    $expectedFull = [System.IO.Path]::GetFullPath($ExpectedPath).TrimEnd("\")
    if (-not [string]::Equals($activatedPath, $expectedFull, [StringComparison]::OrdinalIgnoreCase)) {
        throw "Conda activated an unexpected environment: $activatedPath"
    }
}

function Write-LauncherFiles {
    param(
        [string]$ProjectRoot,
        [string]$EnvironmentPath,
        [string]$CondaCommand
    )
    $launcherDirectory = Join-Path $ProjectRoot ".autoqy-launchers"
    New-Item -ItemType Directory -Path $launcherDirectory -Force | Out-Null

    $condaHook = Get-CondaHookPath -CondaCommand $CondaCommand

    $terminalScript = Join-Path $launcherDirectory "Open-AutoQY-Terminal.ps1"
    $terminalContent = @"
. '$($condaHook.Replace("'", "''"))'
conda activate '$($EnvironmentName.Replace("'", "''"))'
Set-Location -LiteralPath '$($ProjectRoot.Replace("'", "''"))'
Write-Host 'AutoQY environment activated. Try: autoqy-core --help' -ForegroundColor Green
"@
    Set-Content -LiteralPath $terminalScript -Value $terminalContent -Encoding UTF8

    $analysisScript = Join-Path $launcherDirectory "Analyze-AutoQY-JSON.ps1"
    $pythonCommand = Join-Path $EnvironmentPath "python.exe"
    $analysisContent = @"
param([string]`$ConfigPath)
`$ErrorActionPreference = 'Stop'
if (-not `$ConfigPath) { `$ConfigPath = Read-Host 'Enter or paste the path to analysis.json' }
try {
    `$absoluteConfig = (Resolve-Path -LiteralPath `$ConfigPath).Path
    if ([IO.Path]::GetExtension(`$absoluteConfig) -ne '.json') { throw 'The input must be a JSON file.' }
    Write-Host "Configuration: `$absoluteConfig" -ForegroundColor Cyan
    & '$($pythonCommand.Replace("'", "''"))' -m autoqy_core validate `$absoluteConfig
    if (`$LASTEXITCODE -ne 0) { throw 'Configuration validation failed.' }
    `$answer = (Read-Host 'Validation succeeded. Run the analysis now? [Y/n]').Trim()
    if (-not `$answer -or `$answer -match '^(?i:y|yes)$') {
        & '$($pythonCommand.Replace("'", "''"))' -m autoqy_core run `$absoluteConfig
        if (`$LASTEXITCODE -ne 0) { throw 'AutoQY analysis failed.' }
    }
}
catch { Write-Host "Error: `$(`$_.Exception.Message)" -ForegroundColor Red }
Read-Host 'Press Enter to close'
"@
    Set-Content -LiteralPath $analysisScript -Value $analysisContent -Encoding UTF8

    $desktop = [Environment]::GetFolderPath("Desktop")
    $shortcutDirectory = Join-Path $desktop "AutoQY"
    New-Item -ItemType Directory -Path $shortcutDirectory -Force | Out-Null
    $shell = New-Object -ComObject WScript.Shell
    $coreCommand = Join-Path $EnvironmentPath "Scripts\autoqy-core.exe"
    $iconDirectory = Join-Path $ProjectRoot "autoqy_core\assets\icons"
    $guiIcon = Join-Path $iconDirectory "power-gui.ico"
    $smootherIcon = Join-Path $iconDirectory "spectral-smoother.ico"
    $terminalIcon = Join-Path $iconDirectory "terminal.ico"
    $jsonIcon = Join-Path $iconDirectory "analyze-json.ico"
    foreach ($iconPath in @($guiIcon, $smootherIcon, $terminalIcon, $jsonIcon)) {
        if (-not (Test-Path -LiteralPath $iconPath)) { throw "Desktop icon not found: $iconPath" }
    }

    $guiShortcutPath = Join-Path $shortcutDirectory "AutoQY Power GUI.lnk"
    $guiShortcut = $shell.CreateShortcut($guiShortcutPath)
    $guiShortcut.TargetPath = $coreCommand
    $guiShortcut.Arguments = "power-gui"
    $guiShortcut.WorkingDirectory = $ProjectRoot
    $guiShortcut.Description = "Open the AutoQY power-treatment GUI"
    $guiShortcut.IconLocation = "$guiIcon,0"
    $guiShortcut.Save()

    $legacySmootherShortcut = Join-Path $shortcutDirectory "AutoQY Spectral Smoother.lnk"
    if (Test-Path -LiteralPath $legacySmootherShortcut) {
        Remove-Item -LiteralPath $legacySmootherShortcut -Force
    }
    $smootherShortcutPath = Join-Path $shortcutDirectory "AutoQY Spectral Treatment.lnk"
    $smootherShortcut = $shell.CreateShortcut($smootherShortcutPath)
    $smootherShortcut.TargetPath = $coreCommand
    $smootherShortcut.Arguments = "smoother-gui"
    $smootherShortcut.WorkingDirectory = $ProjectRoot
    $smootherShortcut.Description = "Open the AutoQY spectral treatment GUI"
    $smootherShortcut.IconLocation = "$smootherIcon,0"
    $smootherShortcut.Save()

    $terminalShortcutPath = Join-Path $shortcutDirectory "AutoQY Terminal.lnk"
    $terminalShortcut = $shell.CreateShortcut($terminalShortcutPath)
    $terminalShortcut.TargetPath = "powershell.exe"
    $terminalShortcut.Arguments = "-NoExit -ExecutionPolicy Bypass -File `"$terminalScript`""
    $terminalShortcut.WorkingDirectory = $ProjectRoot
    $terminalShortcut.Description = "Open PowerShell with the AutoQY Conda environment"
    $terminalShortcut.IconLocation = "$terminalIcon,0"
    $terminalShortcut.Save()

    $dropCommand = Join-Path $launcherDirectory "Analyze-AutoQY-JSON.cmd"
    $dropContent = @"
@echo off
powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -File "$analysisScript" "%~1"
"@
    Set-Content -LiteralPath $dropCommand -Value $dropContent -Encoding Ascii

    $jsonShortcutPath = Join-Path $shortcutDirectory "AutoQY Analyze JSON.lnk"
    $jsonShortcut = $shell.CreateShortcut($jsonShortcutPath)
    $jsonShortcut.TargetPath = $dropCommand
    $jsonShortcut.WorkingDirectory = $ProjectRoot
    $jsonShortcut.Description = "Drag and drop an AutoQY analysis JSON file"
    $jsonShortcut.IconLocation = "$jsonIcon,0"
    $jsonShortcut.Save()

    return @($guiShortcutPath, $smootherShortcutPath, $terminalShortcutPath, $jsonShortcutPath)
}

try {
    Write-Host "AutoQY Conda bootstrap installer" -ForegroundColor Blue
    Write-Host "Installer file folder: $InstallerDirectory"

    $InstallDirectory = if ($CheckOnly) {
        [System.IO.Path]::GetFullPath($CurrentDirectory)
    }
    else {
        Select-InstallDirectory -CurrentPath $CurrentDirectory
    }
    $ClonePath = Join-Path $InstallDirectory "AutoQY-Core"

    Write-Host "Installation folder: $InstallDirectory"
    Write-Host "Clone destination: $ClonePath"
    Write-Host "Conda environment: $EnvironmentName"

    if ($EnvironmentName -notmatch "^[A-Za-z0-9_-]+$") {
        throw "The Conda environment name may contain only letters, numbers, underscores, and hyphens."
    }
    $condaCommand = Get-CondaCommand
    if (-not $condaCommand) {
        throw "Conda was not found. Install Anaconda/Miniconda or run this BAT from Anaconda PowerShell Prompt."
    }

    $projectInInstallFolder = Test-AutoQYProject -Path $InstallDirectory
    $projectRoot = if ($projectInInstallFolder) { $InstallDirectory } else { $ClonePath }
    $environmentPath = Get-CondaEnvironmentPath -CondaCommand $condaCommand -Name $EnvironmentName

    if ($CheckOnly) {
        Write-Step "Planned installation"
        Write-Host "Conda: $condaCommand"
        Write-Host "Environment: $(if ($environmentPath) { "ask to remove $environmentPath, then recreate it" } else { "create $EnvironmentName" })"
        Write-Host "Source: $(if ($projectInInstallFolder) { "use existing checkout $projectRoot" } else { "clone $RepositoryUrl into $projectRoot" })"
        Write-Host "Package: editable install with both browser GUIs"
        Write-Host "Desktop folder AutoQY: Power GUI, Spectral Treatment, activated terminal, and JSON runner"
        Write-Host "No files or environments were changed." -ForegroundColor Green
        exit 0
    }

    if ($environmentPath) {
        Write-Step "Existing Conda environment detected"
        Write-Host "   $environmentPath" -ForegroundColor Yellow
        Write-Host "   Removing it deletes every package and file stored in that environment." -ForegroundColor Yellow
        if (-not (Read-Confirmation "Delete '$EnvironmentName' and recreate it cleanly?" $false)) {
            Write-Host "Installation cancelled; the existing environment was left unchanged."
            exit 0
        }
        if ($env:CONDA_DEFAULT_ENV -eq $EnvironmentName) {
            $condaHook = Get-CondaHookPath -CondaCommand $condaCommand
            . $condaHook
            conda deactivate
        }
        Invoke-Checked -Command $condaCommand -Arguments @(
            "env", "remove", "--name", $EnvironmentName, "--yes"
        ) -Activity "Removing the existing AutoQY Conda environment."
        $environmentPath = Get-CondaEnvironmentPath -CondaCommand $condaCommand -Name $EnvironmentName
        if ($environmentPath) {
            throw "Conda still reports the environment after removal: $environmentPath"
        }
        Write-Host "   Existing environment removed cleanly." -ForegroundColor Green
    }

    Write-Step "Creating Conda environment '$EnvironmentName'"
    Invoke-Checked -Command $condaCommand -Arguments @(
        "create", "--name", $EnvironmentName, "python=3.12", "pip", "--yes"
    ) -Activity "Installing Python 3.12 and pip. This may take several minutes."
    Invoke-Checked -Command $condaCommand -Arguments @(
        "install", "--name", $EnvironmentName, "git", "--yes"
    ) -Activity "Installing Git into the AutoQY environment."
    $environmentPath = Get-CondaEnvironmentPath -CondaCommand $condaCommand -Name $EnvironmentName
    if (-not $environmentPath) { throw "The new Conda environment could not be located." }

    Write-Step "Activating Conda environment '$EnvironmentName'"
    Activate-AutoQYEnvironment -CondaCommand $condaCommand `
        -Name $EnvironmentName -ExpectedPath $environmentPath
    Write-Host "   Activated: $env:CONDA_PREFIX" -ForegroundColor Green

    if (-not $projectInInstallFolder) {
        if (Test-Path -LiteralPath $ClonePath) {
            if (-not (Test-AutoQYProject -Path $ClonePath)) {
                throw "The clone destination already exists but is not an AutoQY project: $ClonePath"
            }
            if (-not (Read-Confirmation "An AutoQY checkout already exists at the destination. Use it without overwriting local work?" $true)) {
                Write-Host "Installation cancelled."
                exit 0
            }
        }
        else {
            Write-Step "Cloning AutoQY into the selected installation folder"
            $gitCommand = Get-EnvironmentGit -EnvironmentPath $environmentPath
            if (-not $gitCommand) { throw "Git was installed but git.exe could not be located in the environment." }
            Invoke-Checked -Command $gitCommand -Arguments @(
                "clone", "--branch", $Branch, "--single-branch", $RepositoryUrl, $ClonePath
            ) -Activity "Downloading branch '$Branch' from GitHub."
        }
    }

    $projectRoot = if ($projectInInstallFolder) { $InstallDirectory } else { $ClonePath }
    $environmentPython = Join-Path $environmentPath "python.exe"
    if (-not (Test-Path -LiteralPath $environmentPython)) {
        throw "python.exe was not found in the AutoQY environment."
    }

    Write-Step "Installing AutoQY and the power GUI"
    Push-Location $projectRoot
    try {
        Invoke-Checked -Command $environmentPython -Arguments @(
            "-m", "pip", "install", "--disable-pip-version-check", "--progress-bar", "off",
            "--editable", ".[power-gui]"
        ) -Activity "Running: python -m pip install -e `".[power-gui]`""
    }
    finally {
        Pop-Location
    }

    Write-Step "Validating the installation"
    $analysisConfig = Join-Path $projectRoot "ExampleData\Example-2_AB_455nm-100mA\analysis.json"
    Invoke-Checked -Command $environmentPython -Arguments @(
        "-m", "autoqy_core", "validate", $analysisConfig
    ) -Activity "Validating the bundled 455 nm configuration."

    Write-Step "Creating desktop launchers"
    $desktopFiles = Write-LauncherFiles -ProjectRoot $projectRoot `
        -EnvironmentPath $environmentPath -CondaCommand $condaCommand
    $desktopFiles | ForEach-Object { Write-Host "   $_" -ForegroundColor Green }

    Write-Step "Installation complete"
    Write-Host "Repository: $projectRoot" -ForegroundColor Green
    Write-Host "Environment: $environmentPath"
    Write-Host "Elapsed time: $(Get-ElapsedText)"
    Write-Host "Run in any terminal: conda run --name $EnvironmentName autoqy-core --help"
}
catch {
    Write-Host ""
    Write-Host "Installation failed: $($_.Exception.Message)" -ForegroundColor Red
    exit 1
}

