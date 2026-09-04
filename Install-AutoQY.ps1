[CmdletBinding()]
param(
    [string]$EnvironmentName = "autoqy-core",
    [string]$RepositoryUrl = "https://github.com/CrespiLab/autoQY.git",
    [string]$Branch = "main",
    [switch]$CheckOnly
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"
$CurrentDirectory = (Get-Location).Path
$InstallTimer = [System.Diagnostics.Stopwatch]::StartNew()

function Get-ElapsedText {
    return $InstallTimer.Elapsed.ToString("hh\:mm\:ss")
}

function Write-Step {
    param([string]$Message)
    Write-Host ""
    Write-Host "[$(Get-ElapsedText)] $Message" -ForegroundColor Cyan
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

function Invoke-CapturedCommand {
    param(
        [string]$Command,
        [string[]]$Arguments
    )
    $previousErrorPreference = $ErrorActionPreference
    $output = @()
    $exitCode = 1
    try {
        $ErrorActionPreference = "Continue"
        $output = @(& $Command @Arguments 2>&1)
        $commandSucceeded = $?
        $exitCode = if ($commandSucceeded) {
            0
        }
        elseif ($LASTEXITCODE) {
            $LASTEXITCODE
        }
        else {
            1
        }
    }
    finally {
        $ErrorActionPreference = $previousErrorPreference
    }
    return [PSCustomObject]@{
        Output = $output
        ExitCode = $exitCode
    }
}

function Write-CommandFailure {
    param(
        [string]$Command,
        [string[]]$Arguments,
        [object[]]$Output
    )
    Write-Host ""
    Write-Host "Technical details from the failed step:" -ForegroundColor Yellow
    Write-Host "> $Command $($Arguments -join ' ')" -ForegroundColor DarkGray
    foreach ($line in $Output) {
        Write-Host "  $line" -ForegroundColor DarkGray
    }
}

function Invoke-Checked {
    param(
        [string]$Command,
        [string[]]$Arguments,
        [string]$Activity
    )
    $commandTimer = [System.Diagnostics.Stopwatch]::StartNew()
    Write-Host "   $Activity"
    $result = Invoke-CapturedCommand -Command $Command -Arguments $Arguments
    if ($result.ExitCode -ne 0) {
        Write-CommandFailure -Command $Command -Arguments $Arguments -Output $result.Output
        throw "$Activity failed (exit code $($result.ExitCode))."
    }
    Write-Host "   Done in $($commandTimer.Elapsed.ToString('hh\:mm\:ss'))." -ForegroundColor Green
}

function Invoke-CondaWithOfflineFallback {
    param(
        [string]$CondaCommand,
        [string[]]$Arguments,
        [string]$Activity
    )
    $commandTimer = [System.Diagnostics.Stopwatch]::StartNew()
    Write-Host "   $Activity"
    $onlineResult = Invoke-CapturedCommand -Command $CondaCommand -Arguments $Arguments
    if ($onlineResult.ExitCode -ne 0) {
        $offlineArguments = @($Arguments) + "--offline"
        Write-Host "   Internet download did not work. Trying the local package cache..." -ForegroundColor Yellow
        $offlineResult = Invoke-CapturedCommand `
            -Command $CondaCommand -Arguments $offlineArguments
        if ($offlineResult.ExitCode -ne 0) {
            Write-CommandFailure -Command $CondaCommand -Arguments $Arguments `
                -Output $onlineResult.Output
            Write-CommandFailure -Command $CondaCommand -Arguments $offlineArguments `
                -Output $offlineResult.Output
            throw "Conda could not complete this step online or from its local cache."
        }
    }
    Write-Host "   Done in $($commandTimer.Elapsed.ToString('hh\:mm\:ss'))." -ForegroundColor Green
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

function Update-AutoQYCheckout {
    param(
        [string]$GitCommand,
        [string]$ProjectRoot,
        [string]$Branch
    )
    $changes = & $GitCommand -C $ProjectRoot status --porcelain
    if ($LASTEXITCODE -ne 0) { throw "Git could not inspect the existing AutoQY checkout." }
    if ($changes) {
        throw "The existing AutoQY checkout has local changes. Preserve them and choose another installation folder, or clean the checkout manually: $ProjectRoot"
    }
    Invoke-Checked -Command $GitCommand -Arguments @(
        "-C", $ProjectRoot, "fetch", "origin",
        "$Branch`:refs/remotes/origin/$Branch"
    ) -Activity "Checking GitHub for AutoQY updates..."
    $currentBranch = (& $GitCommand -C $ProjectRoot branch --show-current).Trim()
    if ($LASTEXITCODE -ne 0) { throw "Git could not determine the current branch." }
    if ($currentBranch -ne $Branch) {
        $localBranch = & $GitCommand -C $ProjectRoot branch --list $Branch
        if ($localBranch) {
            Invoke-Checked -Command $GitCommand -Arguments @(
                "-C", $ProjectRoot, "switch", $Branch
            ) -Activity "Selecting AutoQY version '$Branch'..."
        }
        else {
            Invoke-Checked -Command $GitCommand -Arguments @(
                "-C", $ProjectRoot, "switch", "--track", "-c", $Branch,
                "origin/$Branch"
            ) -Activity "Preparing AutoQY version '$Branch'..."
        }
    }
    Invoke-Checked -Command $GitCommand -Arguments @(
        "-C", $ProjectRoot, "merge", "--ff-only", "origin/$Branch"
    ) -Activity "Applying the downloaded update..."
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
        [string]$CondaCommand,
        [string]$EnvironmentName
    )
    $desktop = [Environment]::GetFolderPath("Desktop")
    $shortcutDirectory = Join-Path $desktop "AutoQY"
    New-Item -ItemType Directory -Path $shortcutDirectory -Force | Out-Null

    foreach ($legacyShortcutName in @(
        "AutoQY Terminal.lnk",
        "AutoQY Analyze JSON.lnk",
        "AutoQY Spectral Smoother.lnk"
    )) {
        $legacyShortcut = Join-Path $shortcutDirectory $legacyShortcutName
        if (Test-Path -LiteralPath $legacyShortcut) {
            Remove-Item -LiteralPath $legacyShortcut -Force
        }
    }

    $shell = New-Object -ComObject WScript.Shell
    $coreCommand = Join-Path $EnvironmentPath "Scripts\autoqy-core.exe"
    $iconDirectory = Join-Path $ProjectRoot "autoqy_core\assets\icons"
    $guiIcon = Join-Path $iconDirectory "power-gui.ico"
    $smootherIcon = Join-Path $iconDirectory "spectral-smoother.ico"
    $terminalIcon = Join-Path $iconDirectory "terminal.ico"
    $analysisIcon = Join-Path $iconDirectory "analyze-json.ico"
    foreach ($iconPath in @($guiIcon, $smootherIcon, $terminalIcon, $analysisIcon)) {
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

    $smootherShortcutPath = Join-Path $shortcutDirectory "AutoQY Spectral Treatment.lnk"
    $smootherShortcut = $shell.CreateShortcut($smootherShortcutPath)
    $smootherShortcut.TargetPath = $coreCommand
    $smootherShortcut.Arguments = "smoother-gui"
    $smootherShortcut.WorkingDirectory = $ProjectRoot
    $smootherShortcut.Description = "Open the AutoQY spectral treatment GUI"
    $smootherShortcut.IconLocation = "$smootherIcon,0"
    $smootherShortcut.Save()

    $analysisShortcutPath = Join-Path $shortcutDirectory "AutoQY Analysis GUI.lnk"
    $analysisShortcut = $shell.CreateShortcut($analysisShortcutPath)
    $analysisShortcut.TargetPath = $coreCommand
    $analysisShortcut.Arguments = "analysis-gui"
    $analysisShortcut.WorkingDirectory = $ProjectRoot
    $analysisShortcut.Description = "Build, save, run, and inspect AutoQY analyses"
    $analysisShortcut.IconLocation = "$analysisIcon,0"
    $analysisShortcut.Save()

    $uninstallScript = Join-Path $ProjectRoot "Uninstall-AutoQY.ps1"
    if (-not (Test-Path -LiteralPath $uninstallScript)) {
        throw "Uninstaller not found: $uninstallScript"
    }

    $uninstallShortcutPath = Join-Path $shortcutDirectory "Uninstall AutoQY.lnk"
    $uninstallShortcut = $shell.CreateShortcut($uninstallShortcutPath)
    $uninstallShortcut.TargetPath = "powershell.exe"
    $uninstallShortcut.Arguments = (
        "-NoLogo -NoProfile -ExecutionPolicy Bypass " +
        "-File `"$uninstallScript`" " +
        "-ProjectRoot `"$ProjectRoot`" " +
        "-EnvironmentName `"$EnvironmentName`" " +
        "-CondaCommand `"$CondaCommand`" " +
        "-ShortcutDirectory `"$shortcutDirectory`""
    )
    $uninstallShortcut.WorkingDirectory = $desktop
    $uninstallShortcut.Description = "Remove AutoQY and optionally its Conda environment"
    $uninstallShortcut.IconLocation = "$terminalIcon,0"
    $uninstallShortcut.Save()

    return @($analysisShortcutPath, $smootherShortcutPath, $guiShortcutPath,
             $uninstallShortcutPath)
}

try {
    Write-Host "AutoQY installer" -ForegroundColor Blue
    Write-Host "This window may stay quiet while a step is working. Please leave it open."

    $InstallDirectory = if ($CheckOnly) {
        [System.IO.Path]::GetFullPath($CurrentDirectory)
    }
    else {
        Select-InstallDirectory -CurrentPath $CurrentDirectory
    }
    $ClonePath = Join-Path $InstallDirectory "AutoQY-Core"

    Write-Host "AutoQY will be installed in: $InstallDirectory"

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
        Write-Host "Environment: $(if ($environmentPath) { "reuse $environmentPath by default; offer clean recreation" } else { "create $EnvironmentName; retry from local cache if online access fails" })"
        Write-Host "Source: $(if ($projectInInstallFolder) { "use existing checkout $projectRoot" } else { "clone branch '$Branch' from $RepositoryUrl into $projectRoot" })"
        Write-Host "Package: editable install with all three browser GUIs"
        Write-Host "Desktop folder AutoQY: Analysis GUI, Spectral Treatment, Power GUI, and uninstaller"
        Write-Host "No files or environments were changed." -ForegroundColor Green
        exit 0
    }

    $reuseEnvironment = $false
    if ($environmentPath) {
        Write-Step "Existing Conda environment detected"
        Write-Host "   $environmentPath" -ForegroundColor Yellow
        $environmentPython = Join-Path $environmentPath "python.exe"
        if (Test-Path -LiteralPath $environmentPython) {
            $reuseEnvironment = Read-Confirmation "Reuse this environment and update AutoQY?" $true
        }
        else {
            Write-Host "   The registered environment is incomplete because python.exe is missing." -ForegroundColor Yellow
        }
        if (-not $reuseEnvironment) {
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
                "env", "remove", "--name", $EnvironmentName, "--yes", "--offline"
            ) -Activity "Removing the old AutoQY environment..."
            $environmentPath = Get-CondaEnvironmentPath -CondaCommand $condaCommand -Name $EnvironmentName
            if ($environmentPath) {
                throw "Conda still reports the environment after removal: $environmentPath"
            }
        }
    }

    if (-not $reuseEnvironment) {
        Write-Step "Preparing a clean Python environment"
        Invoke-CondaWithOfflineFallback -CondaCommand $condaCommand -Arguments @(
            "create", "--name", $EnvironmentName, "python=3.12", "pip", "--yes"
        ) -Activity "Installing Python. This can take several minutes..."
    }
    $environmentPath = Get-CondaEnvironmentPath -CondaCommand $condaCommand -Name $EnvironmentName
    if (-not $environmentPath) { throw "The AutoQY Conda environment could not be located." }
    $environmentPython = Join-Path $environmentPath "python.exe"
    if (-not (Test-Path -LiteralPath $environmentPython)) {
        throw "python.exe was not found in the AutoQY environment."
    }
    $gitCommand = Get-EnvironmentGit -EnvironmentPath $environmentPath
    if (-not $gitCommand) {
        Invoke-CondaWithOfflineFallback -CondaCommand $condaCommand -Arguments @(
            "install", "--name", $EnvironmentName, "git", "--yes"
        ) -Activity "Adding the tools needed to download AutoQY..."
        $gitCommand = Get-EnvironmentGit -EnvironmentPath $environmentPath
        if (-not $gitCommand) { throw "Git was installed but git.exe could not be located in the environment." }
    }

    Activate-AutoQYEnvironment -CondaCommand $condaCommand `
        -Name $EnvironmentName -ExpectedPath $environmentPath

    if (-not $projectInInstallFolder) {
        if (Test-Path -LiteralPath $ClonePath) {
            if (-not (Test-AutoQYProject -Path $ClonePath)) {
                throw "The clone destination already exists but is not an AutoQY project: $ClonePath"
            }
            if (-not (Read-Confirmation "An AutoQY checkout already exists. Update its clean working tree to '$Branch'?" $true)) {
                Write-Host "Installation cancelled."
                exit 0
            }
            Write-Step "Updating the existing AutoQY files"
            Update-AutoQYCheckout -GitCommand $gitCommand -ProjectRoot $ClonePath -Branch $Branch
        }
        else {
            Write-Step "Downloading AutoQY"
            Invoke-Checked -Command $gitCommand -Arguments @(
                "clone", "--branch", $Branch, "--single-branch", $RepositoryUrl, $ClonePath
            ) -Activity "Downloading AutoQY from GitHub..."
        }
    }

    $projectRoot = if ($projectInInstallFolder) { $InstallDirectory } else { $ClonePath }
    Write-Step "Installing AutoQY"
    Push-Location $projectRoot
    try {
        Invoke-Checked -Command $environmentPython -Arguments @(
            "-m", "pip", "install", "--disable-pip-version-check", "--progress-bar", "off",
            "--editable", ".[gui]"
        ) -Activity "Installing AutoQY and its three GUIs. This can take several minutes..."
    }
    finally {
        Pop-Location
    }

    Write-Step "Checking the installation"
    $analysisConfig = Join-Path $projectRoot "ExampleData\Example-2_AB_455nm-100mA\generic_inputs\analysis.json"
    Invoke-Checked -Command $environmentPython -Arguments @(
        "-m", "autoqy_core", "validate", $analysisConfig
    ) -Activity "Checking that AutoQY works..."

    Write-Step "Creating Desktop shortcuts"
    $shortcutTimer = [System.Diagnostics.Stopwatch]::StartNew()
    $null = Write-LauncherFiles -ProjectRoot $projectRoot `
        -EnvironmentPath $environmentPath -CondaCommand $condaCommand `
        -EnvironmentName $EnvironmentName
    Write-Host "   Done in $($shortcutTimer.Elapsed.ToString('hh\:mm\:ss'))." -ForegroundColor Green

    Write-Step "Installation complete"
    Write-Host "AutoQY is ready." -ForegroundColor Green
    Write-Host "Open the AutoQY folder on your Desktop to start."
    Write-Host "Installed in: $projectRoot"
    Write-Host "Total time: $(Get-ElapsedText)"
}
catch {
    Write-Host ""
    Write-Host "Installation failed: $($_.Exception.Message)" -ForegroundColor Red
    exit 1
}

