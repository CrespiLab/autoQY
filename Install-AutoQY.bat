@echo off
setlocal
cd /d "%~dp0"

echo AutoQY installer
echo.

set "AUTOQY_BRANCH="
set /p "AUTOQY_BRANCH=Git branch to install [main]: "
if not defined AUTOQY_BRANCH set "AUTOQY_BRANCH=main"

set "AUTOQY_BRANCH_TO_VALIDATE=%AUTOQY_BRANCH%"
powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -Command "$branch=$env:AUTOQY_BRANCH_TO_VALIDATE; if ($branch -notmatch '^[A-Za-z0-9][A-Za-z0-9._/-]*$' -or $branch.Contains('..') -or $branch.Contains('//') -or $branch.Contains('@{') -or $branch.EndsWith('/') -or $branch.EndsWith('.') -or $branch.EndsWith('.lock')) { Write-Error 'Invalid Git branch name.'; exit 1 }"
if errorlevel 1 (
    echo.
    echo The branch name is not valid: %AUTOQY_BRANCH%
    echo Use a branch name such as main, develop, or feature/my-change.
    set /p "AUTOQY_CLOSE=Press Enter to close: "
    exit /b 1
)
echo Selected branch: %AUTOQY_BRANCH%
echo.

set "AUTOQY_PS1_TARGET=%~dp0Install-AutoQY.ps1"
set "AUTOQY_PS1_URL=https://raw.githubusercontent.com/CrespiLab/autoQY/%AUTOQY_BRANCH%/Install-AutoQY.ps1"

echo PowerShell installer URL for branch "%AUTOQY_BRANCH%":
echo %AUTOQY_PS1_URL%
echo.
echo Getting the installer from the selected branch...
powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -Command "$ErrorActionPreference='Stop'; [Net.ServicePointManager]::SecurityProtocol=[Net.SecurityProtocolType]::Tls12; $ProgressPreference='SilentlyContinue'; $target=$env:AUTOQY_PS1_TARGET; $temporary=$target+'.download'; try { Invoke-WebRequest -Uri $env:AUTOQY_PS1_URL -OutFile $temporary -UseBasicParsing; if ((Get-Item -LiteralPath $temporary).Length -lt 1000) { throw 'The downloaded installer is unexpectedly small.' }; Move-Item -LiteralPath $temporary -Destination $target -Force } catch { Remove-Item -LiteralPath $temporary -Force -ErrorAction SilentlyContinue; Write-Error $_; exit 1 }"
if errorlevel 1 (
    echo.
    echo Could not download Install-AutoQY.ps1 from:
    echo %AUTOQY_PS1_URL%
    echo GitHub may be unavailable, the branch name may be incorrect, or the installer location may have changed.
    set /p "AUTOQY_CLOSE=Press Enter to close: "
    exit /b 1
)

echo Ready. Starting installation of branch "%AUTOQY_BRANCH%"...
echo When asked for an installation folder, you can copy and paste its full path.
echo.
powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -File "%AUTOQY_PS1_TARGET%" -Branch "%AUTOQY_BRANCH%" -InstallerSourceUrl "%AUTOQY_PS1_URL%" -NoClosePrompt %*
set "AUTOQY_EXIT=%ERRORLEVEL%"

echo.
if not "%AUTOQY_EXIT%"=="0" (
    echo The installer stopped with exit code %AUTOQY_EXIT%.
)
set /p "AUTOQY_CLOSE=Press Enter to close: "
exit /b %AUTOQY_EXIT%
