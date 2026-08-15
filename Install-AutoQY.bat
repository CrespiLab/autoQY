@echo off
setlocal
cd /d "%~dp0"

echo AutoQY Conda bootstrap installer
echo.

set "AUTOQY_PS1_TARGET=%~dp0Install-AutoQY.ps1"
set "AUTOQY_PS1_URL=https://raw.githubusercontent.com/CrespiLab/autoQY/main/Install-AutoQY.ps1"

echo Downloading the current installer logic beside this BAT file...
powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -Command "$ErrorActionPreference='Stop'; [Net.ServicePointManager]::SecurityProtocol=[Net.SecurityProtocolType]::Tls12; $ProgressPreference='SilentlyContinue'; $target=$env:AUTOQY_PS1_TARGET; $temporary=$target+'.download'; try { Invoke-WebRequest -Uri $env:AUTOQY_PS1_URL -OutFile $temporary -UseBasicParsing; if ((Get-Item -LiteralPath $temporary).Length -lt 1000) { throw 'The downloaded installer is unexpectedly small.' }; Move-Item -LiteralPath $temporary -Destination $target -Force } catch { Remove-Item -LiteralPath $temporary -Force -ErrorAction SilentlyContinue; Write-Error $_; exit 1 }"
if errorlevel 1 (
    echo.
    echo Could not download Install-AutoQY.ps1 from:
    echo %AUTOQY_PS1_URL%
    echo GitHub may be unavailable, or the installer location may have changed.
    pause
    exit /b 1
)

echo Download complete.
echo.
powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -File "%AUTOQY_PS1_TARGET%" %*
set "AUTOQY_EXIT=%ERRORLEVEL%"

echo.
if not "%AUTOQY_EXIT%"=="0" (
    echo The installer stopped with exit code %AUTOQY_EXIT%.
)
pause
exit /b %AUTOQY_EXIT%
