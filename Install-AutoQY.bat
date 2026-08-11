@echo off
setlocal
cd /d "%~dp0"

echo AutoQY Conda bootstrap installer
echo.

powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -File "%~dp0Install-AutoQY.ps1" %*
set "AUTOQY_EXIT=%ERRORLEVEL%"

echo.
if not "%AUTOQY_EXIT%"=="0" (
    echo The installer stopped with exit code %AUTOQY_EXIT%.
)
pause
exit /b %AUTOQY_EXIT%
