@echo off
REM UCMuon installer for WINDOWS (on macOS/Linux use: bash setup.sh)
REM UCLouvain Muography Group | Hamid Basiri <hamid.basiri@uclouvain.be>
REM
REM Double-click this file to install UCMuon. It runs install.ps1 with the
REM right execution policy so no PowerShell knowledge is needed.
REM
REM After installation, double-click run_gui.bat to launch the GUI.

cd /d "%~dp0"
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0install.ps1"
echo.
pause
