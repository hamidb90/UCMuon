@echo off
REM UCMuon GUI launcher for Windows
REM UCLouvain Muography Group | Hamid Basiri <hamid.basiri@uclouvain.be>
REM
REM Double-click this file, or run from a terminal:
REM   run_gui.bat
REM
REM Requires: Python 3.9+ with streamlit installed (run install.ps1 first)

cd /d "%~dp0"

REM Set OMP thread count to number of logical cores (falls back to 4)
for /f "tokens=2 delims==" %%i in ('wmic cpu get NumberOfLogicalProcessors /value ^| findstr =') do set OMP_NUM_THREADS=%%i
if "%OMP_NUM_THREADS%"=="" set OMP_NUM_THREADS=4
echo OMP_NUM_THREADS=%OMP_NUM_THREADS%

python -m streamlit run gui\ucmuon_gui.py
if errorlevel 1 (
    echo.
    echo ERROR: Could not start the GUI.
    echo Make sure streamlit is installed:  pip install streamlit
    echo.
    pause
)
