@echo off
REM UCMuon GUI launcher for Windows
REM UCLouvain Muography Group | Hamid Basiri <hamid.basiri@uclouvain.be>
REM
REM Double-click this file, or run from a terminal:
REM   run_gui.bat
REM
REM Requires: Python 3.9+ with streamlit installed (run install.ps1 first)

cd /d "%~dp0"

REM Set OMP thread count to number of logical cores (falls back to 4).
REM NUMBER_OF_PROCESSORS is a Windows built-in; wmic is removed on Windows 11 24H2+.
set OMP_NUM_THREADS=%NUMBER_OF_PROCESSORS%
if "%OMP_NUM_THREADS%"=="" set OMP_NUM_THREADS=4
echo OMP_NUM_THREADS=%OMP_NUM_THREADS%

REM The first Python that has streamlit, searched as install.ps1 searches:
REM python on PATH, then C:\Python31x (install.ps1 also uses those when python
REM is not on PATH), then the py launcher.
set "PY="
python -c "import streamlit" >nul 2>&1 && set "PY=python"
for %%P in (C:\Python312\python.exe C:\Python311\python.exe C:\Python310\python.exe) do (
    if not defined PY if exist %%P (
        %%P -c "import streamlit" >nul 2>&1 && set "PY=%%P"
    )
)
if not defined PY (
    py -3 -c "import streamlit" >nul 2>&1 && set "PY=py -3"
)
if not defined PY (
    echo.
    echo ERROR: no Python with streamlit found.
    echo Run install_windows.bat first, or:  pip install streamlit
    echo.
    pause
    exit /b 1
)
echo Python: %PY%

%PY% -m streamlit run gui\ucmuon_gui.py
if errorlevel 1 (
    echo.
    echo ERROR: Could not start the GUI.
    echo Make sure streamlit is installed:  pip install streamlit
    echo.
    pause
)
