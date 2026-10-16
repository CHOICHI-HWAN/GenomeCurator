@echo off
rem GenomeCurator launcher for Windows: creates .venv on first run, installs, starts the app.
cd /d "%~dp0.."
where py >nul 2>nul && (set "PY=py -3") || (set "PY=python")
if not exist ".venv\Scripts\python.exe" (
  echo Creating virtual environment (.venv) ...
  %PY% -m venv .venv || (echo Python 3.10+ is required: https://www.python.org/downloads/ & pause & exit /b 1)
  ".venv\Scripts\python.exe" -m pip install --upgrade pip >nul
  ".venv\Scripts\python.exe" -m pip install -e . || (pause & exit /b 1)
)
".venv\Scripts\python.exe" -m genomecurator serve --open %*
