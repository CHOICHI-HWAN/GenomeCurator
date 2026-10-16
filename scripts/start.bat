@echo off
rem GenomeCurator launcher for Windows.
rem Prefers the conda/mamba environment "genomecurator" (created from environment.yml on first run);
rem without conda/mamba it falls back to a plain .venv. Extra arguments go to `serve` (e.g. --port 8001).
cd /d "%~dp0.."
set "TOOL="
where mamba >nul 2>nul && set "TOOL=mamba"
if not defined TOOL where conda >nul 2>nul && set "TOOL=conda"
if defined TOOL (
  set "ENV_PY="
  for /f "delims=" %%i in ('%TOOL% run -n genomecurator python -c "import sys; print(sys.executable)" 2^>nul') do set "ENV_PY=%%i"
  if not defined ENV_PY (
    echo Creating conda environment 'genomecurator' from environment.yml ...
    call %TOOL% env create -f environment.yml || (pause & exit /b 1)
    for /f "delims=" %%i in ('%TOOL% run -n genomecurator python -c "import sys; print(sys.executable)"') do set "ENV_PY=%%i"
  )
  "%ENV_PY%" -m genomecurator serve --open %*
  goto :eof
)
where py >nul 2>nul && (set "PY=py -3") || (set "PY=python")
if not exist ".venv\Scripts\python.exe" (
  echo No conda/mamba found: creating a virtual environment (.venv) ...
  %PY% -m venv .venv || (echo Python 3.10+ is required: https://www.python.org/downloads/ & pause & exit /b 1)
  ".venv\Scripts\python.exe" -m pip install --upgrade pip >nul
  ".venv\Scripts\python.exe" -m pip install -e . || (pause & exit /b 1)
)
".venv\Scripts\python.exe" -m genomecurator serve --open %*
