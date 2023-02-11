@echo off

REM Settings
set CONDA_ENV=community_tutorial
set MINICONDA_URL=https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe
set "CONDA=%UserProfile%\miniconda3\condabin\conda.bat"
set "BASE=%UserProfile%\miniconda3"
set "file=src\canno_cells_corr.txt"
set PATH=%PATH%;%UserProfile%\miniconda3;%UserProfile%\miniconda3\Scripts\

if not exist "%CONDA%" set "CONDA=%UserProfile%\miniconda3\condabin\conda.bat"
if not exist "%CONDA%" set "BASE=%UserProfile%\miniconda3"
echo %CONDA%
echo %BASE%


set ACTIVATE=%BASE%\Scripts\activate.bat
set PATH=%PATH%;%UserProfile%\miniconda3\condabin\

:download_data
REM download preprocessed data
if not exist %file% (
  echo File does not exist, downloading...
  curl -o %file% https://zenodo.org/record/7565938/files/anno_cells_corr.txt
) else (
  echo File already exists, skipping download...
)

:install_conda
REM install Miniconda
if not exist "%CONDA%" (
curl -L %MINICONDA_URL% -o miniconda.exe
start /wait /min "Installing Miniconda3..." miniconda.exe /InstallationType=JustMe /S /D="%BASE%"
)

:create_env
REM create or update conda environment
if exist "%CONDA%" (
echo "Conda found at %CONDA%"
if conda.exe info --envs | findstr /C:%CONDA_ENV% > NUL (
echo "Environment %CONDA_ENV% exists, updating..."
call conda activate base
call mamba env update -n %CONDA_ENV% -f environment.yml
) else (
echo "Environment %CONDA_ENV% not found, creating..."
%CONDA% install -n base -c conda-forge mamba
call %ACTIVATE% base
mamba env create -n %CONDA_ENV% -f environment.yml
)
call %ACTIVATE% %CONDA_ENV%
call Rscript "install.R"
) else (
echo "Conda not found, please install and try again."
)

:run_jupyter
REM run jupyter notebooks
call %ACTIVATE% %CONDA_ENV%
jupyter notebook
