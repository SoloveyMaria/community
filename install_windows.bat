@echo off

setlocal enabledelayedexpansion

for /F "delims=" %%a in ('echo %username%') do (
    set "name=%%a"
)

for /L %%n in (0,1,255) do (
    for /F "delims=áß" %%a in ("!name:~%%n,1!") do (
        if "!name:~%%n,1!" neq "%%a" (
            echo Your username contains non-ASCII characters. Please switch to a user with a username that contains only ASCII characters such as Latin.
            exit /b
        )
    )
)

rem Check if Conda is installed
where conda >nul 2>&1
if %errorlevel% neq 0 (
    rem If Conda is not installed, download and install it
    echo Conda is not installed. Downloading and installing...
    powershell.exe -Command "Invoke-WebRequest https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe -OutFile miniconda.exe"
    start /wait miniconda.exe /S /D="%UserProfile%\Miniconda3"
    del miniconda.exe
    call "%UserProfile%\Miniconda3\Scripts\activate.bat" "%UserProfile%\Miniconda3"
    set "PATH=%UserProfile%\Miniconda3;%UserProfile%\Miniconda3\Scripts;%UserProfile%\Miniconda3\Library\bin;%PATH%"
    call conda init cmd.exe
)

rem Check if community environment exists
conda env list | findstr /R /C:"\bcommunity\b" >nul 2>&1
if %errorlevel% equ 0 (
    rem If community environment exists, run Jupyter Notebook
    echo community environment exists. Running Jupyter Notebook...
    call conda activate community
    jupyter notebook
    call conda deactivate
) else (
    rem If community environment does not exist, create it from environment.yml
    echo community environment does not exist. Creating...
    call conda activate
    conda install -c conda-forge mamba
    mamba env create -f environment.yml
    call conda deactivate
    call conda activate community
    R -e 'options(repos = c(CRAN = "https://cloud.r-project.org/")); install.packages("BiocManager"); BiocManager::install("metaboliteIDmapping"); q()'
    R -e 'options(repos = c(CRAN = "https://cloud.r-project.org/")); install.packages("prettyunits"); devtools::install_github("saezlab/OmnipathR@v3.7.0"); q()'
    R -e 'options(timeout=200); devtools::install_github("SoloveyMaria/community", upgrade = "always"); q()'
    jupyter notebook
    call conda deactivate
)

pause
