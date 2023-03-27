@echo off

rem Check if Conda is installed
where conda >nul 2>&1
if %errorlevel% neq 0 (
    rem If Conda is not installed, download and install it
    echo Conda is not installed. Downloading and installing...
    powershell.exe -Command "Invoke-WebRequest https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe -OutFile miniconda.exe"
    start /wait miniconda.exe /S /D="%UserProfile%\Miniconda3"
    set "PATH=%UserProfile%\Miniconda3;%UserProfile%\Miniconda3\Scripts;%UserProfile%\Miniconda3\Library\bin;%PATH%"
)

rem Check if community_tutorial environment exists
conda env list | findstr /C:"community_tutorial" >nul 2>&1
if %errorlevel% equ 0 (
    rem If community_tutorial environment exists, run Jupyter Notebook
    echo community_tutorial environment exists. Running Jupyter Notebook...
    conda activate community_tutorial
    jupyter notebook
    conda deactivate
) else (
    rem If community_tutorial environment does not exist, create it from env.yaml
    echo community_tutorial environment does not exist. Creating...
    conda env create -f environment.yml
    conda activate community_tutorial
    jupyter notebook
    conda deactivate
)

pause
