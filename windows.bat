@ECHO OFF

REM Settings
SET CONDA_ENV=community_tutorial
SET MINICONDA_URL=https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe

IF "%CONDA%" == "" (
    SET BASE=%USERPROFILE%\miniconda3\condabin\conda.bat
) ELSE (
    SET BASE=%~dp0..\conda.bat
)

SET ACTIVATE=%BASE%\Scripts\activate.bat

REM Targets
:create-env
    SETLOCAL EnableDelayedExpansion
    FOR /F %%i IN ('where conda.exe') DO (
        SET CONDA=%%~dpi
        SET BASE=!CONDA!..\condabin\conda.bat
    )
    IF "%CONDA%" == "" (
        curl -L %MINICONDA_URL% -o miniconda.exe
        START /WAIT "" miniconda.exe /InstallationType=JustMe /RegisterPython=0 /S /D=%USERPROFILE%\miniconda3
        DEL /Q miniconda.exe
        FOR /F %%i IN ('where conda.exe') DO (
            SET CONDA=%%~dpi
            SET BASE=!CONDA!..\condabin\conda.bat
        )
    )
    CALL %BASE% activate %CONDA_ENV%
    IF %ERRORLEVEL% NEQ 0 (
        CALL %BASE% create -n %CONDA_ENV%
        CALL %ACTIVATE% %CONDA_ENV%
        CALL %BASE% install -n %CONDA_ENV% -c conda-forge mamba
        CALL mamba env update -n %CONDA_ENV% -f environment.yml
        CALL R -e "devtools::install_github('SoloveyMaria/community', upgrade = 'always'); q()"
    )
    CALL %ACTIVATE% %CONDA_ENV%
    CALL mamba env update -n %CONDA_ENV% -f environment.yml
    CALL R -e "devtools::install_github('SoloveyMaria/community', upgrade = 'always'); q()"


:download-data
    curl https://zenodo.org/record/7565938/files/anno_cells_corr.txt -o src\anno_cells_corr.txt
    curl https://zenodo.org/record/7565938/files/anno_samples_corr.txt -o src\anno_samples_corr.txt
    curl https://zenodo.org/record/7565938/files/counts_corr.csv.gz -o src\counts_corr.csv.gz

:run-jupyter
    CALL %BASE% activate %CONDA_ENV%
    START jupyter notebook

:help
    ECHO Usage:
    FOR /F "tokens=1,* delims=:?" %%A IN ('FINDSTR /B /R /C:"^[a-zA-Z_-][^=]*:.*?##" "%~dpnx0"') DO (
        SET "target=%%A"
        SET "help=%%B"
        CALL :leftPad target padded 30
        ECHO   %padded%  %help%
    )
    EXIT /B

:leftPad
    SETLOCAL
    SET "str=%~1"
    SET /A len=%~2
    SET "padChar= "
    FOR /L %%i IN (1,1,%len%) DO SET "padChar=!padChar! "
    SET "padded=!padChar!!str!"
    SET "padded=!padded:~-%len%!"
    ENDLOCAL & SET "%~3=%padded%"
    EXIT /B
