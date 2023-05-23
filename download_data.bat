@echo off
powershell -Command "Invoke-WebRequest -Uri https://zenodo.org/record/7962808/files/anno_cells_corr.txt -OutFile docs\showcase_notebooks\Lasry\input_data\anno_cells_corr.txt"
powershell -Command "Invoke-WebRequest -Uri https://zenodo.org/record/7962808/files/anno_samples_corr.txt -OutFile docs\showcase_notebooks\Lasry\input_data\anno_samples_corr.txt"
powershell -Command "Invoke-WebRequest -Uri https://zenodo.org/record/7962808/files/counts_corr.csv.gz -OutFile docs\showcase_notebooks\Lasry\input_data\counts_corr.csv.gz"
