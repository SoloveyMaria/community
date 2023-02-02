# Settings
CONDA_ENV=community_tutorial
SHELL=bash
MINICONDA_URL=https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

ifeq ($(OS),Windows_NT)
    CONDA := $(strip $(shell where.exe conda))
else
    CONDA := $(strip $(shell which conda))
endif

ifeq ($(CONDA),"")
    BASE := ${HOME}/miniconda3/bin/conda
else
    BASE := $(shell dirname $(shell dirname ${CONDA}))
endif

ACTIVATE=${BASE}/bin/activate


default: help

install-conda: ## install Miniconda
	curl -L $(MINICONDA_URL) -o miniconda.sh
	bash miniconda.sh -b
.PHONY: install-conda

create-env: ## create conda environment
	if ${CONDA} env list | grep ${CONDA_ENV}; then \
	   mamba env update -n ${CONDA_ENV} -f environment.yml; \
	   source ${ACTIVATE} ${CONDA_ENV} && R -e 'devtools::install_github("SoloveyMaria/community", upgrade = "always"); q()'; \
	else \
	    conda install -n base -c conda-forge mamba && \
	    source ${ACTIVATE} base && \
	    mamba env create -f environment.yml && \
	    source ${ACTIVATE} ${CONDA_ENV} && R -e 'devtools::install_github("SoloveyMaria/community", upgrade = "always"); q()'; \
	fi
.PHONY: create-env

download-data: ## download preprocessed data
	curl https://zenodo.org/record/7565938/files/anno_cells_corr.txt -o src/anno_cells_corr.txt;
	curl https://zenodo.org/record/7565938/files/anno_samples_corr.txt -o src/anno_samples_corr.txt;
	curl https://zenodo.org/record/7565938/files/counts_corr.csv.gz -o src/counts_corr.csv.gz
.PHONY: download-data

run-jupyter: ## run jupyter notebooks
	 source ${ACTIVATE} ${CONDA_ENV} && \
		jupyter notebook

help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
.PHONY: help
