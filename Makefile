# Settings
CONDA_ENV=community_tutorial
SHELL=bash
MINICONDA_URL=https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh


UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)
ifeq ($(UNAME_S),Linux)
    MINICONDA_URL=https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    BASE := ${HOME}/miniconda3/bin/conda
endif
ifeq ($(UNAME_S),Darwin)
    ifeq ($(UNAME_M),x86_64)
        MINICONDA_URL=https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
    endif
    ifeq ($(UNAME_M),arm64)
        MINICONDA_URL=https://repo.anaconda.com/miniconda/Miniconda3-py310_23.3.1-0-MacOSX-arm64.sh
    endif
endif

ifeq ($(OS),Windows_NT)
    CONDA := $(strip $(shell where.exe conda))
else
    CONDA := $(strip $(shell which conda))
endif


BASE := $(shell dirname $(shell dirname ${CONDA}))

ACTIVATE=${BASE}/bin/activate


default: help

install-conda: ## install Miniconda
	echo "installing conda"
	echo base $(BASE)
	echo conda $(CONDA)
	echo activate $(ACTIVATE)
	curl -L $(MINICONDA_URL) -o miniconda.sh
	bash miniconda.sh -b
.PHONY: install-conda

create-env: ## create conda environment
	if ${CONDA} env list | grep ${CONDA_ENV}; then \
	   mamba env update -n ${CONDA_ENV} -f environment.yml; \
	   source ${ACTIVATE} ${CONDA_ENV} && R -e 'options(timeout=200); devtools::install_github("SoloveyMaria/community", upgrade = "always"); q()'; \
	else \
	    ${CONDA} install -n base -c conda-forge mamba && \
	    source ${ACTIVATE} base && \
	    mamba env create -f environment.yml && \
	    source ${ACTIVATE} ${CONDA_ENV} && R -e 'options(timeout=200); devtools::install_github("SoloveyMaria/community", upgrade = "always"); q()'; \
	fi
.PHONY: create-env

download-lasry: ## download preprocessed data
	curl https://zenodo.org/records/10513005/files/anno_cells_norm.txt -o docs/showcase_notebooks/Lasry/input_data/anno_cells_norm.txt;
	curl https://zenodo.org/records/10513005/files/anno_samples_norm.txt -o docs/showcase_notebooks/Lasry/input_data/anno_samples_norm.txt;
	curl https://zenodo.org/records/10513005/files/anno_genes_norm.txt -o docs/showcase_notebooks/Lasry/input_data/anno_genes_norm.txt;
	curl https://zenodo.org/records/10513005/files/counts_norm.csv.gz -o docs/showcase_notebooks/Lasry/input_data/counts_norm.csv.gz
.PHONY: download-lasry

download-vangalen_oetjen: ## download preprocessed data
	curl https://zenodo.org/records/10013368/files/anno_cells_corr.txt -o docs/showcase_notebooks/vanGalen_Oetjen/input_data/anno_cells_corr.txt;
	curl https://zenodo.org/records/10013368/files/anno_samples_corr.txt -o docs/showcase_notebooks/vanGalen_Oetjen/input_data/anno_samples_corr.txt;
	curl https://zenodo.org/records/10013368/files/anno_genes_corr.txt -o docs/showcase_notebooks/vanGalen_Oetjen/input_data/anno_genes_corr.txt;
	curl https://zenodo.org/records/10013368/files/counts_corr.csv.gz -o docs/showcase_notebooks/vanGalen_Oetjen/input_data/counts_corr.csv.gz
.PHONY: vangalen_oetjen

download-simillie: ## download preprocessed data
	curl https://zenodo.org/records/10512663/files/anno_cells_norm.txt -o docs/showcase_notebooks/Simillie/input_data/anno_cells_norm.txt;
	curl https://zenodo.org/records/10512663/files/anno_samples_norm.txt -o docs/showcase_notebooks/Simillie/input_data/anno_samples_norm.txt;
	curl https://zenodo.org/records/10512663/files/anno_genes_norm.txt -o docs/showcase_notebooks/Simillie/input_data/anno_genes_norm.txt;
	curl https://zenodo.org/records/10512663/files/counts_norm.csv.gz -o docs/showcase_notebooks/Simillie/input_data/counts_norm.csv.gz
.PHONY: download-simillie

run-jupyter: ## run jupyter notebooks
	 source ${ACTIVATE} ${CONDA_ENV} && \
		jupyter notebook

help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
.PHONY: help

