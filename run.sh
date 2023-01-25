#!/bin/bash

# Settings
CONDA_ENV=community_tutorial
MINICONDA_URL="to be updated"

CONDA=$(which conda)

if [ -z "$CONDA" ]; then
    CONDA=${HOME}/miniconda3/bin/conda
fi

BASE=$(dirname $(dirname $CONDA))
ACTIVATE=$BASE/bin/activate

install_conda() {
    # install Miniconda
    curl -L $MINICONDA_URL -o miniconda.sh
    bash miniconda.sh -b
}

create_env() {
    if conda env list | grep ${CONDA_ENV}; then
        mamba env update -n ${CONDA_ENV} -f environment.yml;
    else
        conda install -n base -c conda-forge mamba &&
        source ${ACTIVATE} base &&
        mamba env create -f environment.yml &&
        source ${ACTIVATE} ${CONDA_ENV} && R -e 'devtools::install_github("SoloveyMaria/community", upgrade = "always"); q()';
    fi
}

run_jupyter() {
    conda $CONDA_ENV && jupyter notebook
}

help() {
    echo "Available options:"
    echo "install-conda  -  install Miniconda"
    echo "create-env  -  create conda environment"
    echo "run-jupyter  -  run jupyter notebooks"
}

case "$1" in
    "install-conda")
        install_conda
        ;;
    "create-env")
        create_env
        ;;
    "run-jupyter")
        run_jupyter
        ;;
    *)
        help
        ;;
esac

