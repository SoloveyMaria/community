# community
Community is an R package for analyzing single-cell RNA sequencing data to infer interactions between cell types. The package includes functions for preprocessing and quality control, as well as functions for inferring interactions and for analyzing the inferred interactions.

## Installation
```{r df-drop-ok, class.source="bg-success"}
# install.packages("devtools")
devtools::install_github("SoloveyMaria/community")
```
## Tutorial

To create an environment named 'community_tutorial', and install Jupyter Notebook with the R kernel as well as devtools and all the other dependencies, use the command below. 

### Requirements

- Install [conda](https://conda.io/miniconda.html)

    ```
    /.run.sh install-conda
    ```

- Create the conda environment

    ```
    /.run.sh create-env
    ```
- Launch [Jupyter](https://jupyter.org/) to access the notebooks to generate graphs

    ```
    $ /.run.sh run-jupyter
    ```

- Go to [http://localhost:8888](http://localhost:8888) (a page should open automatically in your browser)
- Open:
    - [`src/calculate_comminocaton.ipynb` Notebook](http://localhost:8888/notebooks/src/extract_data_from_website.ipynb) to extract information from the website
## Functionality
The Community package includes the following functions:

- test_diff(): This function performs a differential expression analysis on the inferred interactions.
- plot_meanLig_vs_meanRec(): This function plots the mean ligand expression vs mean receptor expression.
- plot_pca(): This function performs a PCA analysis on the input data.
- rescale(): This function rescales the input data.
- rho_cellType(): This function calculates the relative fraction of cell types.
- rho(): This function calculates the rho values.
- w(): This function calculates the weight of the interactions between two nodes based on fraction of cells and expression based probability of interaction

## Documentation

Link to the notebooks

## Contribution

We welcome contributions from other developers. To contribute to community, please fork the repository, make your changes and submit a pull request.

## Contact

If you have any questions or suggestions, please contact me at MariaTheBoss
