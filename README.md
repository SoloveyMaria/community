# community
Community is an R package for analyzing single-cell RNA sequencing data to infer interactions between cell types. The package includes functions for preprocessing and quality control, as well as functions for inferring interactions and for analyzing the inferred interactions.

# Installation
```{r df-drop-ok, class.source="bg-success"}
# install.packages("devtools")
devtools::install_github("SoloveyMaria/community")
```
# Tutorial

To create an environment named 'community', and install Jupyter Notebook with the R kernel, use the command below. 
If you already have Jupyter Notebook with the R kernel installed, you can skip this step.

```conda create -n community -c r -c conda-forge jupyter r-essentials r-irkernel```

# Functionality?
The Community package includes the following functions:

- test_diff(): This function performs a differential expression analysis on the inferred interactions.
- plot_meanLig_vs_meanRec(): This function plots the mean ligand expression vs mean receptor expression.
- plot_pca(): This function performs a PCA analysis on the input data.
- rescale(): This function rescales the input data.
- rho_cellType(): This function calculates the relative fraction of cell types.
- rho(): This function calculates the rho values.
- w(): This function calculates the weight of the interactions between two nodes based on fraction of cells and expression based probability of interaction

# Documentation

Link to the notebooks

# Contribution

We welcome contributions from other developers. To contribute to community, please fork the repository, make your changes and submit a pull request.

# Contact

If you have any questions or suggestions, please contact me at MariaTheBoss
