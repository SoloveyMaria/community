# community
Community is an R package for analyzing single-cell RNA sequencing data to infer interactions between cell types. The package includes functions for preprocessing and quality control, as well as functions for inferring interactions and for analyzing the inferred interactions.

## Installation
```{r df-drop-ok, class.source="bg-success"}
# install.packages("devtools")
devtools::install_github("SoloveyMaria/community")
```
## Tutorial

To create an environment named 'community_tutorial', and install Jupyter Notebook with the R kernel as well as devtools and all the other dependencies, follow the steps below. All the commands are below run in terminal/cmd. 

### Steps to open up terminal for Windows users:
Below is how you can open up the terminal on the current folder, this would be the folder where you want to store all the files for the workflows and run community.

Once you navigate to the desired folder. Press <kbd>CTRL</kbd>+<kbd>D</kbd> then type "cmd" and press <kbd>ENTER</kbd>. 

Welcome to dark side!..

![](https://media.licdn.com/dms/image/C5612AQGoiz-xymrlFg/article-cover_image-shrink_600_2000/0/1564461754342?e=2147483647&v=beta&t=ismxInUlGeyIjYdoPItiCPVUCM_KMcY8801GoLdIVjk)

You can follow the instrusctions by copy pasting the below commands to the terminal.

### Requirements

- Clone the repo with the below command. This will create a folder named "community" in your current directory.

**Note:** You need to have Git installed on your system. If Git is not installed, please download this repository in zip format, decompress it, navigate to the directory, and then open up terminal in this folder. Here is the link to the zip file for this repository. [Click here to download](https://github.com/SoloveyMaria/community/archive/refs/heads/main.zip)

    ```
    git clone https://github.com/SoloveyMaria/community.git
    ``` 

- Navigate to the `community` filder by executing the following commands.
    ```
    cd community/
    ```

- If you don't have conda installed yet, install [conda](https://conda.io/miniconda.html) by running the command below

    ```
    make install-conda
    ```

- Create a conda environment named "community_tutorial" and install all necessary packages by using the following command:

    ```
    make create-env
    ```
- Launch [Jupyter](https://jupyter.org/) to access the notebooks to generate graphs

    ```
    make run-jupyter
    ```

- Go to [http://localhost:8888](http://localhost:8888) (a page should open automatically in your browser) 

    or
    
- Open:
    - [`src/calculate_communication.ipynb` Notebook](http://localhost:8888/notebooks/src/extract_data_from_website.ipynb) to run the demo workflow.
    
### Getting preprocessed data

You can download the preprocessed data to /src folder by running the below command. You can also visit the link here https://zenodo.org/record/7565938#.Y9FHVxzMJhE and download manually. In order to run the notebook workflow with the preprocessed data, the files should be under the /src directory. 

- Download preprocessed data into /src directory

    ```
    make download-data
    ```
    
    
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
