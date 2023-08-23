# community
`community` is an R package designed to explore the differences in communication between various case and control samples using single-cell RNA sequencing (scRNAseq). With its user-friendly output, `community` lets you easily follow the overall shifts in communication within the cohorts. You can visualize and delve into the most significant differences in interactions, and even investigate what's driving these changes. It's a handy tool for anyone interested in a deeper understanding of cell-to-cell communication.

The `community` tool is a powerful tool with several benefits:

1.  **Compensatory Mechanism Analysis:** Another important aspect offered by the tool is the ability to analyze compensatory mechanisms. Compensatory mechanisms refer to the ways in which cells or organisms compensate for changes or disruptions in a biological process. By using the tool, you can explore and gain insights into the compensatory mechanisms employed by the components in your dataset.
2.  **LogFC Visualization:** The tool provides the ability to visualize LogFC (log-fold change) for each component. LogFC is a measure used to assess the magnitude of change in gene expression between different conditions or groups. With this feature, you can easily analyze and compare the expression levels of different components in your dataset.
3.  .... more


## Installation
If you are not fimilar with the below code, please follow the next steps to install. We have created a script for Windows users to automatically install all the necesary components (`install_windows.bat`). 

```r
# install.packages("devtools")
devtools::install_github("SoloveyMaria/community")
```
## Step by step installation

To create an environment named 'community_tutorial', and install Jupyter Notebook with the R kernel as well as devtools and all the other dependencies, follow the steps below. All the commands below run in terminal/cmd. 

### For Windows Users:
If you have ```git``` installed on your system, you can directly clone the repository by running the following command in your command prompt or Git Bash:

    `
    git clone https://github.com/SoloveyMaria/community.git
    `
**Note:** You need to have Git installed on your system. If Git is not installed, please download this repository in zip format, decompress it, navigate to the directory/folder. Here is the link to the zip file for this repository. [Click here to download](https://github.com/SoloveyMaria/community/archive/refs/heads/main.zip)


Once you navigated to cloned/downloaded folder of the repository `(C:\Users\UserName\Downloads\community)`, you can use the provided **```install_windows```** script (you might need administrative rights) to automate the setup process. You can run this sctips with a double click. Running this script will download and install Conda if it is not already installed on your system. Afterwards, it will create a new Conda environment and install all the dependencies required for the community tool. Upon completion, Jupyter Notebook will open automatically. If you already have all the necessary dependencies installed, the script will only launch Jupyter Notebook.

Please ensure that you have the necessary permissions to run the ```install_windows``` script if you choose to use it.

### Linux/macOS Users:

- Clone the repo with the below command. This will create a folder named "community" in your current directory.

**Note:** Please download and extract the file in here if you do not have `Git` installed on your system. [Click here to download](https://github.com/SoloveyMaria/community/archive/refs/heads/main.zip)   

    ```
    git clone https://github.com/SoloveyMaria/community.git
    ```

- Navigate to the `community` filder by executing the following command.
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

You can download the preprocessed data by running the below command. You can also visit the link here https://zenodo.org/record/7962808 and download manually. In order to run the notebook workflow with the preprocessed data, the files should be under the /docs/showcase_notebooks directory. 

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

If you have any questions or suggestions, please contact me at maria.solovey@bmc.med.lmu.de
