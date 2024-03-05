# community
`community` is an R package designed to explore the differences in communication between various case and control samples using single-cell RNA sequencing (scRNAseq). With its user-friendly output, `community` lets you easily follow the overall shifts in communication within the cohorts. You can visualize and delve into the most significant differences in interactions, and even investigate what's driving these changes. It's a handy tool for anyone interested in a deeper understanding of cell-to-cell communication.

# What's Inside This Repo
In this repository, you'll find:

- Guidelines on how to leverage the `community` tool effectively
- Methods to visualize output results on preprocessed data
- Insights into the construction of our database
- Instructions for preprocessing raw data across multiple datasets

# The `community` tool is a powerful tool with several benefits:

1. **Differential Communication Analysis:** One of the key capabilities of the "community" is its ability to perform differential communication analysis between cohorts of case and control samples using scRNAseq data. This feature enables the identification of communication differences that may underlie various biological conditions or disease states.

2.  **Compensatory Mechanism Analysis:** Another important aspect offered by the tool is the ability to analyze compensatory mechanisms. Compensatory mechanisms refer to the ways in which cells or organisms compensate for changes or disruptions in a biological process. 

3.  **Speed and Efficiency:** `community` is built for speed, making it suitable for large-scale scRNAseq datasets without high computational demands.

4.  **Robustness:** Our tool employs rigorous statistical methods along with multi-factor analysis to ensure the reliability of the results. This makes it less sensitive to outlier samples, providing more accurate differential communication analysis.

5.  **Intuitive Visualization:** The tool offers straightforward graphical outputs, simplifying the understanding of complex cell interactions.

6.  **User-Friendly:** With an easy-to-use workflow and detailed documentation, the tool is designed for users of varying technical backgrounds.


# Installation
If you are not fimilar with the below code, please follow the next steps to install. We have created a script for Windows users to automatically install all the necesary components (`install_windows.bat`). 

```r
# install.packages("devtools")
devtools::install_github("SoloveyMaria/community")
```
## Step by step installation

To create an environment named 'community_tutorial', and install Jupyter Notebook with the R kernel as well as devtools and all the other dependencies, follow the steps below. All the commands below run in terminal/cmd. 

### For Windows Users:
If you have ```git``` installed on your system, you can directly clone the repository by running the following command in your command prompt or Git Bash:

    
    git clone https://github.com/SoloveyMaria/community.git
    
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

# Custom Database
In `community`, database is managed through a comma-separated file, requiring three mandatory columns(`Pair.Name`, `Ligand` and `Receptor`). Users can easily incorporate their own pairs by editing the provided [.csv file.](https://github.com/SoloveyMaria/community/blob/main/extdata/LR_database.csv) For instructions on updating the database or customization using our workflow, detailed information is available [here.](https://github.com/SoloveyMaria/community/tree/main/docs/DB_building)

# Getting preprocessed data

You can download the preprocessed data by running the below command. You can also visit the link here [https://zenodo.org/record/10619771](https://zenodo.org/records/10619771) and download manually. In order to run the notebook workflow with the preprocessed data, the files should be under the /docs/showcase_notebooks directory. 

- Download preprocessed data by running the following code on the terminal. Windows users can double click on `download_data.bat` to download all the preprocessed data to run the showcase notebooks. The downloaded data will be located at `docs/showcase_notebooks/$dataset/input_files/` 
        - to download Lasry dataset directly into the corresponding directory, you can run `make download-lasry` or you can simply visit [the zenodo page here](https://zenodo.org/records/10619771) to download manually and place them into the related directory.
        - For simillie, `make download-simillie` or visit [zenodo link](https://zenodo.org/records/10512663). 
        - For integrated VanGalen-Oetjen, `make download-vangalen_oetjen` or visit [zenodo link](https://zenodo.org/records/10013368).


    
### Step-by-Step Preprocessing

For each dataset, we provide tailored preprocessing workflows. These workflows encompass:

- Initial Data Cleaning and Annotation
- Filtering
- Normalization
- Batch Correction (if applicable)
- Data Visualization

For a comprehensive overview of our preprocessing steps, please follow this [link](https://github.com/colomemaria/community-paper/tree/main/src/data_preprocessing).
    
# Documentation

To learn using `community`, read one of the following vignettes:

Following vignette contains the explanation on how to perform a basic communication, QC, and differential communication analysis. This includes guidelines on data preprocessing, running the main analysis functions, and interpreting the results. This analysis takes only a few minutes to run: 
- [community analysis on Lasry dataset: ](https://github.com/SoloveyMaria/community/blob/main/docs/showcase_notebooks/Lasry/calculate_communication.ipynb)

If you want to make a comprehensive plots visualization of the community output, you can check following vignettes:

- [community visualization on Lasry dataset: ](https://github.com/SoloveyMaria/community/blob/main/docs/showcase_notebooks/Lasry/visualization.ipynb)

# FAQ:

#### How to Choose the Right Seurat Data Slot for Input Count Matrix?


**Q:** I have X number of samples in a Seurat object. Which data slot should I use as the input count matrix?

**A:** If you're working with a Seurat object, you have multiple options for your input count matrix:
- **Normalized Counts:** The preferred choice, especially for setting meaningful expression levels (`threshold_expr`).
- **Log-Normalized Counts:** A viable alternative, but requires extra caution with `threshold_expr`.
- **Integrated Counts:** Can be used but may introduce artifacts, particularly with lowly expressed genes. Ideal if you're dealing with significant batch effects.

#### Can I Use Community Analysis for Samples with Low Cell Type Frequencies?

**Q:** My scRNA dataset has limited cell type diversity. Is it suitable for Community analysis?

**A:** Absolutely, `community` analysis is not solely dependent on the size of each cell type population. Here's what you should consider:

-    **Presence Across Samples:** A cell type should ideally be present in almost all samples. Missing in one or two samples is generally acceptable.
-    **Minimum Cell Count:** Decide on a threshold for the minimum number of cells of a particular type in a sample. The default in Community is six cells per type per sample.

For more details on filtering based on these criteria that we used, refer to the [Data Pre-Processing Notebook](https://github.com/colomemaria/community-paper/blob/main/src/data_preprocessing/Lasry/2.filtering.ipynb).


#### Is Annotation of Healthy vs. Malignant Cells Necessary?

**Q:** Do I need to differentiate between healthy and malignant B cells when annotating my samples?

**A:** No special annotation is required for malignant B cells. You can categorize them under the general "B cell" type. The algorithm will distinguish between the two populations through the "active fraction" component, effectively capturing the nuances.


## Contribution to community

We welcome contributions from the community! If you're interested in contributing to the project, there are a couple of ways you can get involved:

- Fork and Merge: Feel free to fork the repository, make your changes, and then submit a pull request to merge your changes back into the main branch.


- Open an Issue: If you have ideas for improvements, feature requests, or encounter any bugs, please open an issue on the repository. This helps us keep track of community feedback and prioritize development efforts.


## Cite the tool

Citation will be added soon.
