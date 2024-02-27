A brief description of the workflow for building the Ligand-Receptor database for community tool
========

### TL;DR
The community tool requires three columns: "Pair.Name", "Ligand", and "Receptor". Users can add their own pairs to these respective columns. The original database in CSV format is also provided in the [extdata directory.](https://github.com/SoloveyMaria/community/tree/main/extdata)

To follow our procedures for building or updating the database, you need to install `OmnipathR` and the `mygene` library. You can install them using the following code snippets: 

```R 
devtools::install_github("saezlab/OmnipathR")
```

```R 
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("mygene")
```

Auto-updating the database

```R 
library(community) # load community package
library(mygene)
library(OmnipathR)


LR_database <- auto_update_db("both") 
```
**Disclaimer:** Our database excludes swapped pair duplicates (e.g., L1_R1 and R1_L1) to maintain statistical integrity. We include only one of these pairs. Our method involves checking for the existence of a pair in the Protein-Protein Interaction (PPI) database and cross-referencing other databases for consensus direction. If these steps fail, we default to alphabetical order and designate the pair as an adhesive molecule under the "True_LR" column.

If you wish to update the database and intervene during the processes where you provide your annotations, you can follow the guidelines outlined in the [Basic.ipynb](./Basic.ipynb) file. For more advanced users, detailed code functions are provided in the [Advanced section](Advanced.ipynb), along with comprehensive explanations of their purposes.

In the following section, we delve into our workflow logic, providing detailed insights.

### Construction of the database

There are various databases for biological interactions due to the diverse nature of molecular interactions such as Ligand-Receptor, Enzyme-Substrate, Protein Protein Interactions (PPI) and cell-cell adhesion molecules. Each of these interactions plays a crucial role in biological processes. However, the challenge arises from the differing data formats, coverage, and quality among these databases. We carefully considered several factors to address these challenges to build a comprehensive and reliable database. We needed a database that could offer a wide range of biological interactions, integrate data from various sources and maintain high standards of data quality and coverage. For example, we cross-referenced the interactions with various sources to remove the inaccurate, ensuring that only interactions supported by multiple sources were included into the database (consistency). Additionally, prioritizing interactions that are particularly relevant to intercellular communication processes, filtering out noise and irrelevant entries. 

We picked [OmniPath](https://omnipathdb.org/) ([Türei et al., 2021](https://www.embopress.org/doi/full/10.15252/msb.20209923)) as our database of choice which integrates a vast array of biological interactions and features from 103 databases, including well known databases such as Signor, Reactome and Inact. This integration includes protein-protein interactions (PPI), signaling pathways, regulatory interactions and other molecular interactions as well as annotations of protein function, structure and expression. This unified resource enables researchers to carry out cross database analysis and data interoperability. 

We utilized a combination of datasets from OmniPath to build our database, specifically focusing on interactions from the “[`ligrec extra`](https://r.omnipathdb.org/reference/import_ligrecextra_interactions.html)” , containing ligand-receptor interactions without literature reference (total of 8234 interactions as of July 2023), [`curated ligand-receptor interactions`](https://r.omnipathdb.org/reference/curated_ligand_receptor_interactions.html) (total of 5057 interactions as of July 2023). Besides these, we also make use of intercell annotation dataset that provides information on the roles of individual proteins whether a protein is a ligand, a receptor or an extracellular matrix (ECM) component (total number of annotations 323572 as of July 2023). These datasets were adopted based on their relevance to intercellular communication and their ability to enhance the accuracy and coverage. We aim to create a robust and reliable resource that accurately captures the complex landscape of cellular interactions. 

In our database construction, we considered both binary and complex interactions within the “ligrec extra” and curated ligand-receptor interactions datasets. Cellular interactions can occur between individual proteins, where one acts as a transmitter (Ligand, L1) and the other as a receiver (Receptor, R1). This binary interaction can be represented as L1_R1. However, in some cases, these components can be composed of multiple proteins, forming a complex molecule for the transmitter or the receiver side of the interaction. These members of a complex may share certain roles or functions, such as all members being ligands or receptors. Alternatively, the complex may be composed of a mixture of components, such as one ligand and one receptor or one ligand and two receptors, and so on. 

All communication tools have a way of scoring interactions. For example, the community tool assigns weights by calculating rho, phi and mean values. However, if a weight value is assigned to a complex, the user would not know how much each component contributes. In other words, the user would not be able to determine which component was affected by any changes observed. Furthermore, the composition of complexes can vary naturally between different cell types, and if the database does not capture this, none of the components will be present. Thus, we do not directly include the interactions where one of the components is designated as a complex. Instead, we break down these into a binary interaction, identify their functions as transmitters (ligand) or receivers (receptor) and detect their possible connections. These steps are discussed below in detail. 

##### Step 1: Break down complexes and detect pairs

Let's assume a ligand, L1, is linked to a receptor complex, R1_R2, this pair can be shown as L1_R1_R2. We break it down by producing all the possible binary pairwise combinations (Table X: ). After obtaining the binary pairs of these complex molecules, the pairs are checked against the PPI interaction network of OmniPath, the largest dataset of its kind, containing 142 interaction resources (as of July 2023).  However, we are aware that it might contain a number of false positives. Not only the pairs that are detected through PPI, but also the binary pairs that are coming from the original datasets (ligrec extra and curated ligand-receptor interactions).


| Ligand | Receptor | Pair  | Complex_origin |
|--------|----------|-------|----------------|
| L1     | R1       | L1_R1 | L1_R1_R2       |
| L1     | R2       | L1_R2 | L1_R1_R2       |
| R1     | L1       | R1_L1 | L1_R1_R2       |
| R1     | R2       | R1_R2 | L1_R1_R2       |
| R2     | L1       | R2_L1 | L1_R1_R2       |
| R2     | R1       | R2_R1 | L1_R1_R2       |

In addition to considering the presence of false positives in the PPI network and the original datasets, there is another factor that needs to be taken into account. Some ligand-receptor pairs (L1_R1) may have the components in swapped order, meaning that the receptor is listed before the ligand (R1_L1). However, in our database, the pairs are consistently indicated as L1_R1 with a specific directionality. 

Furthermore, there may be cases where swapped duplicates exist, where a ligand-receptor pair appears multiple times with the components in both orders. For example, there may be instances where L1_R1 and R1_L1 are both present as separate entries. 


Considering these factors, it is important to handle the swapped order and swapped duplicate cases appropriately to ensure the accuracy and integrity of your database, while maintaining the consistent directionality of the ligand-receptor pairs as (L1_R1).



To address these possibilities, we utilize the Intercellular communication network dataset from OmniPath. This dataset combines information from 44 resources (as of July 2023), resulting in a total of 301,777 annotations about the intercellular communication roles of proteins. We annotate all the genes in our database based on these annotations. If a gene is categorized as a ligand or receptor by at least two resources, we classify it as a true ligand or receptor. In cases where genes such as integrins or cadherins are involved, we annotate them as adhesion molecules.

By tagging the true ligand-receptor interactions (True_LR) in our database, users can easily subset and identify these specific interactions. Additionally, if a true ligand-receptor interaction is found in a swapped position, we correct the order accordingly.

For adhesion molecules, we handle swapped duplicates by only retaining the interaction present in the True_LR set and discarding the swapped pair from the adhesion molecule database. In situations where there are swapped duplicates within the adhesion molecule classification, we assess the consensus directionality among the resources. We retain the interactions that have a consensus agreement in directionality, and if there is no consensus, we prioritize the interaction based on lexicographical order. The final database includes 6941 pairs, 

To provide additional information about each individual gene in our database, we incorporate gene descriptions using the publicly available gene annotation resource, MyGene.info web service. This resource offers comprehensive annotation information for gene and protein data (Xin et al., 2016). By including gene descriptions alongside gene symbols, we aim to enhance the understanding of the biological processes and pathways associated with each gene, facilitating exploration of the functional implications and potential interactions of the complex molecules and their components in our database.

Lastly, we append all the column information from OmniPath to our database. This includes detailed information such as sources, references, the number of curation efforts, and the number of resources for each interaction. By including this information, we aim to improve the transparency and reliability of the data. Users can easily track and verify the sources and level of curation for each interaction, enhancing their confidence in the database.


Overall, through these steps and the inclusion of additional information, we aim to construct a comprehensive and reliable database that accurately captures the complex landscape of biological interactions.


#### References
1. Türei, D., Korcsmáros, T., & Saez-Rodriguez, J. (2016). OmniPath: guidelines and gateway for literature-curated signaling pathway resources. Nature methods, 13(12), 966–967. https://doi.org/10.1038/nmeth.4077
2. Valdeolivas A, Turei D, Gabor A (2019). “OmnipathR: client for the OmniPath web service.” Bioconductor Package.
3. Türei, D., Valdeolivas, A., Gul, L., Palacio-Escat, N., Klein, M., Ivanova, O., Ölbei, M., Gábor, A., Theis, F., Módos, D., Korcsmáros, T., & Saez-Rodriguez, J. (2021). Integrated intra- and intercellular signaling knowledge for multicellular omics analysis. Molecular systems biology, 17(3), e9923. https://doi.org/10.15252/msb.20209923
4. Xin, J., Mark, A., Afrasiabi, C., Tsueng, G., Juchler, M., Gopal, N., Stupp, G. S., Putman, T. E., Ainscough, B. J., Griffith, O. L., Torkamani, A., Whetzel, P. L., Mungall, C. J., Mooney, S. D., Su, A. I., & Wu, C. (2016). High-performance web services for querying gene and variant annotation. Genome biology, 17(1), 91. https://doi.org/10.1186/s13059-016-0953-9
