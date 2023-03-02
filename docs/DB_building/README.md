A brief description of the workflow for building the Ligand-Receptor database for community tool
========

We utilize two main interaction databases from OmniPath (Türei et al., 2021): the [`ligrec extra` dataset](https://r.omnipathdb.org/reference/import_ligrecextra_interactions.html), and the [`curated ligand-receptor interactions`](https://r.omnipathdb.org/reference/curated_ligand_receptor_interactions.html). To retrieve the databases we make use of OmnipathR package (Valdeolivas et al., 2019) which is a client for the OmniPath web service (https://www.omnipathdb.org).

Cellular interactions can occur between individual genes or proteins, where one acts as a transmitter (ligand, L1) and the other as a receiver (receptor, R1). This binary interaction can be represented as L1_R1. However, in some cases, these components can be composed of multiple genes or proteins. These members of a complex may share certain roles or functions, such as all members being receptors or ligands. Alternatively, the complex may be comprised of a mixture of components, such as one ligand and one receptor or one ligand and two receptors, and so on.

All communication tools have a way of scoring interactions by measuring changes. For example, the communication tool assigns weights by calculating rho, phi, and mean values. However, one of the main reasons we do not include complexes in our database is that if a weight value is assigned to a complex, the user would not know how much each component contributes. In other words, the user would not be able to determine which component was affected by any changes observed. Furthermore, the composition of complexes can vary naturally between different cell types, and if the database does not capture this, none of the components will be present. Conversely, using components separately enables users to check how each component behaves individually.

Therefore, we break down these complexes into their individual components, identify their functions as transmitters (ligand) or receivers (receptor), and detect their possible connections. These steps are discussed below in detail.

#### Step 1: break down complexes
 
 To accomplish this, we rely on the [`OmniPath intercellular communication role annotation database`](https://r.omnipathdb.org/reference/import_omnipath_intercell.html) to indetify the functions of the component and [`Intercellular communication network`](https://r.omnipathdb.org/reference/import_intercell_network.html) to detect their binary connections.

Lets assume a complex G1_G2 _G3 is linked to a complex G4_G5_G6 to make a complex interaction G1_G2_G3_G4_G5_G6. We break down the complexes into their components and make a list of transmitter-receiver  pairs by concatenating all the combinations. ~~produce all the possible binary pairwise combinations~~.


| transmitter | receiver | pair  | complex_origin    |
| :---------- | :------- | :---- | :---------------- |
| G1          | G2       | G1_G2 | G1_G2_G3_G4_G5_G6 |
| G1          | G3       | G1_G3 | G1_G2_G3_G4_G5_G6 |
| G1          | G4       | G1_G4 | G1_G2_G3_G4_G5_G6 |
| G1          | G5       | G1_G5 | G1_G2_G3_G4_G5_G6 |
| G1          | G6       | G1_G6 | G1_G2_G3_G4_G5_G6 |
| G2          | G1       | G2_G1 | G1_G2_G3_G4_G5_G6 |
| G2          | G3       | G2_G3 | G1_G2_G3_G4_G5_G6 |
| ..          | ..       |       | G1_G2_G3_G4_G5_G6 |


#### Step 2: detect pairs

After obtaining the binary interactions of complex molecules, the interactions are checked against a large interaction network, which is the [all post-translational datasets of OmniPath](https://r.omnipathdb.org/reference/import_post_translational_interactions.html), the largest database of its kind, containing [139 interaction databases](https://r.omnipathdb.org/reference/get_interaction_resources.html)  as of March 2023. However, the creators of this network have noted that it may contain a significant number of false positives. Therefore, to mitigate this issue, an additional step is taken to eliminate false positive interactions using an annotation database that categorizes the components into ligand or receptor. Specifically, interactions in which a component is annotated as a receptor but is listed in the "source" column (which is designated for transmitters) of the interaction network, or vice versa for ligands, are discarded. This approach ensures that only valid interactions between ligands and receptors are retained.


#### Step 3: annotate individual components

We annotate each individual component by using the [`OmniPath intercellular communication role annotation database`](https://r.omnipathdb.org/reference/import_omnipath_intercell.html).
If at least two databases categorize a component as a ligand or receptor, it is annotated as such. If not, we check other possible categories such as 
*extracellular matrix*, *secreted*, and *transmembrane*. We categorize *extracellular matrix* and *secreted* as transmitter (ligand) while for transmitters we do a manual annotation through genecards and UniProt.

These steps are done separetely for each datasets, namely, `ligrec extra`, and the `curated ligand-receptor interactions`.

#### Step 4: add gene descriptions
As the last step in our pipeline, we add gene names to each individual gene in our database by using publicly available gene annotation resource, namely, MyGene.info web service, which provides comprehensive annotation information for gene and protein data (Xin et al., 2016). This additional information provides users with a more comprehensive understanding of the biological processes and pathways associated with each gene, enabling them to quickly and efficiently explore the functional implications and potential interactions of the complex molecules and their components in our database. By incorporating gene names alongside gene symbols, we aim to enhance the accessibility and interpretability of our database for researchers across a variety of fields, from computational biology to systems pharmacology.


Once we have mapped the gene symbols to protein descriptions and incorporated this information into the dataset, we reorder the columns and rename them to ensure consistency across all the datasets, `ligrec extra`, and the `curated ligand-receptor interactions`. This results in a clean and organized database of ligand-receptor interactions. Additionally, we append all of the column information that originates from OmniPath to our database. This allows users to track and see detailed information such as the sources, references, number of curation efforts, 
and number of resources for each interaction. By including this information, we hope to improve the transparency and reliability of the data, 
as users can easily verify the sources and level of curation for each interaction.

#### References
1. Türei, D., Korcsmáros, T., & Saez-Rodriguez, J. (2016). OmniPath: guidelines and gateway for literature-curated signaling pathway resources. Nature methods, 13(12), 966–967. https://doi.org/10.1038/nmeth.4077
2. Valdeolivas A, Turei D, Gabor A (2019). “OmnipathR: client for the OmniPath web service.” Bioconductor Package.
3. Türei, D., Valdeolivas, A., Gul, L., Palacio-Escat, N., Klein, M., Ivanova, O., Ölbei, M., Gábor, A., Theis, F., Módos, D., Korcsmáros, T., & Saez-Rodriguez, J. (2021). Integrated intra- and intercellular signaling knowledge for multicellular omics analysis. Molecular systems biology, 17(3), e9923. https://doi.org/10.15252/msb.20209923
4. Xin, J., Mark, A., Afrasiabi, C., Tsueng, G., Juchler, M., Gopal, N., Stupp, G. S., Putman, T. E., Ainscough, B. J., Griffith, O. L., Torkamani, A., Whetzel, P. L., Mungall, C. J., Mooney, S. D., Su, A. I., & Wu, C. (2016). High-performance web services for querying gene and variant annotation. Genome biology, 17(1), 91. https://doi.org/10.1186/s13059-016-0953-9
