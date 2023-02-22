A brief description of the workflow for building the Ligand-Receptor database for community tool
========

We utilize two main interaction databases from Omnipath. To retrieve the databases we make use of OmnipathR package which is a client for the OmniPath web service (https://www.omnipathdb.org) and many other resources (Valdeolivas et al., 2019).
1. Interactions from the `ligrec extra` dataset of OmniPath. [RefLink](https://r.omnipathdb.org/reference/import_ligrecextra_interactions.html) (Türei et al., 2021)
2. Curated ligand-receptor interactions of OmniPath. [RefLink](https://r.omnipathdb.org/reference/curated_ligand_receptor_interactions.html) (Türei et al., 2021)

The community tool cannot handle complex molecules, so it's necessary to break them down into their components and classify each one as a ligand 
or receptor. To accomplish this, we rely on the OmniPath intercellular communication role annotation database.

#### Breaking down complexes:

Example: lets assume complex G1_G2_G3 is linked to another complex G4_G5_G6. We break down into components and produce all the possible pairwise combinations.

| c1 | c2 | complex_origin    |
|----|----|-------------------|
| G1 | G2 | G1_G2_G3_G4_G5_G6 |
| G1 | G3 | G1_G2_G3_G4_G5_G6 |
| G1 | G4 | G1_G2_G3_G4_G5_G6 |
| G1 | G5 | G1_G2_G3_G4_G5_G6 |
| G1 | G6 | G1_G2_G3_G4_G5_G6 |
| G2 | G1 | G1_G2_G3_G4_G5_G6 |
| G2 | G3 | G1_G2_G3_G4_G5_G6 |
| .. | .. | G1_G2_G3_G4_G5_G6 |

#### Annotation

The complexes are decomposed into their individual components. The Omnipath Intercell annotation database is imported and used to annotate each component
[RefLink](https://r.omnipathdb.org/reference/import_omnipath_intercell.html) (Türei et al., 2021). 
If at least two databases categorize a component as a ligand or receptor, it is annotated as such. If not, we check other possible categories such as 
extracellular matrix, secreted, and transmembrane. 


#### Detecting interaction pairs. 

We are utilizing all of the post-translational datasets from OmniPath, which is the largest network of its kind, to detect interactions rather than 
make predictions. The creators of the network have acknowledged that it may include a significant number of false positives. 
However, we are combining it with an annotations database to identify interactions. The network consists of 98,165 edges, and 
manual curation of interactions is performed once the entire database is built for the community. [RefLink](https://r.omnipathdb.org/reference/import_post_translational_interactions.html)(Türei et al., 2021)

Next, we filter this extensive network to only include the components of the previously decomposed and annotated complexes, 
resulting in a network that comprises only ligand-receptor interactions from complex molecules.

Finally, we can verify whether any of the pairwise combinations exist within this network.

These steps are done separetely for each dataset and merged together. 

#### Mapping Gene Descriptions

**After identifying the ligand and receptor components in each complex and merging them with the single XXX PAIRS**, we want to annotate them with protein descriptions. To achieve this, we utilize `mygene` which is an R package that provides an easy-to-use interface to access the MyGene.info web service, which provides comprehensive annotation information for gene and protein data (Xin et al., 2016). We use the queryMany function to map gene symbols to protein descriptions from the human genome. We then map the protein descriptions to the dataset by matching them with their corresponding gene symbols in each components of each interaction.

Once we have mapped the gene symbols to protein descriptions and incorporated this information into the dataset, we reorder the columns and rename them to ensure consistency across all the datasets. This results in a clean and organized dataset that includes not only information about the interactions but  also the names of the proteins involved. Additionally, we append all of the column information that originates from Omnipath to the ligand-receptor 
interaction data. This allows users to track and see detailed information such as the sources, references, number of curation efforts, 
and number of resources for each interaction. By including this information, we hope to improve the transparency and reliability of the data, 
as users can easily verify the sources and level of curation for each interaction.

#### References
1. Türei, D., Korcsmáros, T., & Saez-Rodriguez, J. (2016). OmniPath: guidelines and gateway for literature-curated signaling pathway resources. Nature methods, 13(12), 966–967. https://doi.org/10.1038/nmeth.4077
2. Valdeolivas A, Turei D, Gabor A (2019). “OmnipathR: client for the OmniPath web service.” Bioconductor Package.
3. Türei, D., Valdeolivas, A., Gul, L., Palacio-Escat, N., Klein, M., Ivanova, O., Ölbei, M., Gábor, A., Theis, F., Módos, D., Korcsmáros, T., & Saez-Rodriguez, J. (2021). Integrated intra- and intercellular signaling knowledge for multicellular omics analysis. Molecular systems biology, 17(3), e9923. https://doi.org/10.15252/msb.20209923
4. Xin, J., Mark, A., Afrasiabi, C., Tsueng, G., Juchler, M., Gopal, N., Stupp, G. S., Putman, T. E., Ainscough, B. J., Griffith, O. L., Torkamani, A., Whetzel, P. L., Mungall, C. J., Mooney, S. D., Su, A. I., & Wu, C. (2016). High-performance web services for querying gene and variant annotation. Genome biology, 17(1), 91. https://doi.org/10.1186/s13059-016-0953-9
