The input datasets have undergone preprocessing, including filtration and normalization. For detailed insights into the preprocessing steps, you can review the information provided at our paper repository: [preprocessing](https://github.com/colomemaria/community-paper/tree/main/src/data_preprocessing)


### Brief information on datasets
**Lasry:** 

Single-cell RNA sequencing (scRNAseq) data obtained from bone marrow samples of healthy and AML (acute myeloid leukemia) individuals. The raw data used in this study was sourced from publicly available datasets [Lasry, et al. 2022](https://www.nature.com/articles/s43018-022-00480-0) along with annotation files, were downloaded from [GSE185381](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185381).

**VanGalen_Hourigan:** 

To investigate alterations in cell-to-cell communication in AML, we constructed an integrated dataset containing the single-cell RNA sequencing (scRNAseq) profiles of bone marrow from healthy individuals ( [GSE120221 Oetjen et al., 2018](https://doi.org/10.1172/jci.insight.124928)) and AML patients at diagnosis ([GSE116256 vanGalen et al., 2019](https://doi.org/10.1016/j.cell.2019.01.031),). The raw read counts and annotation files were acquired from [GSE116256](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116256) and [GSE120221](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120221) datasets. Duplicated samples S1, Sk1, S2, Ck, and C2 were removed from the GSE120221 dataset. Additionally, the dataset GSE116256 had sample BM5-34p removed due to the absence of total bone marrow. Samples with less than 50% blasts were also excluded, resulting in 35 samples in the joint dataset.
