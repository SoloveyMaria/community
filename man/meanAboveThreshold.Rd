% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meanAboveThreshold.R
\name{meanAboveThreshold}
\alias{meanAboveThreshold}
\title{meanAboveThreshold}
\usage{
meanAboveThreshold(
  cell_type_name = "all",
  counts,
  anno_cells,
  threshold_expr = 0.05
)
}
\arguments{
\item{cell_type_name}{character string for the specific cell type to subset the data for. Default is "all" to return data for all cell types.}

\item{counts}{numeric dataframe: normalized expression data frame containing ligands and receptors in the rows and cells in the columns. Row names and column names should be defined.}

\item{anno_cells}{dataframe: concatenated data frame of the cell annotation from all samples (rows are cell IDs, columns annotation columns: must contain "cell_ID", "cell_type" and "sample_ID" columns).}

\item{threshold_expr}{numeric: expression threshold below which the expression of the ligand (receptor) is not considered. The default value is 0.05.}
}
\value{
a dataframe with the mean expression values above the threshold for each cell type and ligand/receptor.
}
\description{
calculates the mean expression values above the threshold for each cell type and ligand/receptor.
}
\examples{
# create a toy example
counts <- data.frame(
     gene1 = c(0,1,0,1,0,1)
 , gene2 = c(1,0,1,0,1,0)
, gene3 = c(0,1,0,1,0,1)
)
anno_cells <- data.frame(
    cell_ID = c("cell1","cell2","cell3","cell4","cell5","cell6")
, cell_type = c("cell_type1","cell_type1","cell_type1","cell_type2","cell_type2","cell_type2")
, sample_ID = c("sample1","sample1","sample1","sample1","sample1","sample1")
)
# calculate the mean expression of each gene in all cell types for a given sample and threshold of expression
meanAboveThreshold(counts = counts
                 ,anno_cells = anno_cells
                ,threshold_expr = 0.05)
# calculate the mean expression of each gene in cell type 1 for a given sample and threshold of expression
meanAboveThreshold(cell_type_name = "cell_type1"
                ,counts = counts
               ,anno_cells = anno_cells
              ,threshold_expr = 0.05)

}
