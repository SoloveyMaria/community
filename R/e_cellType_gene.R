#' @title e_cellType_gene
#' 
#' @description The function calculates the mean expression of each gene in all cell types for a given sample and threshold of expression. It also checks if there are any missing cell types in the annotation of interactions and gives a warning if there are.
#' 
#' @param counts numeric dataframe [0,Inf]: normalized expression data frame containing ligands and receptors in the rows and cells in the columns. Row names and column names should be defined.
#' @param anno_cells dataframe: concatenated data frame of the cell annotation from all samples (rows are cell IDs, columns annotation columns: must contain "cell_ID", "cell_type" and "sample_ID" columns).
#' @param anno_interactions dataframe containing the annotation of interactions.
#' @param sample: string indicating the sample ID.
#' @param threshold_expr numeric: expression threshold below which the expression of the ligand (receptor) is not considered. The default value is 0.05.
#' 
#' @return dataframe (numeric [0,Inf]) containing the mean expression of genes in cell types of the sample of interest.
#' 
#' @export
#' @examples
#' # create a toy example
#' counts <- data.frame(
#'       gene1 = c(0,1,0,1,0,1)
#'    , gene2 = c(1,0,1,0,1,0)
#'   , gene3 = c(0,1,0,1,0,1)
#' )
#' anno_cells <- data.frame(
#'      cell_ID = c('cell1','cell2','cell3','cell4','cell5','cell6')
#'  , cell_type = c('cell_type1','cell_type1','cell_type1','cell_type2','cell_type2','cell_type2')
#' , sample_ID = c('sample1','sample1','sample1','sample1','sample1','sample1')
#' )
#' anno_interactions <- data.frame(
#'     ligand_gene = c('gene1','gene2')
#'  , receptor_gene = c('gene2','gene3')
#' , cell_type_s = c('cell_type1','cell_type1')
#' , cell_type_r = c('cell_type2','cell_type2')
#' , sample_ID = c('sample1','sample1')
#' )
#' # calculate the mean expression of each gene in all cell types for a given sample and threshold of expression
#' e_cellType_gene(counts = counts
#'               ,anno_cells = anno_cells
#'              ,anno_interactions = anno_interactions
#'             ,sample = 'sample1'
#'           ,threshold_expr = 0.5
#' )
#' 
e_cellType_gene <- function(counts, anno_cells, anno_interactions, sample, threshold_expr) {
        means <- meanAboveThreshold(counts = counts, anno_cells = anno_cells, threshold_expr = threshold_expr)
        # check if any cell types are missing
        missing_cell_types <- unique(anno_interactions$sending_cell_type[!anno_interactions$sending_cell_type %in%
                                                                                 colnames(means)])
        if (length(missing_cell_types) != 0) {
                means <- cbind(means, {
                        missing_means <- as.data.frame(matrix(0, ncol = length(missing_cell_types),
                                                              nrow = nrow(means)))
                        colnames(missing_means) <- missing_cell_types
                        missing_means
                })
                warning(paste("WARNING: sample", sample, "does not contain cell type", unique(missing_cell_types),
                              "-- interactions for this cell type in this sample will get zero values."))
        }
        return(means)
}