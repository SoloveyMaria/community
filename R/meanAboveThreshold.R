#' @title meanAboveThreshold
#' 
#' @description calculates the mean expression values above the threshold for each cell type and ligand/receptor.
#' 
#' @param cell_type_name character string for the specific cell type to subset the data for. Default is "all" to return data for all cell types.
#' @param counts numeric dataframe: normalized expression data frame containing ligands and receptors in the rows and cells in the columns. Row names and column names should be defined.
#' @param anno_cells dataframe: concatenated data frame of the cell annotation from all samples (rows are cell IDs, columns annotation columns: must contain "cell_ID", "cell_type" and "sample_ID" columns).
#' @param threshold_expr numeric: expression threshold below which the expression of the ligand (receptor) is not considered. The default value is 0.05.
#' 
#' @return a dataframe with the mean expression values above the threshold for each cell type and ligand/receptor.
#' 
#' @export
#' @examples
#' # create a toy example
#' counts <- data.frame(
#'      gene1 = c(0,1,0,1,0,1)
#'  , gene2 = c(1,0,1,0,1,0)
#' , gene3 = c(0,1,0,1,0,1)
#' )
#' anno_cells <- data.frame(
#'     cell_ID = c("cell1","cell2","cell3","cell4","cell5","cell6")
#' , cell_type = c("cell_type1","cell_type1","cell_type1","cell_type2","cell_type2","cell_type2")
#' , sample_ID = c("sample1","sample1","sample1","sample1","sample1","sample1")
#' )
#' # calculate the mean expression of each gene in all cell types for a given sample and threshold of expression
#' meanAboveThreshold(counts = counts
#'                  ,anno_cells = anno_cells
#'                 ,threshold_expr = 0.05)
#' # calculate the mean expression of each gene in cell type 1 for a given sample and threshold of expression
#' meanAboveThreshold(cell_type_name = "cell_type1"
#'                 ,counts = counts
#'                ,anno_cells = anno_cells
#'               ,threshold_expr = 0.05)
#' 
meanAboveThreshold <- function(cell_type_name = "all"
                               ,counts
                               ,anno_cells
                               ,threshold_expr = 0.05){
        
        
        # check if cell IDs are identical in my_counts (column names) and my_anno_cells (row names)
        if(!identical(rownames(anno_cells)
                      ,colnames(counts))){
                stop({
                        "cell IDs are NOT identical in counts (column names) and anno_cells (row names)"
                })
        }
        
        # define all cell typess in the data
        cell_types <- unique(anno_cells[,"cell_type"])
        
        # create empty dataframe nr_ligands_receptors x nr_cell_types
        # we will populate it the the values of the mean expression values
        df_mean_above_threshold <- matrix(,nrow = nrow(counts)
                                          ,ncol = length(cell_types)
        )
        rownames(df_mean_above_threshold) <- rownames(counts)
        colnames(df_mean_above_threshold) <- cell_types
        
        # populate the df_mean_above_threshold
        for(cell_type in cell_types){
                # identify cells that belong to the cell_type
                idx_cell_type <- anno_cells[,"cell_type"] == cell_type
                
                # subset counts to the cell type of interest
                counts_sub <- counts[,idx_cell_type]
                
                # identify counts that pass the threshold
                above_threshold <- counts_sub > threshold_expr
                
                # set all counts which did not pass the threshold to zero
                counts_sub_thr <- counts_sub * above_threshold
                
                # populate the df_mean_above_threshold
                ifelse(is.null(dim(counts_sub_thr))
                       ,df_mean_above_threshold[,cell_type] <- counts_sub_thr
                       ,{
                               # substitute zeros with NAs
                               counts_sub_thr[counts_sub_thr == 0] <- NA
                               df_mean_above_threshold[,cell_type] <- rowMeans(counts_sub_thr
                                                                               ,na.rm = TRUE)
                       }
                       
                )
                # subsitute NAs with zeroes
                df_mean_above_threshold[is.na(df_mean_above_threshold)] <- 0
                
        }
        
        ifelse(cell_type_name == "all"
               ,return(df_mean_above_threshold)
               ,return(df_mean_above_threshold[,cell_type_name]))
}