#' @title a_cellType_gene
#'
#' @description calculates the fraction of cells expressing the ligand (receptor) above the threshold in one sample.
#'
#' @param cell_types character vector: all cell types in the data.
#' @param counts numeric dataframe: normalized expression data frame containing ligands and receptors in the rows and cells in the columns. Row names and column names should be defined.
#' @param anno_cells dataframe: concatenated data frame of the cell annotation from all samples (rows are cell IDs, columns annotation columns: must contain "cell_ID", "cell_type" and "sample_ID" columns).
#' @param threshold_expr numeric: expression threshold below which the expression of the ligand (receptor) is not considered. The default value is 0.05.
#' @param threshold_nr_active_cells numeric: threshold for minimal number of active cells, default is 0.
#'
#' @return 
#' list of 2:
#' - df_frac_above_threshold: numeric [0,1] dataframe: rows are genes, columns are cell types. The values are the fraction of cells expressing the ligand(receptor) above the threshold. If a cell type is missing in a sample, all the values for this cell type are set to zero.
#' - df_nr_above_threshold: numeric [0,Inf] dataframe: rows are genes, columns are cell types. The values are the number of cells expressing the ligand(receptor) above the threshold. If a cell type is missing in a sample, all the values for this cell type are set to zero.
#'
#' @export
#' @examples
#' # create a toy example
#' cell_types <- c("cell_type1","cell_type2")
#' counts <- data.frame(
#'        gene1 = c(0,1,0,1,0,1)
#'      ,gene2 = c(0,0,1,1,0,0)
#'     ,gene3 = c(0,0,0,0,1,1)
#'    ,gene4 = c(0,0,0,0,0,0)
#'    ,gene5 = c(1,1,1,1,1,1)
#'   ,gene6 = c(0,0,0,0,0,0)
#'  ,gene7 = c(1,1,1,1,1,1)
#' ,gene8 = c(0,0,0,0,0,0)
#' )
#' rownames(counts) <- paste0("gene",1:nrow(counts))
#' colnames(counts) <- paste0("cell",1:ncol(counts))
#' anno_cells <- data.frame(
#'       cell_ID = colnames(counts)
#'    , cell_type = c("cell_type1","cell_type1","cell_type1","cell_type1","cell_type2","cell_type2")
#'   , sample_ID = c("sample1","sample1","sample1","sample1","sample1","sample1")
#' )
#' # run the function
#' a_cellType_gene(cell_types = cell_types
#'                ,counts = counts
#'               ,anno_cells = anno_cells
#'              ,threshold_expr = 0.05
#'            ,threshold_nr_active_cells = 0
#' )
#'
a_cellType_gene <- function(cell_types
                            ,counts
                            ,anno_cells # rawnames should be identical to colnames of counts
                            ,threshold_expr = 0.05
                            ,threshold_nr_active_cells = 0
){

        
        # check if cell IDs are identical in my_counts (column names) and my_anno_cells (row names)
        if(!identical(rownames(anno_cells)
                      ,colnames(counts))){
                stop({
                        "cell IDs are NOT identical in counts (column names) and anno_cells (row names)"
                })
        }
        
        # create empty dataframe nr_ligands_receptors x nr_cell_types
        # we will populate it the values of the fraction of cells expressing the ligand(receptor(above threshold))
        df_frac_above_threshold <- matrix(0
                                          ,nrow = nrow(counts)
                                          ,ncol = length(cell_types)
        )
        
        rownames(df_frac_above_threshold) <- rownames(counts)
        colnames(df_frac_above_threshold) <- cell_types
        #print(str(df_frac_above_threshold))
        
        # create empty dataframe nr_ligands_receptors x nr_cell_types
        # we will populate it the values of the number of cells expressing the ligand(receptor(above threshold))
        df_nr_above_threshold <- matrix(0
                                        ,nrow = nrow(counts)
                                        ,ncol = length(cell_types)
        )
        
        rownames(df_nr_above_threshold) <- rownames(counts)
        colnames(df_nr_above_threshold) <- cell_types
        
        # populate the df_frac_above_threshold
        for(cell_type in cell_types){
                #print(cell_type)
                
                # identify cells that belong to the cell_type
                idx_cell_type <- anno_cells[,"cell_type"] == cell_type
                #print("sum(idx_cell_type)")
                #print(sum(idx_cell_type))
                
                if(sum(idx_cell_type) != 0){ # check if the cell type exists in the sample
                        
                        # subset counts to the cell type of interest
                        ifelse(sum(idx_cell_type) == 1
                               ,{
                                       counts_sub <- as.data.frame(counts[,idx_cell_type])
                                       colnames(counts_sub) <- colnames(counts)[idx_cell_type]
                                       rownames(counts_sub) <- rownames(counts)
                               }
                               ,counts_sub <- counts[,idx_cell_type])
                        
                        #print("str(counts_sub)")
                        #print(str(counts_sub))
                        
                        # identify counts that pass the threshold
                        above_expr_threshold <- counts_sub > threshold_expr
                        #print("str(above_threshold)")
                        #print(str(above_threshold))
                        #print(dim(above_threshold))
                        
                        # calculate number of cell in the cell type
                        nr_cells_in_cell_type <- sum(idx_cell_type)
                        #print("nr_cells_in_cell_type")
                        #print(nr_cells_in_cell_type)
                        
                        # populate the df_frac_above_threshold
                        ifelse(is.null(dim(above_expr_threshold))
                               ,nr_above_expr_threshold <- as.numeric(above_expr_threshold)
                               ,nr_above_expr_threshold <- rowSums(above_expr_threshold) 
                        )
                        # check if nr active cells (i.e. nr_above_threshold) is higher than the threshold_nr_active_cells
                        idx_above_nr_active_cells_threshold <- nr_above_expr_threshold > threshold_nr_active_cells
                        # set values that didn't pass the threshold to 0
                        nr_above_expr_threshold[!idx_above_nr_active_cells_threshold] <- 0
                        
                        df_frac_above_threshold[,cell_type] <- nr_above_expr_threshold / nr_cells_in_cell_type
                        df_nr_above_threshold[,cell_type] <- nr_above_expr_threshold
                        
                } else {
                        df_frac_above_threshold[,cell_type] <- 0 # if the cell type is missing, phi is equal to 0
                }
        }
        
        #print("str(df_frac_above_threshold)")
        #print(str(df_frac_above_threshold))
        
        return(list(df_frac_above_threshold = df_frac_above_threshold
                    ,df_nr_above_threshold = df_nr_above_threshold
        )
        )
}