#' @title f_ct_allSamples 
#' 
#' @description calculates fraction of each cell type in each sample.
#'
#' @param anno_cells dataframe: concatenated data frame of the cell annotation from all samples (rows are cell IDs, columns annotation columns: must contain "cell_ID", "cell_type" and "sample_ID" columns).
#' @param anno_interactions_allSamples a list of data frames containing annotation information for interactions in all samples.
#' @param threshold_celltype_size numeric [0,Inf]: threshold of minimum number of cells in a cell type in a sample. The default value is 4.
#' @param verbose a logical value indicating whether to print progress messages (default: FALSE).
#' 
#' @return a dataframe (numeric [0,1]) with fractions of cell types where columns are samples and rows are cell types.
#' 
#' @export
#' @examples
#' # calculate fraction of each cell type in each sample
#' anno_interactions_allSamples <- f_ct_allSamples(anno_cells
#'                                         ,anno_interactions_allSamples
#'                                        ,threshold_celltype_size = 4
#'                                       ,verbose = TRUE
#' )
#' 
f_ct_allSamples <- function(anno_cells, anno_interactions_allSamples, threshold_celltype_size = 4,
                            verbose = FALSE) {
        
        if (verbose) {
                print("calculate fraction of each cell type in each sample")
        }
        
        
        # calculate fraction for each cell type in each sample
        for (sample in names(anno_interactions_allSamples)) {
                
                # print(sample)
                
                # subset anno_interactions
                anno_interactions_sub <- anno_interactions_allSamples[[sample]]
                
                # subset anno_cells
                anno_cells_sub <- anno_cells[anno_cells$sample_ID == sample, ]
                
                # create a zero-filled f_s and f_r columns
                anno_interactions_sub$f_s <- rep(0, nrow(anno_interactions_sub))
                anno_interactions_sub$f_r <- rep(0, nrow(anno_interactions_sub))
                
                # check fractions for each cell type
                for (my_cell_type in unique(anno_cells$cell_type)) {
                        
                        # print(my_cell_type)
                        
                        # find indices for this cell type in the sending_cell_type column
                        # and the receiving_cell_type column
                        idx_s <- anno_interactions_sub$sending_cell_type == my_cell_type
                        idx_r <- anno_interactions_sub$receiving_cell_type == my_cell_type
                        
                        f <- f_cellType(anno_interactions = anno_interactions_sub, cell_type = my_cell_type,
                                        anno_cells = anno_cells_sub, threshold_celltype_size = threshold_celltype_size)
                        
                        # fill in the f values to the correct places in the f_s and f_r
                        # columns
                        if (sum(idx_s != 0)) {
                                anno_interactions_sub$f_s[idx_s] <- f
                        }
                        if (sum(idx_r != 0)) {
                                anno_interactions_sub$f_r[idx_r] <- f
                        }
                }
                
                # write the values back to the original anno_interactions_allSamples
                anno_interactions_allSamples[[sample]] <- anno_interactions_sub
        }
        
        return(anno_interactions_allSamples)
}