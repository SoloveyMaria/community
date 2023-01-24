#' @title f_cellType_max
#' 
#' @description calculates the maximum fraction of each cell type in all samples, and stores it in the anno_interactions_allSamples list.
#' 
#' @param anno_cells dataframe: concatenated data frame of the cell annotation from all samples (rows are cell IDs, columns annotation columns: must contain "cell_ID", "cell_type" and "sample_ID" columns).
#' @param anno_interactions_allSamples a list of data frames containing annotation information for interactions in all samples.
#' @param threshold_celltype_size numeric: threshold of minimum number of cells in a cell type in a sample. The default value is 4.
#' @param verbose a logical value indicating whether to print progress messages (default: FALSE).
#' 
#' @return anno_interactions_allSamples: list of dataframes (numeric [0,1]) with additional columns for the maximum cell type abundance of the sending cell type, and of the receiving cell type.
#' 
#' @export
#' @examples
#' # calculate the maximum fraction of each cell type in all samples
#' anno_interactions_allSamples <- f_cellType_max(anno_cells
#'                                           ,anno_interactions_allSamples
#'                                         ,verbose = TRUE
#'                                      ,threshold_celltype_size = 4
#'                                   )
#' 
f_cellType_max <- function(anno_cells, anno_interactions_allSamples, threshold_celltype_size = 4,
                           verbose = FALSE) {
        
        if (verbose) {
                print("calculate maximum (over all samples) fraction for each cell type")
        }
        
        my_samples <- names(anno_interactions_allSamples)
        my_cell_types <- unique(anno_cells$cell_type)
        
        # make a dummy df with fractions = 0, row being cell types and columns
        # being samples
        frac_cell_types <- as.data.frame(matrix(0, nrow = length(my_cell_types), ncol = length(my_samples)))
        colnames(frac_cell_types) <- my_samples
        rownames(frac_cell_types) <- my_cell_types
        
        # populate the dataframe with fraction values
        for (my_sample in my_samples) {
                for (my_cell_type in my_cell_types) {
                        
                        # calculate nr of cells in the sample
                        nr_cells_total <- sum(anno_cells$sample_ID == my_sample)
                        
                        # calculate nr if cells of this cell type in this sample
                        nr_cell_celltype <- sum((anno_cells$sample_ID == my_sample) & (anno_cells$cell_type ==
                                                                                               my_cell_type))
                        
                        # calculate fraction of the cell type in the sample check if cell
                        # type size passes the threshold
                        ifelse(nr_cells_total > threshold_celltype_size, frac_cell_types[my_cell_type,
                                                                                         my_sample] <- nr_cell_celltype/nr_cells_total, frac_cell_types[my_cell_type,
                                                                                                                                                        my_sample] <- 0)
                        
                }
        }
        
        # create a vector with max value for each cell type
        max_frac <- apply(frac_cell_types, 1, max)
        
        names(max_frac) <- rownames(frac_cell_types)
        
        # store the values in the anno_interactions_allSamples
        for (sample in names(anno_interactions_allSamples)) {
                # print(sample)
                
                # subset anno_interactions
                anno_interactions_sub <- anno_interactions_allSamples[[sample]]
                
                # create a zero-filled f_s_max and f_r_max columns
                anno_interactions_sub$f_s_max <- rep(0, nrow(anno_interactions_sub))
                anno_interactions_sub$f_r_max <- rep(0, nrow(anno_interactions_sub))
                
                for (my_cell_type in unique(anno_cells$cell_type)) {
                        
                        # find indices for this cell type in the sending_cell_type column
                        # and the receiving_cell_type column
                        idx_s <- anno_interactions_sub$sending_cell_type == my_cell_type
                        idx_r <- anno_interactions_sub$receiving_cell_type == my_cell_type
                        
                        # fill in the f_max values to the correct places in the f_s and f_r
                        # columns
                        if (sum(idx_s != 0)) {
                                anno_interactions_sub$f_s_max[idx_s] <- max_frac[my_cell_type]
                        }
                        if (sum(idx_r != 0)) {
                                anno_interactions_sub$f_r_max[idx_r] <- max_frac[my_cell_type]
                        }
                }
                
                
                # write the values back to the original anno_interactions_allSamples
                anno_interactions_allSamples[[sample]] <- anno_interactions_sub
        }
        
        return(anno_interactions_allSamples)
}