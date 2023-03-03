#' @title e_ct_g_allSamples
#' 
#' @description The function calculates the mean expression of each gene in all cell types for all samples and a given threshold of expression. It also checks if there are any missing cell types in the annotation of interactions and gives a warning if there are.
#' 
#' @param counts numeric dataframe [0,Inf]: normalized expression data frame containing ligands and receptors in the rows and cells in the columns. Row names and column names should be defined.
#' @param anno_cells dataframe: concatenated data frame of the cell annotation from all samples (rows are cell IDs, columns annotation columns: must contain "cell_ID", "cell_type" and "sample_ID" columns).
#' @param anno_interactions_allSamples list containing the interactions for each sample.
#' @param threshold_expr numeric: expression threshold below which the expression of the ligand (receptor) is not considered. The default value is 0.05.
#' @param verbose a logical value indicating whether to print progress messages (default: FALSE).
#' 
#' @return anno_interactions_allSamples: list containing the interactions for each sample, with additional columns e_s_l and e_r_r indicating the mean expression in the active cell fraction for the ligand and receptor in the sending and receiving cell types respectively.
#' 
#' @export
#'
e_ct_g_allSamples <- function(counts, anno_cells, anno_interactions_allSamples, threshold_expr = 0.05,
                              verbose = FALSE) {
        # calculate mean expression
        if (verbose) {
                print("calculate mean expression in the active fraction")
        }
        
        # check if cell IDs are identical in my_counts (column names) and my_anno_cells (row names)
        if (!identical(rownames(anno_cells), colnames(counts))) {
                stop({
                        "cell IDs are NOT identical in counts (column names) and anno_cells (row names)"
                })
        }
        
        # calculate the e_cellType_gene values for each sample
        for (sample in names(anno_interactions_allSamples)) {
                
                idx_sample <- anno_cells$sample_ID == sample
                
                # subset anno_interactions
                anno_interactions_sub <- anno_interactions_allSamples[[sample]]
                
                # subset counts
                ifelse(sum(idx_sample) == 1, {
                        counts_sub <- as.data.frame(counts[, idx_sample])
                        colnames(counts_sub) <- colnames(counts)[idx_sample]
                        rownames(counts_sub) <- rownames(counts)
                }, counts_sub <- counts[, idx_sample])
                
                # subset anno_cells
                anno_cells_sub <- anno_cells[idx_sample, ]
                
                means <- e_cellType_gene(counts = counts_sub, anno_cells = anno_cells_sub,
                                         anno_interactions = anno_interactions_sub, sample = sample, threshold_expr = threshold_expr)

                means <- as.matrix(means)
                
                matrix_melt <- reshape2::melt(means, value.name = "value")

                anno_interactions_sub=dplyr::left_join(anno_interactions_sub, matrix_melt, 
                                by = c("ligand_gene_name" = "Var1", "sending_cell_type" = "Var2")) %>% 
                                dplyr::rename(e_s_l = value)

                anno_interactions_sub=dplyr::left_join(anno_interactions_sub, matrix_melt, 
                                by = c("receptor_gene_name" = "Var1", "receiving_cell_type" = "Var2")) %>% 
                                dplyr::rename(e_r_r = value)

                #anno_interactions_sub$e_s_l <- diag(as.matrix(means[anno_interactions_sub$ligand_gene_name,
                #                                                    anno_interactions_sub$sending_cell_type]))
                #anno_interactions_sub$e_r_r <- diag(as.matrix(means[anno_interactions_sub$receptor_gene_name,
                #                                                    anno_interactions_sub$receiving_cell_type]))
                
                # write the values back to the original anno_interactions_allSamples
                anno_interactions_allSamples[[sample]] <- anno_interactions_sub
        }
        return(anno_interactions_allSamples)
}