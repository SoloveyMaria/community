#' @title a_ct_g_allSamples
#' 
#' @description calculates the relative active fraction for each gene in each cell type for all samples.
#'
#' @param counts numeric dataframe: normalized expression data frame containing ligands and receptors in the rows and cells in the columns. Row names and column names should be defined.
#' @param anno_cells dataframe: concatenated data frame of the cell annotation from all samples (rows are cell IDs, columns annotation columns: must contain 'cell_ID', 'cell_type' and 'sample_ID' column).
#' @param anno_interactions_allSamples list: list containing the interactions for each sample.
#' @param threshold_expr numeric: expression threshold below which the expression of the ligand (receptor) is not considered, default is 0.05.
#' @param threshold_nr_active_cells numeric: threshold for minimal number of active cells, default is 0.
#' @param verbose logical: whether to print progress messages or not, default is FALSE.
#'
#' @return 
#' list:
#' - anno_interactions_allSamples: list containing the interactions for each sample, with additional columns 'a_s_l', 'nr_s_l_active', 'a_r_r' and 'nr_r_r_active' indicating the active fraction and active number of cells for the ligand and receptor in the sending and receiving cell types respectively.
#'
#' @export
#' @examples
#' # load data
#' data("counts")
#' data("anno_cells")
#' data("anno_interactions_allSamples")
#' 
#' 
a_ct_g_allSamples <- function(counts
                              , anno_cells
                              , anno_interactions_allSamples
                              , threshold_expr = 0.05
                              , threshold_nr_active_cells = 0
                              , verbose = FALSE
                              ) {
        
        # calculate a
        if (verbose) {
                print("calculate active fraction for each gene in each cell type")
        }
        
        for (sample in names(anno_interactions_allSamples)) {
                # print(sample)
                anno_interactions_sub <- anno_interactions_allSamples[[sample]]
                
                
                idx_sample <- anno_cells$sample_ID == sample
                # print(idx_sample)
                
                cell_types <- unique(anno_cells$cell_type)
                # print(cell_types)
                
                lig_rec_names <- unique(c(anno_interactions_sub$ligand_gene_name, anno_interactions_sub$receptor_gene_name))
                # keep only lig-rec present in the count matrix
                lig_rec_names <- lig_rec_names[lig_rec_names %in% rownames(counts)]
                # print(lig_rec_names)
                
                # subset counts and anno_cells
                if (sum(idx_sample) == 1) {
                        
                        counts_sub <- counts[lig_rec_names, idx_sample]
                        counts_sub <- as.data.frame(counts_sub)
                        colnames(counts_sub) <- colnames(counts)[idx_sample]
                        rownames(counts_sub) <- rownames(counts[lig_rec_names, ])
                        
                        anno_cells_sub <- anno_cells[idx_sample, ]
                        
                } else {
                        # print(str(counts)) print(str(idx_sample))
                        counts_sub <- counts[lig_rec_names, idx_sample]
                        anno_cells_sub <- anno_cells[idx_sample, ]
                }
                
                # print(str(counts_sub)) print(str(anno_cells_sub))
                
                # calculate a in each cell type for each gene n this sample
                # print('calculate a_cellType_gene')
                a <- a_cellType_gene(cell_types = cell_types, counts = counts_sub, anno_cells = anno_cells_sub,
                                     threshold_expr = threshold_expr, threshold_nr_active_cells = threshold_nr_active_cells)
                # print('str(a)') print(str(a))
                
                # write the values to the anno_interactions list
                for (cell_type in unique(anno_cells$cell_type)) {
                        
                        # print(cell_type)
                        
                        idx_cell_type <- anno_cells$cell_type == cell_type
                        # print(idx_cell_type)
                        
                        idx_send <- anno_interactions_sub$sending_cell_type == cell_type
                        anno_interactions_sub$a_s_l[idx_send] <- a$df_frac_above_threshold[anno_interactions_sub$ligand_gene_name[idx_send],
                                                                                           cell_type]
                        anno_interactions_sub$nr_s_l_active[idx_send] <- a$df_nr_above_threshold[anno_interactions_sub$ligand_gene_name[idx_send],
                                                                                                 cell_type]
                        
                        idx_rec <- anno_interactions_sub$receiving_cell_type == cell_type
                        anno_interactions_sub$a_r_r[idx_rec] <- a$df_frac_above_threshold[anno_interactions_sub$receptor_gene_name[idx_rec],
                                                                                          cell_type]
                        anno_interactions_sub$nr_r_r_active[idx_rec] <- a$df_nr_above_threshold[anno_interactions_sub$receptor_gene_name[idx_rec],
                                                                                                cell_type]
                        
                }
                
                
                anno_interactions_allSamples[[sample]] <- anno_interactions_sub
        }
        anno_interactions_allSamples
}