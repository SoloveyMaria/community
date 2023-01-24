#' @title calculate_weight
#' @description calculates individual components and the weight per sample.
#' 
#' @param anno_interactions dataframe: contains annotations to the interactions in the correspoding sample. The dataframe contains interactions IDs as row names, as well as the following columns: interaction_ID, ligand_gene_name, receptor_gene_name, sending_cell_type, receiving_cell_type, f_s, f_s_max, f_r, f_r_max, a_s_l, a_s_l_max, a_r_r, a_r_r_max, e_s_l, e_s_l_max, e_r_r, e_r_r_max.
#' @param verbose a logical value indicating whether to print progress messages (default: FALSE).
#' 
#' @return anno_interactions dataframe with the following columns added: rho_s, rho_r, rho, phi_s_l, phi_r_r, phi, p_s_l, p_r_r, p, w.
#' 
#' @export
#' @examples
#' anno_interactions <- data.frame(
#'        interaction_ID = c('interaction1','interaction2','interaction3')
#'      ,ligand_gene_name = c('gene1','gene2','gene3')
#'     ,receptor_gene_name = c('gene1','gene2','gene3')
#'   ,sending_cell_type = c('cellType1','cellType2','cellType3')
#' ,receiving_cell_type = c('cellType1','cellType2','cellType3')
#' ,mean_expr_l_s = c(0.1,0.2,0.3)
#' ,mean_expr_r_r = c(0.1,0.2,0.3)
#' ,max_mean_expr_all_samples = c(0.1,0.2,0.3)
#' ,phi_l_s = c(0.1,0.2,0.3)
#' ,phi_r_r = c(0.1,0.2,0.3)
#' )
#' anno_interactions <- calculate_weight(anno_interactions = anno_interactions
#'                                       ,verbose = TRUE
#' )
#' 
calculate_weight <- function(anno_interactions, verbose = FALSE) {
        
        if (verbose) {
                print("calculate interactions weights")
        }
        
        ## rho #### calculate rho (normalized cell type abundance) for the sending
        ## cell types
        anno_interactions$rho_s <- rho_cellType(f_cellType = anno_interactions$f_s, f_cellType_max = anno_interactions$f_s_max,
                                                verbose = verbose)
        
        # calculate rho (normalized cell type abundance) for the receiving cell
        # types
        anno_interactions$rho_r <- rho_cellType(f_cellType = anno_interactions$f_r, f_cellType_max = anno_interactions$f_r_max,
                                                verbose = verbose)
        
        # calculate rho for all interactions
        anno_interactions$rho <- rho(rho_s = anno_interactions$rho_s, rho_r = anno_interactions$rho_r,
                                     verbose = verbose)
        
        # calculate phi (normalized active fraction) for the sending cell types
        anno_interactions$phi_s_l <- phi_cellType_gene(a_cellType_gene = anno_interactions$a_s_l,
                                                       a_cellType_gene_max = anno_interactions$a_s_l_max, verbose = verbose)
        
        # calculate phi (normalized active fraction) for the receiving cell types
        anno_interactions$phi_r_r <- phi_cellType_gene(a_cellType_gene = anno_interactions$a_r_r,
                                                       a_cellType_gene_max = anno_interactions$a_r_r_max, verbose = verbose)
        
        
        ## phi #### calculate phi for all interactions
        anno_interactions$phi <- phi(phi_s_l = anno_interactions$phi_s_l, phi_r_r = anno_interactions$phi_r_r,
                                     verbose = verbose)
        
        ## p #### calculate p_l_s (normalized expression strength) for all
        ## interactions
        anno_interactions$p_s_l <- p_cellType_gene(e_cellType_gene = anno_interactions$e_s_l,
                                                   e_cellType_gene_max = anno_interactions$e_s_l_max, verbose = verbose)
        
        # calculate p_r_r (normalized expression strength) for all interactions
        anno_interactions$p_r_r <- p_cellType_gene(e_cellType_gene = anno_interactions$e_r_r,
                                                   e_cellType_gene_max = anno_interactions$e_r_r_max, verbose = verbose)
        
        # calculate p for all interactions
        anno_interactions$p <- p(p_s_l = anno_interactions$p_s_l, p_r_r = anno_interactions$p_r_r)
        
        
        
        ## w #### calculate w for all interactions
        anno_interactions$w <- w(rho = anno_interactions$rho, phi = anno_interactions$phi,
                                 p = anno_interactions$p, verbose = verbose)
        
        return(anno_interactions)
        
}