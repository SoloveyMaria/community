#' @title e_cellType_gene_max
#' 
#' @description calculates max mean expression of each gene in all cell types in all samples.
#'
#' @param anno_interactions_allSamples a list of data frames containing annotation information for interactions in all samples.
#' @param verbose a logical value indicating whether to print progress messages (default: FALSE).
#' 
#' @return a list of data frames (numeric [0,Inf]) with additional columns for the maximum mean expression for ligand in sending cell type, and for the maximum active fraction for receptor in receiving cell type.
#' 
#' @export
e_cellType_gene_max <- function(anno_interactions_allSamples
                                ,verbose = FALSE
){
        
        # calculate max mean expression for ligand in sending cell type
        if(verbose){print("calculate max mean expression for ligand in sending cell type")}
        max_mean_expr_l_s <- apply(
                #rbind(
                sapply(names(anno_interactions_allSamples) # for each sample
                       ,function(my_sample){
                               anno_interactions_allSamples[[my_sample]][,"e_s_l"] # extract the vector of values of expr_l_s_active (for all interactions)
                       }
                )
                #)
                ,1
                ,max
        )
        for(my_sample in names(anno_interactions_allSamples)){
                anno_interactions_allSamples[[my_sample]]$e_s_l_max <- max_mean_expr_l_s
        }
        
        # calculate max mean expression for receptor in receiving cell type
        if(verbose){print("calculate max mean expression for receptor in receiving cell type")}
        max_mean_expr_r_r <- apply(
                #rbind(
                sapply(names(anno_interactions_allSamples) # iterate over all samples
                       ,function(my_sample){
                               anno_interactions_allSamples[[my_sample]][,"e_r_r"]
                       }
                )
                #)
                ,1
                ,max
        )
        for(my_sample in names(anno_interactions_allSamples)){
                anno_interactions_allSamples[[my_sample]]$e_r_r_max <- max_mean_expr_r_r
        }
        
        return(anno_interactions_allSamples)
        
}
