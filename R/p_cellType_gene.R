#' @title p_cellType_gene
#'
#' @description calculates the probability of interaction of a ligand (or receptor) with its counterpart receptor (or ligand) based on its expression level in the cell type of interest.
#'
#' @param e_cellType_gene numeric [0,Ifn]: mean expression level of the ligand (or receptor) in the active fraction of the cell type of interest.
#' @param e_cellType_gene_max numeric [0,Ifn]: maximum expression level of the ligand (or receptor) in the cell type of interest across all samples.
#' @param verbose a logical value indicating whether to print progress messages (default: FALSE).
#'
#' @return p_expr numeric [0,1]: probability of interaction, normalized to the maximum expression level across all cell types.
#'
#' @export
#' @examples
#' e_cellType_gene <- 0.1
#' e_cellType_gene_max <- 0.1
#' p_cellType_gene(e_cellType_gene = e_cellType_gene
#'               ,e_cellType_gene_max = e_cellType_gene_max
#'              ,verbose = TRUE
#' )
#' 
p_cellType_gene <- function(e_cellType_gene, e_cellType_gene_max, verbose = FALSE) {
        if (verbose) {
                print("calculate p_cellType_gene")
        }
        
        # normalize the means for the cell type of interest to max means
        e_norm <- norm_to_max(param = e_cellType_gene, param_max = e_cellType_gene_max)
        return(e_norm)
        
}