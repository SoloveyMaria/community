#' @title phi_cellType_gene
#' 
#' @description calculates the relative active fraction of a ligand or receptor based on its expression level in the cell type of interest.
#' 
#' @param a_cellType_gene numeric [0,1]: active fraction data for the ligand (or receptor) in the cell type of interest.
#' @param a_cellType_gene_max numeric [0,1]: maximum active fraction for the ligand (or receptor) in the cell type of interest across all samples.
#' @param verbose a logical value indicating whether to print progress messages (default: FALSE).
#' 
#' @return a_norm: numeric [0,1]: normalised active fraction for the ligand (or receptor) in the cell type of interest.
#' 
#' @export
#' @examples
#' a_cellType_gene <- 0.1
#' a_cellType_gene_max <- 0.1
#' phi_cellType_gene(a_cellType_gene = a_cellType_gene
#'                  ,a_cellType_gene_max = a_cellType_gene_max
#'                 ,verbose = TRUE
#' )
phi_cellType_gene <- function(a_cellType_gene, a_cellType_gene_max, verbose = FALSE) {
        if (verbose) {
                print("calculate phi (relative active fraction) for each cell type")
        }
        
        # normalize the active fraction for the cell type of interest to max active
        # fraction
        a_norm <- norm_to_max(param = a_cellType_gene, param_max = a_cellType_gene_max)
        return(a_norm)
}