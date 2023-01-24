#' @title w
#' 
#' @description calculates the weight of the interactions between two nodes based on fraction of cells and expression based probability of interaction.
#' 
#' @param rho numeric [0,1]: relative fraction of cell types.
#' @param phi numeric [0,1]: probability of interaction between two cell types based on fraction of communication cells.
#' @param p numeric [0,1]: probability of interaction between two cell types based on expression levels of the ligand and receptor.
#' @param verbose a logical value indicating whether to print progress messages (default: FALSE).
#' 
#' @return numeric [0,1]: weight of the interaction.
#' 
#' @export
#'
w <- function(rho, phi, p, verbose = FALSE) {
        if (verbose)
                print("calculate interactions weights")
        rho * phi * p
}