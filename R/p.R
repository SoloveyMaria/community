#' @title p
#'
#' @description calculates the expression based probability of communication
#' 
#' @param p_s_l numeric [0,1]: probability of interaction of a ligand with its counterpart receptor based on the relative mean expression of this ligand in the sending cell type. 
#' @param p_r_r numeric [0,1]: probability of interaction of a receptor with its counterpart ligand based on the relative mean expression of this receptor in the receiving cell type.
#' @param verbose a logical value indicating whether to print progress messages (default: FALSE).
#'  
#' @return numeric [0,1]: probability of interaction based on the relative mean expression of the ligand and receptor in the sending and receiving cell types. 
#' 
#' @export
#' @examples
#' p_s_l <- 0.1
#' p_r_r <- 0.1
#' p_expr(p_s_l = p_s_l
#'      ,p_r_r = p_r_r
#'     ,verbose = TRUE
#' )
#' 
p <- function(p_s_l, p_r_r, verbose = FALSE) {
        if (verbose) {
                print("calculate p")
        }
        
        p_s_l * p_r_r
}