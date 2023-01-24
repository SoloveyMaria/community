#' @title phi
#'
#' @description calculates the active fraction based probability of communication
#' 
#' @param phi_s_l numeric [0,1]: probability of interaction of a ligand with its counterpart receptor based on the relative active fraction of this ligand in the sending cell type. 
#' @param phi_r_r numeric [0,1]: probability of interaction of a receptor with its counterpart ligand based on the relative active fraction of this receptor in the receiving cell type.
#' @param verbose a logical value indicating whether to print progress messages (default: FALSE).
#'
#' @return numeric [0,1]: probability of interaction based on the relative mean expression of the ligand and receptor in the sending and receiving cell types. 
#' 
#'
#' @export
#' 
phi <- function(phi_s_l
                ,phi_r_r
                ,verbose = FALSE){
        
        if(verbose){print("calculate phi")}
        
        phi_s_l * phi_r_r
}