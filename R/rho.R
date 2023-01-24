#' @title rho
#' 
#' @description calculates the normalized abundance component of the weight formula.
#' 
#' @param rho_s numeric [0,1]: normalized abundance of the sending cell type within the sample of interest.
#' @param rho_r numeric [0,1]: normalized abundance of the receiving cell type within the sample of interest.
#' @param verbose Logical, whether to print additional information
#' 
#' @return rho numeric [0,1]: the normalized abundance component of the weight formula
#' 
#' @export
#' 
rho <- function(rho_s, rho_r, verbose = FALSE) {
        
        if (verbose) {
                print("calculate rho")
        }
        
        rho_s * rho_r
}