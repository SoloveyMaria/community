#' @title rho_cellType
#' 
#' @description calculates relative abundance of cell types.
#' 
#' @param f_cellType numeric [0,1]: abundance of the cell types (as fractions).
#' @param f_cellType_max numeric [0,1]: maximum abundance of the cell type of interest across all samples.
#' @param verbose a logical value indicating whether to print progress messages (default: FALSE).
#' 
#' @return relative abundance for each cell type.
#' 
#' @export
#' 
rho_cellType <- function(f_cellType, f_cellType_max, verbose = FALSE) {
        
        if (verbose) {
                print("calculate rho (relative abundance) for each cell type")
        }
        
        # calculate relative abundance for each cell type
        f_norm <- norm_to_max(param = f_cellType, param_max = f_cellType_max)
        return(f_norm)
}