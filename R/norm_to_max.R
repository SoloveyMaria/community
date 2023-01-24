#' @title norm_to_max
#'
#' @description normalizes the param for the cell type of interest to param_max.
#'
#' @param param a numeric vector representing a parameter of interest for normalization.
#' @param param_max a numeric vector representing the maximum value of the parameter.
#'
#' @return param_norm a numeric vector of normalized parameter values.
#'
#' @export
#' 
norm_to_max <- function(param, param_max) {
        
        # normalize the param for the cell type of interest to param_max
        param_norm <- sapply(1:length(param_max), function(i) {
                ifelse((param_max[i] == 0) | (is.na(param_max[i])), 0, param[i]/param_max[i])
        })
        return(param_norm)
}