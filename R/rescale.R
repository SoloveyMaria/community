#' @title rescale
#' 
#' @description Rescales a numeric vector to range between 0 and 1
#' 
#' @param x A numeric vector
#' @return A numeric vector with all values rescaled to be between 0 and 1
#' @export
#' 
rescale <- function(x) {
        (x - min(x))/(max(x) - min(x))
}