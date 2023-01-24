#' @title extract_params_annointeractions
#' 
#' @description collects the values of a specific parameter from multiple dataframes of interactions annotation into one single dataframe.
#' 
#' @param anno_interactions a list of dataframes, where each dataframe contains the annotation of interactions for a specific sample. Each dataframe should have a column named 'interaction_ID' that serves as a unique identifier for the interactions, and the column of the parameter of interest.
#' @param param string, name of the column in the anno_interactions dataframes, containing the parameter of interest.
#' 
#' @return dataframe where the rows are interactions and the columns are samples. The values in the dataframe are the values of the parameter of interest in each sample.
#' @export
#' 
#' @examples
#' # load example data
#' data('anno_interactions_allSamples')
#' # extract interactions weights
#' extract_params_annointeractions(anno_interactions = anno_interactions_allSamples
#'                        ,param = 'weight'
#' )
#' 
extract_params_annointeractions <- function(anno_interactions, param) {
        my_param <- as.data.frame(matrix(, nrow = nrow(anno_interactions[[1]]), ncol = length(anno_interactions)))
        rownames(my_param) <- anno_interactions[[1]]$interaction_ID
        colnames(my_param) <- names(anno_interactions)
        
        for (sample in names(anno_interactions)) {
                my_param[anno_interactions[[sample]]$interaction_ID, sample] <- anno_interactions[[sample]][,
                                                                                                             param]
        }
        return(my_param)
}
extract_params_annointeractions <- function(anno_interactions, param) {
        my_param <- as.data.frame(matrix(, nrow = nrow(anno_interactions[[1]]), ncol = length(anno_interactions)))
        rownames(my_param) <- anno_interactions[[1]]$interaction_ID
        colnames(my_param) <- names(anno_interactions)
        
        for (sample in names(anno_interactions)) {
                my_param[anno_interactions[[sample]]$interaction_ID, sample] <- anno_interactions[[sample]][,
                                                                                                             param]
        }
        return(my_param)
}