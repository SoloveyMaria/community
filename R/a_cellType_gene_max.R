#' @title a_cellType_gene_max
#' 
#' @description calculates the maximum active fraction for ligand in sending cell type and receptor in receiving cell type for each gene over all samples. 
#'
#' @param anno_interactions_allSamples a list of data frames containing annotation information for interactions in all samples.
#' @param verbose a logical value indicating whether to print progress messages (default: FALSE).
#' 
#' @return a list of data frames (numeric [0,1]) with an additional column for the maximum active fraction for ligand in sending cell type and another for the maximum active fraction for receptor in receiving cell type.
#' 
#' @export
#' @examples
#' # load example data
#' data(example_anno_interactions_allSamples)
#' 
a_cellType_gene_max <- function(anno_interactions_allSamples, verbose = FALSE) {
        # calculate max active fraction for ligand in sending cell type
        if (verbose) {
                print("calculate max active fraction for ligand in sending cell type")
        }
        a_s_l_max <- apply(sapply(names(anno_interactions_allSamples)  # for each sample
                                  ,
                                  function(my_sample) {
                                          anno_interactions_allSamples[[my_sample]][, "a_s_l"]  # extract the vector of values of a_s_l (for all interactions)
                                  }), 1, max)
        for (my_sample in names(anno_interactions_allSamples)) {
                anno_interactions_allSamples[[my_sample]]$a_s_l_max <- a_s_l_max
        }
        
        # calculate max active fraction for receptor in receiving cell type
        if (verbose) {
                print("calculate max active fraction for receptor in receiving cell type")
        }
        a_r_r_max <- apply(sapply(names(anno_interactions_allSamples)  # iterate over all samples
                                  ,
                                  function(my_sample) {
                                          anno_interactions_allSamples[[my_sample]][, "a_r_r"]
                                  }), 1, max)
        for (my_sample in names(anno_interactions_allSamples)) {
                anno_interactions_allSamples[[my_sample]]$a_r_r_max <- a_r_r_max
        }
        
        return(anno_interactions_allSamples)
        
}