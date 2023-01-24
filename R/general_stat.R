#' @title general_stat
#'
#' @description calculates the log2 fold change between the mean values in the case cohort and the mean values in the control cohort for the following parameters: weights, rho_s, rho_r, rho, phi_l_s, phi_r_r, phi, p_l_s, p_r_r, and p.
#'
#' @param comm_result: list containing information about the interactions, samples, and their annotations.
#' @param verbose: a logical value indicating whether to print progress messages (default: FALSE).
#' 
#' @return comm_results: The dataframe anno_interactions gets following new columns: 
#' - log2FC_weights
#' - log2FC_phi_l_s
#' - log2FC_phi_r_r
#' - log2FC_phi
#' - log2FC_p_l_s
#' - log2FC_p_r_r
#' - log2FC_p
#'
#' @export
#' @examples
#' # load example data
#' data('comm_result')
#' # calculate general statistics
#' comm_result <- general_stat(comm_result)


general_stat <- function(comm_result, verbose = FALSE) {
        
        # calculate mean e_s_l and e_r_r for cases and for controls. These
        # parameters will be used in the QC
        comm_result <- mean_e_ct_g_condition(comm_result)
        
        # calculate log2FC for differential communication and visualization
        for (param in c("rho_s", "rho_r", "rho", "phi_s_l", "phi_r_r", "phi", "p_s_l",
                        "p_r_r", "p", "weights")) {
                # print(param)
                my_log2FC <- log2FC(x = comm_result[[param]], anno_samples = comm_result$anno_samples,
                                    verbose = verbose)
                comm_result$anno_interactions$log2FC <- my_log2FC
                colnames(comm_result$anno_interactions)[colnames(comm_result$anno_interactions) ==
                                                                "log2FC"] <- paste0("log2FC_", param)
        }
        
        return(comm_result)
}