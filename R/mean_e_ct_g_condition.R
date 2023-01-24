#' @title mean_e_ct_g_condition
#'
#' @description calculates the mean for the parameters e_s_l and e_r_r for the case and control cohorts.
#'
#' @param comm_result: list containing information about the interactions, samples, and their annotations.
#'
#' @return comm_result: The dataframe anno_interactions gets following new columns: 
#' - mean_e_s_l_control,
#' - mean_e_s_l_case,
#' - mean_e_r_r_control,
#' - mean_e_r_r_case,
#'
#' @export
#' @examples
#' # load example data
#' data('comm_result')
#' # calculate mean_e_ct_g_condition
#' comm_result <- mean_e_ct_g_condition(comm_result)
#' 
mean_e_ct_g_condition <- function(comm_result) {
        
        anno_samples <- comm_result$anno_samples
        
        # define idx for control and for case
        if (!("case_or_control" %in% colnames(anno_samples))) {
                stop("anno_samples does not contain a column named 'case_or_control'. Please add this column to anno_samples and make sure that it contains the values 'case' and 'control'.")
        }
        
        idx_control <- anno_samples$case_or_control == "control"
        idx_case <- anno_samples$case_or_control == "case"
        
        # calculate mean and SD for expression and nr cells
        for (param in c("e_s_l", "e_r_r")) {
                
                
                if (!("sample_ID" %in% colnames(anno_samples))) {
                        stop("anno_samples does not contain a column named 'sample_ID'. Please add this column to anno_samples and make sure that it contains the values identical in column names in x.")
                }
                
                # check correctness
                if (!all(identical(colnames(comm_result[[param]]), anno_samples$sample_ID))) {
                        stop("Sample IDs are NOT identical in x (column names) and anno_samples$sample_ID.")
                }
                
                # check if case or control group contains only one sample
                ifelse((sum(idx_control) == 1), {
                        data_control <- as.data.frame(comm_result[[param]][, idx_control])
                        rownames(data_control) <- rownames(comm_result[[param]])
                        colnames(data_control) <- colnames(comm_result[[param]])[idx_control]
                }, {
                        data_control <- comm_result[[param]][, idx_control]
                })
                ifelse((sum(idx_case) == 1), {
                        data_case <- as.data.frame(comm_result[[param]][, idx_case])
                        rownames(data_case) <- rownames(comm_result[[param]])
                        colnames(data_case) <- colnames(comm_result[[param]])[idx_case]
                }, {
                        data_case <- comm_result[[param]][, idx_case]
                })
                
                # mean conrtol
                comm_result$anno_interactions$mean <- rowMeans(data_control)
                colnames(comm_result$anno_interactions)[colnames(comm_result$anno_interactions) ==
                                                                "mean"] <- paste0("mean_", param, "_control")
                
                # mean case
                comm_result$anno_interactions$mean <- rowMeans(data_case)
                colnames(comm_result$anno_interactions)[colnames(comm_result$anno_interactions) ==
                                                                "mean"] <- paste0("mean_", param, "_case")
                
        }
        
        return(comm_result)
}