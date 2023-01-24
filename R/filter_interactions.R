#' @title filter_interactions
#'
#' @description The `filter_interactions` function produces one plot for the **quality filter** and two plots for the **discrepancy filter** (one for the controls and one for the cases). It writes the selected threshold values in the `thresholds` list of the interaction object and stores the filtering results as boolean vectors (one per threshold) in the `anno_interactions` list.
#' 
#' @param comm_result: list containing information about the interactions, samples, and their annotations.
#' @param threshold_log10_cum_weight numeric: the threshold for the log10 cumulative weight of interactions. The default value is 0.05.
#' @param threshold_nr_expr_samples numeric: the threshold for the number of samples in which an interaction has a non-zero value. The default value is 2.
#' @param threshold_log10_meanexpr_per_condition numeric: the threshold for the log10 mean expression per condition. The default value is 0.01.
#' @param verbose a logical value indicating whether to print progress messages (default: FALSE).
#'
#' @return comm_results: The dataframe anno_interactions gets following new columns: 
#' - anno_interactions
#' -- log10_cum_weight: numeric: the log10 cumulative weight of an interaction
#' -- nr_expr_samples: numeric: the number of samples in which an interaction has a non-zero value
#' -- passed_log10_cum_weight_filter: a logical value indicating whether an interaction passed the log10_cum_weight filter
#' -- passed_nr_expr_samples_filter: a logical value indicating whether an interaction passed the nr_expr_samples filter
#' -- passed_log10_meanexpr_control_filter: a logical value indicating whether an interaction passed the threshold_log10_meanexpr_per_condition filter in the control samples
#' -- passed_log10_meanexpr_case_filter: a logical value indicating whether an interaction passed the threshold_log10_meanexpr_per_condition filter in the case samples
#' -- passed_log10_meanexpr_per_condition_filter: a logical value indicating whether an interaction passed the threshold_log10_meanexpr_per_condition filter. An interaction passes this filter if both its ligand and receptor pass the threshold either in control samples or in case samples or in both.
#' -- passed_QC_filter: a logical value indicating whether an interaction passed all the QC filters.
#' 
#' - thresholds
#' -- threshold_log10_cum_weight: numeric: the threshold for the log10 cumulative weight of interactions. The default value is 0.05.
#' -- threshold_nr_expr_samples: numeric: the threshold for the number of samples in which an interaction has a non-zero value. The default value is 2.
#' -- threshold_log10_meanexpr_per_condition: numeric: the threshold for the log10 mean expression per condition. The default value is 0.01.
#'   
#' @export

filter_interactions <- function(comm_result, threshold_log10_cum_weight = 0.05, threshold_nr_expr_samples = 2,
                                threshold_log10_meanexpr_per_condition = 0.1, verbose = TRUE) {
        # plot the distribution of log10 cumucative interactions weight over the
        # number of samples in which the interactions is expressed
        
        # calculate log10 cumulative interactions weights
        comm_result$anno_interactions$log10_cum_weight <- log10(rowSums(comm_result$weights) + 1)
        
        # calculate number of samples exresssion the interactions
        comm_result$anno_interactions$nr_expr_samples <- rowSums(comm_result$weights != 0)
        
        # set thresholds
        comm_result$thresholds$threshold_log10_cum_weight <- threshold_log10_cum_weight
        comm_result$thresholds$threshold_nr_expr_samples = threshold_nr_expr_samples
        comm_result$thresholds$threshold_log10_meanexpr_per_condition = threshold_log10_meanexpr_per_condition
        
        
        plot_cumW_vs_nrSamples(df = comm_result$anno_interactions, threshold_log10_cum_weight = threshold_log10_cum_weight,
                               threshold_nr_expr_samples = threshold_nr_expr_samples)
        
        
        plot_meanLig_vs_meanRec(comm_result$anno_interactions, threshold_log10_meanexpr_per_condition = threshold_log10_meanexpr_per_condition)
        
        # filter interactions which did not pass the threshold in any sample
        comm_result$anno_interactions$passed_log10_cum_weight_filter <- comm_result$anno_interactions$log10_cum_weight >
                threshold_log10_cum_weight
        comm_result$anno_interactions$passed_nr_expr_samples_filter <- comm_result$anno_interactions$nr_expr_samples >
                threshold_nr_expr_samples
        comm_result$anno_interactions$passed_log10_meanexpr_control_filter <- (log10(comm_result$anno_interactions$mean_e_s_l_control +
                                                                                   1) > threshold_log10_meanexpr_per_condition) & (log10(comm_result$anno_interactions$mean_e_r_r_control +
                                                                                                                                                 1) > threshold_log10_meanexpr_per_condition)
        
        
        comm_result$anno_interactions$passed_log10_meanexpr_case_filter <- (log10(comm_result$anno_interactions$mean_e_s_l_case +
                                                                                1) > threshold_log10_meanexpr_per_condition) & (log10(comm_result$anno_interactions$mean_e_r_r_case +
                                                                                                                                              1) > threshold_log10_meanexpr_per_condition)
        
        comm_result$anno_interactions$passed_log10_meanexpr_per_condition_filter <- comm_result$anno_interactions$passed_log10_meanexpr_control_filter |
                comm_result$anno_interactions$passed_log10_meanexpr_case_filter
        
        # filter anno_interactions
        comm_result$anno_interactions$passed_QC_filter <- (comm_result$anno_interactions$passed_log10_cum_weight_filter &
                                                                   comm_result$anno_interactions$passed_nr_expr_samples_filter & comm_result$anno_interactions$passed_log10_meanexpr_per_condition_filter)
        
        samples <- names(comm_result$per_sample_anno_interactions)
        
        if (verbose) {
                print(paste(sum(!(comm_result$anno_interactions$passed_log10_cum_weight_filter & comm_result$anno_interactions$passed_nr_expr_samples_filter)),
                            "out of", nrow(comm_result$weights), "interactions do not pass the thresholds for log10 cumulative interactions weight >",
                            threshold_log10_cum_weight, "and number of expression samples >", threshold_nr_expr_samples,
                            ". Also ", sum(!comm_result$anno_interactions$passed_log10_meanexpr_per_condition_filter),
                            " interactions didn't pass the discrepancy filter.", " In total,", sum(!comm_result$anno_interactions$passed_QC_filter),
                            " bad quality interactions will be removed and", sum(comm_result$anno_interactions$passed_QC_filter),
                            "good quality interactions will remain."))
        }
        
        return(comm_result)
}