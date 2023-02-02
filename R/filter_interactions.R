#' @title filter_interactions
#'
#' @description The `filter_interactions` function produces one plot for the **quality filter** and two plots for the **discrepancy filter** (one for the controls and one for the cases). It writes the selected threshold values in the `thresholds` list of the interaction object and stores the filtering results as boolean vectors (one per threshold) in the `anno_interactions` list.
#' 
#' @param comm_result: list containing information about the interactions, samples, and their annotations.
#' @param threshold_log10_cum_weight numeric: the threshold for the log10 cumulative weight of interactions. The default value is 0.05.
#' @param threshold_frac_samples_per_condition numeric: the threshold for the fraction of samples in which an interaction has a non-zero value. The threshold is applied separately for controls and the cases. An interaction passes the filter if it passes this threshold either for the control samples, or for the case samples, or for both. The default value is 0.8.
#' @param threshold_log10_meanexpr_per_condition numeric: the threshold for the log10 mean expression per condition. The default value is 0.01.
#' @param verbose a logical value indicating whether to print progress messages (default: FALSE).
#'
#' @return comm_results: The dataframe anno_interactions gets following new columns: 
#' - anno_interactions
#' -- log10_cum_weight: numeric: the log10 cumulative weight of an interaction
#' -- frac_samples_controls: numeric: the fraction of control samples in which an interaction has a non-zero value
#' -- frac_samples_cases: numeric: the fraction of case samples in which an interaction has a non-zero value
#' -- passed_log10_cum_weight_filter: a logical value indicating whether an interaction passed the log10_cum_weight filter
#' -- passed_frac_samples_filter: a logical value indicating whether an interaction passed the frac_samples filter
#' -- passed_log10_meanexpr_control_filter: a logical value indicating whether an interaction passed the threshold_log10_meanexpr_per_condition filter in the control samples
#' -- passed_log10_meanexpr_case_filter: a logical value indicating whether an interaction passed the threshold_log10_meanexpr_per_condition filter in the case samples
#' -- passed_log10_meanexpr_per_condition_filter: a logical value indicating whether an interaction passed the threshold_log10_meanexpr_per_condition filter. An interaction passes this filter if both its ligand and receptor pass the threshold either in control samples or in case samples or in both.
#' -- passed_QC_filter: a logical value indicating whether an interaction passed all the QC filters.
#' 
#' - thresholds
#' -- threshold_log10_cum_weight: numeric: the threshold for the log10 cumulative weight of interactions. The default value is 0.05.
#' -- threshold_frac_samples_per_condition: numeric: the threshold for the fraction of samples in which an interaction has a non-zero value. The threshold os applied separately for controls and the cases. The default value is 0.8.
#' -- threshold_log10_meanexpr_per_condition: numeric: the threshold for the log10 mean expression per condition. The default value is 0.01.
#'   
#' @export

filter_interactions <- function(comm_result, threshold_log10_cum_weight = 0.05, threshold_frac_samples_per_condition = 0.8,
                                threshold_log10_meanexpr_per_condition = 0.1, verbose = TRUE) {
        # plot the distribution of log10 cumulative interactions weight over the
        # fraction of samples in which the interactions is expressed
        
        # calculate log10 cumulative interactions weights
        comm_result$anno_interactions$log10_cum_weight <- log10(rowSums(comm_result$weights) + 1)
        
        # identify control samples
        idx_control <- comm_result$anno_samples$case_or_control == "control"
        idx_case <- comm_result$anno_samples$case_or_control == "case"
        
        # calculate the fraction of samples expressing the interactions
        comm_result$anno_interactions$frac_samples_controls <- rowSums(comm_result$weights[,idx_control] != 0) / sum(idx_control) 
        comm_result$anno_interactions$frac_samples_cases <- rowSums(comm_result$weights[,idx_case] != 0) / sum(idx_case) 
        
        # set thresholds
        comm_result$thresholds$threshold_log10_cum_weight <- threshold_log10_cum_weight
        comm_result$thresholds$threshold_frac_samples_per_condition = threshold_frac_samples_per_condition
        comm_result$thresholds$threshold_log10_meanexpr_per_condition = threshold_log10_meanexpr_per_condition
        
        
        cumW <- plot_cumW(df = comm_result$anno_interactions, threshold_log10_cum_weight = threshold_log10_cum_weight)
        fracSamp <- plot_fracSamples(df = comm_result$anno_interactions, threshold_frac_samples_per_condition = threshold_frac_samples_per_condition)
        
        p <- arrangeGrob(fracSamp$ydensity
                                     ,fracSamp$QC_plot
                                     ,fracSamp$blankPlot
                                     ,fracSamp$xdensity
                                     ,ncol=2
                                     ,nrow=2
                                     ,widths=c(2.5, 5.5)
                                     ,heights=c(6.5, 1.5)
                        )
    
        # arrange plots
        grid.arrange(cumW
                     , p
                     , ncol=2
                     , widths = c(3.5,4.5)
                    )
        
        plot_meanLig_vs_meanRec(comm_result$anno_interactions, threshold_log10_meanexpr_per_condition = threshold_log10_meanexpr_per_condition)
        
        # filter interactions which did not pass the threshold in any sample
        comm_result$anno_interactions$passed_log10_cum_weight_filter <- comm_result$anno_interactions$log10_cum_weight >
                threshold_log10_cum_weight
        comm_result$anno_interactions$passed_frac_samples_filter <- (comm_result$anno_interactions$frac_samples_controls >
                                                                                  threshold_frac_samples_per_condition) | (comm_result$anno_interactions$frac_samples_cases > threshold_frac_samples_per_condition)
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
                                                                   comm_result$anno_interactions$passed_frac_samples_filter & comm_result$anno_interactions$passed_log10_meanexpr_per_condition_filter)
        
        samples <- names(comm_result$per_sample_anno_interactions)
        
        if (verbose) {
                print(paste(sum(!(comm_result$anno_interactions$passed_log10_cum_weight_filter & comm_result$anno_interactions$passed_frac_samples_filter)),
                            "out of", nrow(comm_result$weights), "interactions do not pass the thresholds for log10 cumulative interactions weight >",
                            threshold_log10_cum_weight, "and fraction of expressing samples >", threshold_frac_samples_per_condition,
                            ". Also ", sum(!comm_result$anno_interactions$passed_log10_meanexpr_per_condition_filter),
                            " interactions didn't pass the discrepancy filter.", " In total,", sum(!comm_result$anno_interactions$passed_QC_filter),
                            " bad quality interactions will be removed and", sum(comm_result$anno_interactions$passed_QC_filter),
                            "good quality interactions will remain."))
        }
        
        return(comm_result)
}