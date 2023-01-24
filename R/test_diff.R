#' @title test_diff
#' 
#' @description tests if an interaction is statistically differential between the cases and the controls.
#' 
#' @param comm_result: list containing information about the interactions, samples, and their annotations.
#' @param which_adjustment method used for multiple testing correction by the p.adjust function. The default is 'fdr'.
#' @param threshold_fdr the p-value threshold for adjusted p-values. The default value is 0.1
#' @param threshold_log2FC the Log2 fold-change threshold for significant interactions. The default value is 1.
#' @param verbose a logical value indicating whether to print progress messages (default: FALSE).
#' @param ... Additional parameters passed to the t.test or the wilcoxon test function
#' 
#' @return A comm_result object with the following columns added:
#' - anno_interactions:
#' -- log2FC_weights: numeric: log2 fold change of the mean interaction weight between the cases and the controls
#' -- p.value: numeric: p -value
#' -- p.adj: numeric: adjusted p -value
#' -- passed_log2FC_threshold: a logical value indicating whether an interaction passed the passed_log2FC_threshold
#' -- passed_FDR_threshold: a logical value indicating whether an interaction passed the passed_FDR_threshold
#' -- sign: a logical value indicating whether an interaction passed both the passed_log2FC_threshold and the passed_FDR_threshold thresholds
#' 
#' - thresholds:
#' -- threshold_log2FC: numeric: the Log2 fold-change threshold for significant interactions.
#' -- threshold_fdr: numeric: the p-value threshold for adjusted p-values.
#' 
#' @export
#' @examples
#' # load example data
#' data("comm_result")
#' # calculate general statistics
#' comm_result <- general_stat(comm_result)
#' # test for differential expression
#' comm_result <- test_diff(comm_result
#'                        ,which_adjustment = "fdr"
#'                       ,verbose = TRUE
#'                      ,threshold_fdr = 0.1
#'                    ,threshold_log2FC = 1
#' )
#'
test_diff <- function(comm_result
                      ,which_adjustment = "fdr"
                      ,threshold_fdr = 0.1 
                      ,threshold_log2FC = 1
                      ,verbose = FALSE
                      ,... # params for the t.test of wilcoxon function
){
        anno_samples <- comm_result$anno_samples
        
        # define idx for control and for case
        if(!("case_or_control" %in% colnames(anno_samples))
        ){stop(
                "anno_samples does not contain a column named 'case_or_control'. Please add this column to anno_samples and make sure that it contains the values 'case' and 'control'."
        )
        }
        
        idx_control <- anno_samples$case_or_control == 'control'
        idx_case <- anno_samples$case_or_control == 'case'
        
        interactions <- comm_result$weights
        anno_interactions <- comm_result$anno_interactions
        
        # calculate p.values
        anno_interactions <- cbind(anno_interactions
                                   ,do.call(rbind.data.frame
                                            ,lapply(1:nrow(interactions)
                                                    ,function(i){
                                                            
                                                            # do test
                                                            if(anno_interactions$passed_QC_filter[i]){
                                                                    test <- wilcox.test(as.numeric(interactions[i,idx_case])
                                                                                        ,y = as.numeric(interactions[i,idx_control])
                                                                                        #,conf.int = TRUE
                                                                                        ,exact=FALSE)
                                                                    data.frame(p.value = test$p.value)
                                                            } else data.frame(p.value = NA)
                                                            
                                                    })
                                   )
        )
        
        # FDR adjustment
        anno_interactions$p.adj <- p.adjust(anno_interactions$p.value
                                            ,method = which_adjustment)
        
        # check significance
        anno_interactions$passed_FDR_threshold <- anno_interactions$p.adj < threshold_fdr
        anno_interactions$passed_log2FC_threshold <- abs(anno_interactions$log2FC_weights) > threshold_log2FC
        anno_interactions$sign <- anno_interactions$passed_FDR_threshold & anno_interactions$passed_log2FC_threshold
        
        if(verbose){print(paste("We have"
                                ,sum(anno_interactions$sign
                                     ,na.rm = TRUE)
                                ,"dignificantly differential interactions"
        )
        )
        }
        
        comm_result$anno_interactions <- anno_interactions
        comm_result$thresholds$threshold_fdr <- threshold_fdr
        comm_result$thresholds$threshold_log2FC <- threshold_log2FC
        
        
        return(comm_result)
}