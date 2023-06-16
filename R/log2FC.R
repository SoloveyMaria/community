# NEW:
# the minimum value over the dataset, which is used for dummy in cases of zero, is not NOT devided by 10, but kept as it is

#' @title log2FC
#' 
#' @description calculates log2FC for a given parameter between the case and the control level. To avoid division by zero, all zero values are substituted by a pseudo-count. The pseudo-count is calculated as the minimum value in the x divided by 10.
#'
#' @param x numeric dataframe: can be a dataframe with interactions weights, phi, p values, etc. The rows are interaction_IDs, the columns are sample IDs. The column names should be identical to the row names in anno_samples.
#' @param anno_samples dataframe: should contain the columns: sample_ID and case_or_control. The row names should be defined and equal to sample_ID.
#' @param verbose a logical value indicating whether to print progress messages (default: FALSE).
#' 
#' @return numeric [-Inf, Inf]: log2FC for each value of the submitted parameter.
#' 
#' @export
#' @examples
#' # load example data
#' data("comm_result")
#' # calculate general statistics
#' comm_result <- general_stat(comm_result)
#' # calculate log2FC for interactions weights
#' log2FC(x = comm_result$weights
#'       ,anno_samples = comm_result$anno_samples
#'      ,verbose = TRUE
#' )
#' # calculate log2FC for phi
#' log2FC(x = comm_result$phi
#'      ,anno_samples = comm_result$anno_samples
#'     ,verbose = TRUE
#' )
log2FC <- function(x
                   ,anno_samples
                   ,verbose = FALSE
){
        
        if(!("sample_ID" %in% colnames(anno_samples))
        ){stop(
                "anno_samples does not contain a column named 'sample_ID'. Please add this column to anno_samples and make sure that it contains the values identical in column names in x."
        )
        }
        
        # check correctness
        if(!all(identical(colnames(x)
                          ,anno_samples$sample_ID
        )
        )
        ){stop(
                "Sample IDs are NOT identical in x (column names) and anno_samples$sample_ID."
        )
        }
        
        if(!("case_or_control" %in% colnames(anno_samples))
        ){stop(
                "anno_samples does not contain a column named 'case_or_control'. Please add this column to anno_samples and make sure that it contains the values 'case' and 'control'."
        )
        }
        
        # find min value
        min_value <- min(x[x != 0])
        if(verbose){
                print("min value in x is")
                print(min_value)
        }
        
        # define min count 
        min_count <- min_value 
        if(verbose){
                print("min pseudo-count to substitute zeros is")
                print(min_count)
        }
        
        # substitute zero with min_count
        x[x == 0] <- min_count
        
        # define idx control and case
        idx_control <- anno_samples$case_or_control == "control"
        idx_case <- anno_samples$case_or_control == "case"
        
        # calculate mean values for control and case
        df_means <- data.frame(control = rowMeans(x[,idx_control])
                               ,case = rowMeans(x[,idx_case])
        )
        log2FC <- log2(df_means$case / df_means$control)
        
        return(log2FC)
}