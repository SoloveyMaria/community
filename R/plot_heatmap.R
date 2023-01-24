#' @title plot_heatmap
#' 
#' @description plots heatmap of chosen parameter
#' 
#' @param comm_result A dataframe of common interactions
#' @param which_interactions A logical vector with length of all interactions, indicating which interactions to include in the heatmap. NULL if top and by_param are set or "all"
#' @param by_param  A string indicating which column of the common interactions dataframe to use for filtering interactions, e.g. "passed_QC_filter", "passed_FDR_threshold", "passed_log2FC_threshold", "sign"
#' @param values_to_plot A string indicating which column of the common interactions dataframe to use as the values in the heatmap. e.g. "interactions_weights", "expr_l_s_active", "expr_r_r_active", "nr_l_s_active", "nr_r_r_active", "phi", "phi_l_s", "phi_r_r", "p", "p_l_s", "p_r_r"
#' @param row_font_size Font size for the row labels
#' @param color_case A string indicating the color for case samples in the heatmap
#' @param color_control A string indicating the color for control samples in the heatmap
#' 
#' @export
#' @examples
#' # load example data
#' data("comm_result")
#' # calculate general statistics
#' comm_result <- general_stat(comm_result)
#' # plot heatmap
#' plot_heatmap(comm_result = comm_result
#'             ,which_interactions = NULL
#'            ,by_param = "passed_QC_filter"
#'           ,values_to_plot = "interactions_weights"
#'         ,row_font_size = 8
#'      ,color_case = "#7C001F" # "darkred"
#'  ,color_control = "#7AC5CD" # "CadetBlue3"
#' )
#'
plot_heatmap <- function(comm_result
                         ,which_interactions = NULL 
                         # NULL if top and by_param are set 
                         # or "all" 
                         # or boolean vector with length of all interactions
                         # e.g. "ADAM10" %in% Single_Cell_Result$anno_interactions$ligand_gene_name 
                         # or "T" %in% Single_Cell_Result$anno_interactions$sending_cell_type
                         #,top = NULL 
                         ,by_param 
                         # passed_QC_filter
                         # passed_FDR_threshold
                         # passed_log2FC_threshold
                         # sign
                         ,values_to_plot 
                         # weights, 
                         # expr_l_s_active, 
                         # expr_r_r_active, 
                         # nr_l_s_active, 
                         # nr_r_r_active, 
                         # phi, 
                         # phi_l_s, 
                         # phi_r_r, 
                         # p, 
                         # p_l_s, 
                         # p_r_r
                         ,row_font_size = 8
                         ,color_case = "#7C001F" # "darkred"
                         ,color_control = "#7AC5CD" # "CadetBlue3"
                         ,color_values = NULL # or 
){
        
        df <- as.matrix(comm_result[[values_to_plot]])
        #print(str(df))
        if(is.null(which_interactions)){
                idx_interactions <- comm_result$anno_interactions[,by_param] & (!is.na(comm_result$anno_interactions[,by_param]))
        } else if(class(which_interactions) == 'logical'){
                idx_interactions <- which_interactions
        } else if(which_interactions == "all"){
                idx_interactions <- rep(TRUE
                                        ,nrow(df))
        } else {
                stop("ERROR: parameter which_interactions can be either NULL or 'all' or a logical vector.")
        }
        #print(str(idx_interactions))
        
        col_samples <- sapply(comm_result$anno_samples$case_or_control ## add to object; create column with "case" and "control"
                              ,function(i){
                                      ifelse(grepl("case"
                                                   ,i)
                                             ,color_case
                                             ,color_control
                                      )
                              })
        #print(str(col_samples))
        
        if(is.null(color_values)){
                my_color <- colorRamp2(seq(0, 1, length = 5), c("white", "red", "red4",  "darkred", "black"))
        } else my_color <-  color_values
        
        h <- Heatmap(df[idx_interactions,]
                     ,name = values_to_plot
                     ,row_names_gp = grid::gpar(fontsize = row_font_size)
                     ,col = my_color
                     ,column_names_side = "top"
                     ,column_names_gp = gpar(col = col_samples)
                     ,heatmap_legend_param = list(direction = "horizontal"
                     )
        )
        plot(h
             ,heatmap_legend_side = "bottom")
}