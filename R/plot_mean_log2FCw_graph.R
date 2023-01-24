#' @title plot_mean_log2FCw_graph
#' 
#' @description creates a network graph of cell types interactions, where each interactions represents the mean of log2FC weights of interactions between two cell types.
#' 
#' @param my_interactions A list containing information about the interactions, samples, and their annotations.
#' @param my_idx A logical vector of the same length as number of interactions in my_interactions. Only interactions with TRUE values will be plotted.
#' @param direction A character string, "down" or "up", indicating the direction of the plot (up=upregulate, down=downregulated).
#' 
#' @return A graph of cell types interactions.
#' 
#' @export
#' @examples
#' # plot_mean_log2FCw_graph
#' plot_mean_log2FCw_graph(my_interactions
#'                  ,my_idx
#'                 ,"up"
#' )
#' 
plot_mean_log2FCw_graph <- function(my_interactions
                                    ,my_idx
                                    ,direction){
        
        cell_types <- unique(my_interactions$anno_cells$cell_type)
        health_status <- unique(my_interactions$anno_samples$health_status)
        
        dummy_min <- 0.00000001
        
        # mean weight of interactions. Sign 
        # prepare adjusency matrix
        mean_mat <- matrix(NA
                           ,nrow = length(cell_types)
                           ,ncol = length(cell_types))
        rownames(mean_mat) <- cell_types
        colnames(mean_mat) <- cell_types
        
        # calculate mean weight per pair of sending and receiving cell types and per condition
        for(i in cell_types){
                for(j in cell_types){
                        idx_send <- my_interactions$anno_interactions$sending_cell_type == i
                        idx_rec <- my_interactions$anno_interactions$receiving_cell_type == j
                        
                        idx_interactions <- my_idx & idx_send & idx_rec
                        
                        ifelse(sum(idx_interactions) != 0
                               ,{means <- mean(my_interactions$anno_interactions$log2FC_weights[idx_interactions]
                               )
                               mean_mat[i,j]  <- mean(means[means != 0])
                               }
                               ,mean_mat[i,j]  <- 0
                        )
                        
                        
                        
                }
        }
        
        interactions.width <- unlist(lapply(1:nrow(mean_mat),function(i) mean_mat[i,]))
        
        
        max_value <- max(ceiling(abs(interactions.width)))
        
        ifelse(direction == "down"
               ,color_palette  <-  c("white","blue","blue3","black")
               ,color_palette  <-  c("white","red","red3","black")
        )
        
        col <- colorRampPalette(color_palette)(max_value)
        
        # construct graph
        g <- graph(c(t(my_interactions$anno_interactions[ ,c("sending_cell_type","receiving_cell_type")]
        )))
        
        g_simp <- igraph::simplify(g, remove.loops = F)
        
        options(repr.plot.width = 7
                ,repr.plot.height = 7)
        
        
        my_interactions.width <- abs(interactions.width) + dummy_min
        
        p <- plot(g_simp
                  ,layout = layout_in_circle(g_simp)
                  ,interactions.curved=0.05
                  ,vertex.shape="none"
                  ,vertex.label.cex = 2
                  ,edge.color = col[round(my_interactions.width,digits = 0)+1]
                  ,edge.width = my_interactions.width 
        )
        
        p
        
        return(p)
}