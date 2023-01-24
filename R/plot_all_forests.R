#' @title plot_all_forests
#'
#' @description creates 7 forest plots for different parameters.
#'
#' @param my_idx Index of the interactions to be plotted
#' @param my_anno_interactions dataframe: Dataframe containing the interactions, having the rows as interactions and columns as the different parameters and the interaction_ID.
#' @param keep_order a logical value. If TRUE then keeps oder of my_anno_interactions ignoring the my_idx, default is FALSE
#' @param show_labels a logical value. If true interaction_IDs are plotted, default is FALSE
#' @param plot_legend a logical value. If true the legend is plotted, default is TRUE
#'
#' @return a list of 7 forest plots
#'
#' @export
#' @examples
#' # plot_all_forests
#' plot_all_forests(my_idx
#'                ,my_anno_interactions
#' )
#' 
plot_all_forests <- function(my_idx
                             ,my_anno_interactions
                             ,keep_order=FALSE
                             ,show_labels=FALSE
                             ,plot_legend=TRUE
){
        
        # define order
        if(!keep_order){
                my_anno_interactions <- order_interactions_for_forests(my_idx)  
        }
        
        params <- c("log2FC_weights"
                    ,"log2FC_rho_s"
                    ,"log2FC_phi_s_l"
                    ,"log2FC_p_s_l"
                    ,"log2FC_rho_r"
                    ,"log2FC_phi_r_r"
                    ,"log2FC_p_r_r"
        )
        my_data <- lapply(params
                          ,function(i){
                                  
                                  test_df <- my_anno_interactions[,c(i,"interaction_ID")]
                                  colnames(test_df) <- c("log2FC","interaction_ID")
                                  
                                  # cut the too large values
                                  idx_less_minus_5 <- test_df$log2FC < -4
                                  test_df$log2FC[idx_less_minus_5] <- -4
                                  test_df
                          }
        )
        
        names(my_data) <- params
        
        p0 <- plot_forest(my_data$log2FC_weights
                          ,my_title = "w"
                          ,min = min(my_data$log2FC_weights$log2FC)
                          ,max = max(my_data$log2FC_weights$log2FC)
                          ,plot_legend = plot_legend
        )
        
        p1 <- plot_forest(my_data$log2FC_rho_s
                          ,my_title = "rho_s"
                          ,min = min(my_data$log2FC_rho_s$log2FC)
                          ,max = max(my_data$log2FC_rho_s$log2FC)
                          ,plot_legend = plot_legend
        )
        
        p2 <- plot_forest(my_data$log2FC_phi_s_l
                          ,my_title = "phi_s_l"
                          ,min = min(my_data$log2FC_phi_s_l$log2FC)
                          ,max = max(my_data$log2FC_phi_s_l$log2FC)
                          ,plot_legend = plot_legend
        )
        p3 <- plot_forest(my_data$log2FC_p_s_l
                          ,my_title = "p_s_l"
                          ,min = min(my_data$log2FC_p_s_l$log2FC)
                          ,max = max(my_data$log2FC_p_s_l$log2FC)
                          ,plot_legend = plot_legend
        )
        p4 <- plot_forest(my_data$log2FC_rho_r
                          ,my_title = "rho_r"
                          ,min = min(my_data$log2FC_rho_r$log2FC)
                          ,max = max(my_data$log2FC_rho_r$log2FC)
                          ,plot_legend = plot_legend
        )
        
        
        p5 <- plot_forest(my_data$log2FC_phi_r_r
                          ,my_title = "phi_r_r"
                          ,min = min(my_data$log2FC_phi_r_r$log2FC)
                          ,max = max(my_data$log2FC_phi_r_r$log2FC)
                          ,plot_legend = plot_legend
        )
        p6 <- plot_forest(my_data$log2FC_p_r_r
                          ,my_title = "p_r_r"
                          ,min = min(my_data$log2FC_p_r_r$log2FC)
                          ,max = max(my_data$log2FC_p_r_r$log2FC)
                          ,plot_legend = plot_legend
        )
        
        my_data$empty_values <- data.frame(log2FC = rep(0,(nrow(my_data[[1]])+1))
                                           ,interaction_ID = c("",as.character(my_data[[1]]$interaction_ID))
        )
        my_data$empty_values$interaction_ID <- factor(my_data$empty_values$interaction_ID
                                               ,levels = c("",as.character(my_data[[1]]$interaction_ID))
                                               ,ordered = TRUE)
        
        #print(str(my_data))
        p_IDs <- ggplot(my_data$empty_values
                        ,aes(y = interaction_ID
                             ,x = log2FC
                        )
        )+ theme_classic()
        p_IDs <- p_IDs + theme(axis.line.x=element_blank(),
                               axis.text.x=element_blank(),
                               axis.ticks.x=element_blank(),
                               axis.title.x=element_blank(),
                               panel.grid.minor.x=element_blank(),
                               panel.grid.major.x=element_blank()
                               ,axis.title.y=element_blank()
                               ,axis.ticks.y=element_blank()
                               ,axis.line.y=element_blank())
        
        
        
        
        margin = theme(plot.margin = unit(c(0
                                            ,-0.25
                                            ,0
                                            ,-0.25
        )
        , "cm")
        )
        margin_middle = theme(plot.margin = unit(c(1
                                                   ,2.5
                                                   ,1
                                                   ,1
        )
        , "cm")
        )
        
        
        #options(repr.plot.width = 18)
        
        if(!show_labels){
                grid.arrange(p0
                             ,ggplot()+ theme_void()+margin
                             ,p1
                             ,p2
                             ,p3
                             ,ggplot()+ theme_void()+margin
                             ,p4
                             ,p5
                             ,p6
                             ,nrow = 1)
        }else{
                grid.arrange(p0
                             ,ggplot()+ theme_void()+margin
                             ,p1
                             ,p2
                             ,p3
                             ,ggplot()+ theme_void()+margin
                             ,p4
                             ,p5
                             ,p6
                             ,p_IDs
                             ,nrow = 1)
        }
        
        
        
}