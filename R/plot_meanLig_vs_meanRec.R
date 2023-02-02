#' @title plot_meanLig_vs_meanRec
#' 
#' @description auxillary function for interactions QC: plot log10 mean ligand expression vs log10 mean receptor expression (in controls and cases separately).
#' 
#' @param df the df contains following columns: interactions_ID, mean_e_s_l_control, mean_e_s_l_case, mean_e_r_r_control, mean_e_r_r_case.
#' @param threshold_log10_meanexpr_per_condition numeric: threshold for log10 mean expression per condition. The default values is 0.1.
#' 
#' @return plotting function
#' 
#' @export
#' @examples
#' plot_meanLig_vs_meanRec(df = x$anno_interactions
#'                        ,threshold_log10_meanexpr_per_condition = threshold_log10_meanexpr_per_condition)
#' 
plot_meanLig_vs_meanRec <- function(df,threshold_log10_meanexpr_per_condition = 0.1
){
        
        max_value <- log10(max(df[,c("mean_e_s_l_control"
                                     ,"mean_e_s_l_case"
                                     ,"mean_e_r_r_control"
                                     ,"mean_e_r_r_case")]
                               ,na.rm=TRUE)+1)
        
        df$passed_discrepancy_threshold_control <- (log10(df$mean_e_s_l_control+1) > threshold_log10_meanexpr_per_condition) & (
                log10(df$mean_e_r_r_control+1) > threshold_log10_meanexpr_per_condition)
        
        df$passed_discrepancy_threshold_case <- (log10(df$mean_e_s_l_case+1) > threshold_log10_meanexpr_per_condition) & (
                log10(df$mean_e_r_r_case+1) > threshold_log10_meanexpr_per_condition)
        
        df$passed_discrepancy_threshold <- df$passed_discrepancy_threshold_control | df$passed_discrepancy_threshold_case
        
        plot_list <- lapply(c("control"
                              ,"case")
                            ,function(condition){
                                    ifelse(condition == "control"
                                           ,{my_df <- data.frame(interaction_ID = df$interaction_ID
                                                                 ,mean_expr_s_l_active = df$mean_e_s_l_control
                                                                 ,mean_expr_r_r_active = df$mean_e_r_r_control
                                                                 ,passed_discrepancy_threshold = df$passed_discrepancy_threshold
                                           )
                                           }
                                           ,{my_df <- data.frame(interaction_ID = df$interaction_ID
                                                                 ,mean_expr_s_l_active = df$mean_e_s_l_case
                                                                 ,mean_expr_r_r_active = df$mean_e_r_r_case
                                                                 ,passed_discrepancy_threshold = df$passed_discrepancy_threshold
                                           )
                                           }
                                    )
                                    
                                    # y density plot
                                    ydensity <- ggplot(data = my_df
                                                       ,aes(x = log10(mean_expr_s_l_active+1))
                                    )+
                                            geom_density()+
                                            xlim(c(0,max_value))+
                                            xlab("")+
                                            geom_vline(xintercept = threshold_log10_meanexpr_per_condition
                                                       , color = "red")+
                                            theme(text = element_text(size=0))+
                                            ggtitle("")+
                                            coord_flip()+
                                            theme_bw()+
                                            theme(plot.margin = unit(c(0.5, 1, 2.5, 2), "cm"))
                                    
                                    # main plot
                                    QC_plot <- ggplot(data = my_df
                                                      ,aes(x = log10(mean_expr_r_r_active+1)
                                                           ,y = log10(mean_expr_s_l_active+1)
                                                           ,color = passed_discrepancy_threshold
                                                           ,shape = passed_discrepancy_threshold
                                                      )
                                    )+
                                            geom_point(alpha = 0.5
                                                       #,show.legend = FALSE
                                            )+
                                            scale_color_manual(labels = c("FALSE"="failed threshold in\nboth conditions"
                                                                          , "TRUE"="passed threshold in\nat least one condition")
                                                               ,values=c("FALSE"="gray80"
                                                                         ,"TRUE"="black")
                                                               ,guide = guide_legend()
                                            )+
                                            scale_shape_manual(labels = c("FALSE"="failed threshold in\nboth conditions"
                                                                          , "TRUE"="passed threshold in\nat least one condition")
                                                               ,values=c("FALSE"=16,"TRUE"=1)
                                                               ,guide = guide_legend()
                                            )+
                                            ggtitle(condition)+
                                            ylim(c(0,max_value))+
                                            xlim(c(0,max_value))+
                                            xlab("receptor in receiving cell type\n[log10 mean expression]")+
                                            ylab("ligand in sending cell type\n[log10 mean expression]")+
                                            geom_vline(xintercept = threshold_log10_meanexpr_per_condition
                                                       , color = "red")+
                                            geom_hline(yintercept = threshold_log10_meanexpr_per_condition
                                                       , color = "red")+
                                            theme_bw()+
                                            theme(text = element_text(size=20)
                                                  ,legend.position="bottom"
                                                  ,legend.title=element_blank()
                                            )
                                    #guides(color = guide_legend(label.position = "bottom")
                                    #       ,shape = guide_legend(label.position = "bottom")
                                    #      )+
                                    
                                    
                                    # x density plot
                                    xdensity <- ggplot(data = my_df
                                                       ,aes(x = log10(mean_expr_r_r_active+1)))+
                                            geom_density()+
                                            xlim(c(0,max_value))+
                                            xlab("")+
                                            geom_vline(xintercept = threshold_log10_meanexpr_per_condition
                                                       , color = "red")+
                                            theme(text = element_text(size=0))+
                                            ggtitle("")+
                                            theme_bw()+
                                            theme(plot.margin = unit(c(0, 0, 1.8, 0), "cm"))
                                    
                                    
                                    # Create a blank placeholder plot :
                                    blankPlot <- ggplot(data = df
                                                        ,aes(x = log10_cum_weight
                                                             ,y = nr_expr_samples
                                                        )
                                    )+
                                            geom_jitter(alpha=0) +
                                            theme(axis.line=element_blank()
                                                  ,plot.background = element_blank(), 
                                                  panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(), 
                                                  panel.border = element_blank(),
                                                  panel.background = element_blank(),
                                                  axis.title.x = element_blank(),
                                                  axis.title.y = element_blank(),
                                                  axis.text.x = element_blank(), 
                                                  axis.text.y = element_blank(),
                                                  axis.ticks = element_blank())
                                    
                                    
                                    
                                    
                                    my_l <- list(ydensity=ydensity
                                                 ,QC_plot=QC_plot
                                                 ,blankPlot=blankPlot
                                                 ,xdensity=xdensity
                                    )
                                    
                            })
        
        grid.arrange(plot_list[[1]][["ydensity"]]
                     ,plot_list[[1]][["QC_plot"]]
                     ,plot_list[[2]][["ydensity"]]
                     ,plot_list[[2]][["QC_plot"]]
                     ,plot_list[[1]][["blankPlot"]]
                     ,plot_list[[1]][["xdensity"]]
                     ,plot_list[[2]][["blankPlot"]]
                     ,plot_list[[2]][["xdensity"]]
                     ,nrow = 2
                     ,ncol = 4
                     
                     ,widths=c(2, 6, 2, 6)
                     ,heights=c(6.5, 1.5)
                     
                     ,top=textGrob("log10 mean expression in active fraction"
                                   ,gp = gpar(fontsize = 24))
        )
}