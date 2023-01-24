#' @title plot_cumW_vs_nrSamples
#'
#' @description creates a scatter plot of log10 cumulative weight of interactions vs number of samples in which the interactions are expressed. It also gives a visual representation of whether the interactions passed a certain threshold for log10 cumulative weight and number of samples expressed in.
#'
#' @param df: dataframe: dataframe containing the interactions, including columns "log10_cum_weight" and "nr_expr_samples"
#' @param threshold_log10_cum_weight: numeric: threshold for log10 cumulative weight of interactions. Interactions with a log10 cumulative weight above this threshold will be considered as passing the threshold.
#' @param threshold_nr_expr_samples: numeric: threshold for the number of samples in which the interactions are expressed. Interactions expressed in more than this number of samples will be considered as passing the threshold.
#'
#' @return list of ggplots: list containing the scatter plot, the y-density plot, and the x-density plot
#'
#' @export
#' @examples
#' # plot the cumulative weight vs number of samples in which the interactions are expressed
#' plot_cumW_vs_nrSamples(df = anno_interactions_allSamples[[1]]
#'                       ,threshold_log10_cum_weight = 0
#'                      ,threshold_nr_expr_samples = 0
#' )
####### is this returning a list of ggplots?? because it's not clear from the code below


plot_cumW_vs_nrSamples <- function(df
                                   ,threshold_log10_cum_weight
                                   ,threshold_nr_expr_samples
){
        # check which interactions passed the threshold
        df$passed_threshold <- (df$log10_cum_weight > threshold_log10_cum_weight) & (df$nr_expr_samples>threshold_nr_expr_samples)
        
        # main plot
        QC_plot <- ggplot(data = df
                          ,aes(x = log10_cum_weight
                               ,y = nr_expr_samples
                               ,color = passed_threshold
                               ,shape = passed_threshold
                          )
        )+
                
                geom_point(alpha = 0.5
                )+
                scale_color_manual(labels = c("FALSE"="failed at least one threshold"
                                              , "TRUE"="passed both thresholds")
                                   ,values=c("FALSE"="gray80"
                                             ,"TRUE"="black")
                                   ,guide = guide_legend()
                )+
                scale_shape_manual(labels = c("FALSE"="failed at least one threshold"
                                              , "TRUE"="passed both thresholds")
                                   ,values=c("FALSE"=16,"TRUE"=1)
                                   ,guide = guide_legend()
                )+
                ggtitle("interactions QC")+
                geom_vline(xintercept = threshold_log10_cum_weight
                           , color = "red")+
                geom_hline(yintercept = threshold_nr_expr_samples
                           , color = "red")+
                theme_bw()+
                theme(text = element_text(size=20)
                      ,legend.position="bottom"
                      ,legend.title=element_blank())
        
        # y density plot
        ydensity <- ggplot(data = df
                           ,aes(x = nr_expr_samples)
        )+
                geom_density()+
                geom_vline(xintercept = threshold_nr_expr_samples
                           , color = "red")+
                theme(text = element_text(size=0))+
                ggtitle("")+
                xlab("")+
                coord_flip()+
                theme_bw()+
                theme(plot.margin = unit(c(0.5
                                           , 0.5
                                           , 1
                                           , 3
                )
                , "cm"))
        
        # x density plot
        xdensity <- ggplot(data = df
                           ,aes(x = log10_cum_weight))+
                geom_density()+
                geom_vline(xintercept = threshold_log10_cum_weight
                           , color = "red")+
                theme(text = element_text(size=0))+
                ggtitle("")+
                xlab("")+
                theme_bw()+
                theme(plot.margin = unit(c(0.5
                                           , 9# right
                                           , 1
                                           , 1 # left
                )
                , "cm")
                )
        
        
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
        
        # arrange plots
        grid.arrange(ydensity
                     ,QC_plot
                     ,blankPlot
                     ,xdensity
                     ,ncol=2
                     ,nrow=2
                     ,widths=c(2, 6)
                     ,heights=c(6, 2))
}
