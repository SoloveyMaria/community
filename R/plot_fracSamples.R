#' @title plot_fracSamples
#'
#' @description creates a scatter plot of fraction of samples in which the interactions are expressed. It also gives a visual representation of whether the interactions passed a certain threshold for log10 cumulative weight and number of samples expressed in.
#'
#' @param df: dataframe: dataframe containing the interactions, including columns "frac_samples_controls" and "frac_samples_cases"
#' @param threshold_frac_samples_per_condition: numeric: threshold for the fraction of samples in which the interactions are detected. 
#'
#' @return scatter plot, x axis is fraction of samples in controls, y axis is fraction of samples in cases
#'
#' @export
#' @examples
#' # plot the fraction of samples in which the interactions are detected
#' plot_fracSamples(df = anno_interactions_allSamples[[1]]
#'                       ,threshold_frac_samples_per_condition = 0.8
#' )
####### is this returning a list of ggplots?? because it's not clear from the code below


plot_fracSamples <- function(df,threshold_frac_samples_per_condition,...){
        # check which interactions passed the threshold
        df$passed_threshold <- (df$frac_samples_controls > threshold_frac_samples_per_condition) | (df$frac_samples_cases > threshold_frac_samples_per_condition)
        
        # main plot
        QC_plot <- ggplot(data = df
                          ,aes(x = frac_samples_controls
                               ,y = frac_samples_cases
                               ,color = passed_threshold
                               ,shape = passed_threshold
                          )
        )+
                
                geom_jitter(shape = 1
                            ,width = 0.01
                            ,height = 0.01)+
                scale_color_manual(labels = c("FALSE"="failed both\n thresholds"
                                              , "TRUE"="passed at least\none threshold")
                                   ,values=c("FALSE"="gray80"
                                             ,"TRUE"="black")
                                   ,guide = guide_legend()
                )+
                ggtitle("fraction of samples in which\nan interaction is detected")+
                xlab("controls")+
                ylab("cases")+
                geom_vline(xintercept = threshold_frac_samples_per_condition
                           , color = "red")+
                geom_hline(yintercept = threshold_frac_samples_per_condition
                           , color = "red")+
                theme_bw()+
                theme(text = element_text(size=20)
                      ,legend.position="bottom"
                      ,legend.title=element_blank())
        
        
        # x density plot
        xdensity <- ggplot(data = df
                           ,aes(x = frac_samples_controls))+
                geom_density()+
                geom_vline(xintercept = threshold_frac_samples_per_condition
                           , color = "red")+
                theme(text = element_text(size=0))+
                ggtitle("")+
                xlab("")+
                theme_bw()+
                theme(plot.margin = unit(c(0, 0, 1.8, 0), "cm"))
        
        # y density plot
        ydensity <- ggplot(data = df
                           ,aes(x = frac_samples_cases)
        )+
                geom_density()+
                xlab("")+
                geom_vline(xintercept = threshold_frac_samples_per_condition
                           , color = "red")+
                theme(text = element_text(size=0))+
                ggtitle("")+
                coord_flip()+
                theme_bw()+
                theme(plot.margin = unit(c(0.5, 1, 2.5, 2), "cm"))
        
        # Create a blank placeholder plot :
        blankPlot <- ggplot(data = df
                            ,aes(x = frac_samples_controls
                                 ,y = frac_samples_cases
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
        list(QC_plot = QC_plot
             ,ydensity = ydensity
             ,xdensity = xdensity
             ,blankPlot = blankPlot
            )
}