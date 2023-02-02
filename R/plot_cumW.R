#' @title plot_cumW
#'
#' @description creates a histogram plot of log10 cumulative weight of interactions.
#'
#' @param df: dataframe: dataframe containing the interactions, including a column "log10_cum_weight"
#' @param threshold_log10_cum_weight: numeric: threshold for log10 cumulative weight of interactions. Interactions with a log10 cumulative weight above this threshold will be considered as passing the threshold.
#'
#' @return histogram plot
#'
#' @export
#' @examples
#' # plot the cumulative weight of the interactions
#' plot_cumW(df = anno_interactions_allSamples[[1]], threshold_log10_cum_weight = 0)
#' 
####### is this returning a list of ggplots?? because it's not clear from the code below


plot_cumW <- function(df,threshold_log10_cum_weight){
        
        # check which interactions passed the threshold
        df$passed_threshold <- df$log10_cum_weight > threshold_log10_cum_weight
        
        # main plot
        ggplot(data = df
                          ,aes(x = log10_cum_weight)
        )+
                
                geom_density()+
                ggtitle("interaction weight filter")+
                geom_vline(xintercept = threshold_log10_cum_weight
                           , color = "red")+
                theme_bw()+
                theme(text = element_text(size=20)
                      ,legend.title=element_blank())
        
}