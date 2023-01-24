#' @title plot_forest
#'
#' @description creates a horizontal bar chart of log2FC expression values for a set of interactions. The function also allows to customize the color scale and the range of the log2FC values.
#'
#' @param my_df: a data frame containing the following columns: 'interaction_ID' and 'log2FC'. The interaction_ID column must be unique.
#' @param my_title: character: title for the plot
#' @param plot_legend: logical: should the legend be plotted (default is TRUE)
#' @param min: numeric: minimal log2FC value to be plotted (default is min of my_df$log2FC)
#' @param max: numeric: maximal log2FC value to be plotted (default is max of my_df$log2FC)
#'
#' @return ggplot object
#'
#' @export
#' @examples
#' # plot_forest
#' plot_forest(my_df
#'           ,my_title
#'          ,plot_legend = TRUE
#'         ,min
#'        ,max
#' )
#' 
plot_forest <- function(my_df, my_title, plot_legend = TRUE, min, max) {
        
        my_values <- c(min, -1.1, -0.9, 0, 0.9, 1.1, max)
        my_colors <- c("blue3", "blue3", "aliceblue", "gray90", "lavenderblush", "red3",
                       "red3")
        names(my_colors) <- my_values
        # print(max)
        idx_max <- max <= my_values
        # print(idx_max)
        if (sum(!idx_max) != length(my_values) - 1) {
                ifelse(min == max, {
                        my_colors <- my_colors[c(rep(TRUE, sum(!idx_max) + 2), rep(FALSE, sum(idx_max) -
                                                                                           2))]
                }, {
                        my_colors <- my_colors[c(rep(TRUE, sum(!idx_max) + 1), rep(FALSE, sum(idx_max) -
                                                                                           1))]
                })
                
        }
        my_length <- length(my_colors)
        
        idx_min <- min >= my_values[1:my_length]
        
        if (sum(!idx_min) != my_length - 1) {
                ifelse(min == max, {
                        my_colors <- my_colors[c(rep(FALSE, sum(idx_min) - 2), rep(TRUE, sum(!idx_min) +
                                                                                           2))]
                }, {
                        my_colors <- my_colors[c(rep(FALSE, sum(idx_min) - 1), rep(TRUE, sum(!idx_min) +
                                                                                           1))]
                })
                
        }
        
        resc_values <- rescale(as.numeric(names(my_colors)))
        
        my_p <- ggplot() + geom_bar(data = my_df, aes(x = interaction_ID, y = log2FC,
                                                      color = log2FC, fill = log2FC), stat = "identity", position = "identity") +
                scale_colour_gradientn(colours = my_colors, values = resc_values) + scale_fill_gradientn(colours = my_colors,
                                                                                                         values = resc_values) + theme_void() + theme(axis.text.x = element_text(size = 0), axis.text.y = element_blank(),
                                                                                                                                                      axis.ticks.y = element_blank(), legend.position = "bottom", plot.title = element_text(hjust = 0.5)) +
                xlab("") + ylab("") + ggtitle(my_title) + coord_flip()
        
        ifelse(plot_legend, return(my_p), return(my_p + theme(legend.position = "none")))
}