#' @title f_cellType
#' 
#' @description calculates the fraction of cells of a given cell type in a sample
#' 
#' @param anno_interactions dataframe containing interactions annotation for a given sample.
#' @param cell_type string : the cell type of interest.
#' @param anno_cells dataframe: concatenated data frame of the cell annotation from all samples (rows are cell IDs, columns annotation columns: must contain "cell_ID", "cell_type" and "sample_ID" columns).
#' @param threshold_celltype_size numeric [0,Inf]: threshold of minimum number of cells in a cell type in a sample. The default value is 4.
#' 
#' @return f numeric [0,1]: fraction of cells of the given cell type in the sample.
#' 
#' @export
#' @examples
#' anno_interactions <- data.frame(sending_cell_type = c('A', 'B', 'C', 'D', 'E')
#'                         ,receiving_cell_type = c('A', 'B', 'C', 'D', 'E')
#'                        ,sending_cell = c('A1', 'B1', 'C1', 'D1', 'E1')
#'                       ,receiving_cell = c('A2', 'B2', 'C2', 'D2', 'E2')
#'                      ,sending_gene = c('A', 'B', 'C', 'D', 'E')
#'                    ,receiving_gene = c('A', 'B', 'C', 'D', 'E')
#'                  ,sending_gene_expression = c(1, 2, 3, 4, 5)
#'               ,receiving_gene_expression = c(1, 2, 3, 4, 5)
#'            ,sending_cell_expression = c(1, 2, 3, 4, 5)
#'        ,receiving_cell_expression = c(1, 2, 3, 4, 5)
#'     ,sending_cell_type_expression = c(1, 2, 3, 4, 5)
#' ,receiving_cell_type_expression = c(1, 2, 3, 4, 5)
#' )
#' cell_type <- 'A'
#' anno_cells <- data.frame(cell_type = c('A', 'B', 'C', 'D', 'E')
#'                        ,cell = c('A1', 'B1', 'C1', 'D1', 'E1')
#' )
#' f_cellType(anno_interactions
#'           ,cell_type
#'          ,anno_cells
#'        ,threshold_celltype_size = 4
#' )
#' 
f_cellType <- function(anno_interactions, cell_type, anno_cells, threshold_celltype_size = 4) {
        
        # calculate number of cells in the sample
        nr_cells_total <- nrow(anno_cells)
        
        # calculate number of cells of this cell type in this sample
        nr_cells_celltype <- sum(anno_cells$cell_type == cell_type)
        
        # calculate fraction of the cell type in the sample check if cell type size
        # passes the threshold
        ifelse(nr_cells_total > threshold_celltype_size, f <- nr_cells_celltype/nr_cells_total,
               f <- 0)
        
        return(f)
}