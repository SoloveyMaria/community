#' @title clean_lrp_database
#' 
#' @param counts dataframe: expression data frame containing ligands and receptors in the rows and cells in the columns. Rownames and column names should be defined.
#' @param lrp_database dataframe: contains L and R Pairs, should contain solumns Ligand.ApprovedSymbol, Receptor.ApprovedSymbol and Pair.Name.
#' @param verbose logical: whether to print the progress of the function or not.
#' 
#' @return Dataframe: Extract L and R matched in counts for the analysis without duplicates and NA
#' 
#' @export
#' @examples
#' counts <- data.frame(
#'       gene = c('gene1','gene2','gene3')
#'    ,cell1 = c(1,2,3)
#'  ,cell2 = c(4,5,6)
#' )
#' rownames(counts) <- counts$gene
#' counts$gene <- NULL
#' lrp_database <- read.csv('LRP_database.csv')
#' clean_lrp_database(counts
#'                  ,lrp_database
#'                ,verbose = TRUE
#' )
#' 
clean_lrp_database <- function(counts, lrp_database, verbose = FALSE) {
        if (verbose) {
                print("cleaning lrp database")
        }
        
        # substitute NAs in the database with receptor.names
        idx_rec_na <- is.na(lrp_database$Receptor.ApprovedSymbol)
        lrp_database$Receptor.ApprovedSymbol[idx_rec_na] <- lrp_database$Receptor.Name[idx_rec_na]
        
        # filter out duplicated LRPs from the database
        lrp_database$lrp <- paste0(lrp_database$Ligand.ApprovedSymbol, "_", lrp_database$Receptor.ApprovedSymbol)
        # print(dim(lrp_database))
        lrp_database <- lrp_database[!duplicated(lrp_database$lrp), ]
        # print(dim(lrp_database))
        
        # filter out LRP which contain ligands or receptors that are not present in
        # the count matrix
        if (verbose) {
                print("filtering out ligands and receptor which are not present in the count matrix")
        }
        lrp_database <- lrp_database[(lrp_database$Ligand.ApprovedSymbol %in% rownames(counts)) &
                                             (lrp_database$Receptor.ApprovedSymbol %in% rownames(counts)), ]
        
}