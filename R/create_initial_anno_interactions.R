#' @title create_initial_anno_interactions
#' 
#' @description creates initial anno_interactions for all samples. 
#' 
#' @param counts numeric [0,Inf]: normalized expression data frame contatining ligands and receptors in the rows and cells in the columns. Rownames and column names should be defined.
#' @param anno_samples dataframe: concatenated data frame of the sample annotation from all samples (rows are sample IDs, columns annotation columns: must contains "sample_ID" and case_or_control columns).
#' @param anno_cells dataframe: concatenated data frame of the cell annotation from all samples (rows are cell IDs, columns annotation columns: must contains "cell_ID", "cell_type" and "sample_ID" column).
#' @param lrp_database dataframe: contains L and R Pairs, should contain solumns Ligand.ApprovedSymbol, Receptor.ApprovedSymbol and Pair.Name.
#'
#' @return anno_interactions_allSamples: list of dataframes. The list has a length equal to the number of samples in the dataset, and has names equal to the sample IDs:
#' * Each dataframe is annotation to the interactions in the correspoding sample.
#' * Each dataframe contains interactions IDs in the rows
#' * Each dataframe comtains following columns:
#'  interaction_ID: A unique identifier for each interaction created by concatenating sending cell type, ligand gene name, receiving cell type, and receptor gene name.
#'  ligand_gene_name: The gene name of the ligand in the interaction.
#'  receptor_gene_name: The gene name of the receptor in the interaction.
#'  sending_cell_type: The cell type that expresses the ligand.
#'  receiving_cell_type: The cell type that expresses the receptor.
#'  
#' @export
#' @examples
#'      # load data
#'     data("counts")
#'    data("anno_samples")
#'   data("anno_cells")
#' data("lrp_database")
#'      # run function
#'    anno_interactions_allSamples <- create_initial_anno_interactions(counts
#'                                                    ,anno_samples
#'                                                   ,anno_cells
#'                                                 ,lrp_database
#'   )
#'     # check output
#'   print("str(anno_interactions_allSamples)")
create_initial_anno_interactions <- function(counts 
                                             ,anno_samples
                                             ,anno_cells
                                             ,lrp_database
                                             ,verbose = FALSE
){
        if(verbose){print("create initial anno_interactions")}
        
        # define interaction_IDs
        #print("define interaction_IDs")
        interaction_IDs <- c("initial_value")
        
        for(send in unique(anno_cells[,"cell_type"])
        ){
                #print("send:")
                #  print(send)
                for(rec in unique(anno_cells[,"cell_type"])
                ){
                        #print("rec:")
                        #  print(rec)
                        new_interactions <- paste0(send 
                                                   ,":"
                                                   ,lrp_database$Ligand.ApprovedSymbol 
                                                   ,"_"
                                                   ,rec
                                                   ,":" 
                                                   ,lrp_database$Receptor.ApprovedSymbol
                        )
                        
                        interaction_IDs <- c(interaction_IDs
                                              ,new_interactions
                        )
                }
        }
        interaction_IDs <- interaction_IDs[-1]
        
        # make an empty dataframe for anno_interactions
        initial_params <- c("interaction_ID"
                            ,"ligand_gene_name"
                            ,"receptor_gene_name"
                            ,"sending_cell_type"
                            ,"receiving_cell_type"
        )
        dummy_anno_interactions <- as.data.frame(matrix(,nrow = length(interaction_IDs)
                                                        ,ncol = length(initial_params)
        )
        )
        rownames(dummy_anno_interactions) <- interaction_IDs
        colnames(dummy_anno_interactions) <- initial_params
        
        # fill in interaction_IDs
        dummy_anno_interactions$interaction_ID <- interaction_IDs
        
        # fill in ligand_gene_name
        dummy_anno_interactions$ligand_gene_name <- sapply(strsplit(dummy_anno_interactions$interaction_ID
                                                                    ,"_"
        )
        ,function(interactions) strsplit(interactions[1] 
                                         ,":"
        )[[1]]
        )[2,]
        # fill in receptor_gene_name
        dummy_anno_interactions$receptor_gene_name <- sapply(strsplit(dummy_anno_interactions$interaction_ID
                                                                      ,"_"
        )
        ,function(interactions) strsplit(interactions[2]
                                         ,":"
        )[[1]]
        )[2,]
        # fill in sending_cell_type
        dummy_anno_interactions$sending_cell_type <- sapply(strsplit(dummy_anno_interactions$interaction_ID
                                                                     ,"_"
        )
        ,function(interactions) strsplit(interactions[1]
                                         ,":"
        )[[1]]
        )[1,]
        # fill in receiving_cell_type
        dummy_anno_interactions$receiving_cell_type <- sapply(strsplit(dummy_anno_interactions$interaction_ID
                                                                       ,"_")
                                                              ,function(interactions) strsplit(interactions[2] 
                                                                                               ,":"
                                                              )[[1]]
        )[1,]
        
        # make a list of anno_interactions for all samples
        anno_interactions_allSamples <- lapply(anno_samples$sample_ID
                                               ,function(sample){
                                                       dummy_anno_interactions
                                               })
        names(anno_interactions_allSamples) <- anno_samples$sample_ID
        
        anno_interactions_allSamples
}