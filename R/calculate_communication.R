#' @title calculate_communication
#' 
#' @description calculates communication between cell types based on the provided count data and annotation information.
#' 
#' @param counts numeric dataframe: normalized expression data frame containing ligands and receptors in the rows and cells in the columns. Row names and column names should be defined.
#' @param anno_samples dataframe: concatenated data frame of the sample annotation from all samples (rows are sample IDs, columns annotation columns: must contain "sample_ID" and "case_or_control" columns).
#' @param anno_cells dataframe: concatenated data frame of the cell annotation from all samples (rows are cell IDs, columns annotation columns: must contain "cell_ID", "cell_type" and "sample_ID" columns).
#' @param lrp_database dataframe: contains L and R Pairs, should contain solumns Ligand.ApprovedSymbol, Receptor.ApprovedSymbol and Pair.Name.
#' @param threshold_expr numeric: expression threshold below which the expression of the ligand (receptor) is not considered. The default value is 0.05.
#' @param threshold_nr_active_cells numeric: threshold for minimal number of active cells.The default value is 0.
#' @param threshold_celltype_size numeric: threshold of minimum number of cells in a cell type in a sample. The default value is 4.
#' @param verbose a logical value indicating whether to print progress messages (default: FALSE).
#' 
#' @return list of the following objects:
#' 
#' - per_sample_anno_interactions a list of dataframes (one dataframe per sample). The names of the dataframes are equal to sample_ID. Each dataframe contains the interaction IDs as row names and the following columns:
#' -- interaction_ID: an interaction ID is constructed in the following form: "sending cell type : ligand _ receiving cell type : receptor". E.g. in an interaction "T:TGFB1_B:CXCR4", the T cells are the sending cell type expressing the TGFB1 ligand and the B cells are the receiving cell type expressing the CXCR4 receptor.
#' -- ligand_gene_name: the gene symbol of the ligand, e.g. "TGFB1"
#' -- receptor_gene_name: the gene symbol of the receptor, e.g. "CXCR4"
#' -- sending_cell_type: the sending cell type, e.g. "T"
#' -- receiving_cell_type: the receiving cell type, e.g. "B"
#' -- f_s: numeric [0,1] the abundance (as fraction) of the sending cell type within the sample of interest
#' -- f_r: numeric [0,1] the abundance (as fraction) of the receiving cell type within the sample of interest
#' -- f_s_max: numeric [0,1] the maximum abundance (as fraction) value for the sending cell type over all samples in the dataset
#' -- f_r_max: numeric [0,1] the maximum abundance (as fraction) value for the receiving cell type over all samples in the dataset
#' -- a_s_l: numeric [0,1] the fraction of the sending cells expressing the ligand above the threshold in the sample of interest
#' -- nr_s_l_active: numeric [0,Inf] the number of the sending cells expressing the ligand above the threshold in the sample of interest
#' -- a_r_r: numeric [0,1] the fraction of the receiving cells expressing the receptor above the threshold in the sample of interest
#' -- nr_r_r_active: numeric [0,Inf] the number of the receiving cells expressing the receptor above the threshold in the sample of interest
#' -- a_s_l_max: numeric [0,1] the maximum fraction of the sending cells expressing the ligand above the threshold over all samples
#' -- a_r_r_max: numeric [0,1] the maximum fraction of the receiving cells expressing the receptor above the threshold over all samples
#' -- e_s_l: numeric [0,Inf] the mean expression of the ligand in the active fraction of the sending cells in the sample of interest
#' -- e_r_r: numeric [0,Inf] the mean expression of the receptor in the active fraction of the receiving cells in the sample of interest
#' -- e_s_l_max: numeric [0,Inf] the maximum mean expression of the ligand in the active fraction of the sending cells over all samples
#' -- e_r_r_max: numeric [0,Inf] the maximum mean expression of the receptor in the active fraction of the receiving cells over all samples
#' -- rho_s: numeric [0,1] the normalized abundance of the sending cell type within the sample of interest
#' -- rho_r: numeric [0,1] the normalized abundance of the receiving cell type within the sample of interest
#' -- rho: numeric [0,1] the normalized abundance component of the weight formula
#' -- phi_s_l: numeric [0,1] the normalized fraction of the sending cells expressing the ligand above the threshold in the sample of interest
#' -- phi_r_r: numeric [0,1] the normalized fraction of the receiving cells expressing the receptor above the threshold in the sample of interest
#' -- phi: numeric [0,1] the normalized active fraction component of the weight formula
#' -- p_s_l: numeric [0,1] the normalized mean expression of the ligand in the active fraction of the sending cells in the sample of interest
#' -- p_r_r: numeric [0,1] the normalized mean expression of the receptor in the active fraction of the receiving cells in the sample of interest
#' -- p: numeric [0,1] the normalized mean expression in the active fraction component of the weight formula
#' -- w: numeric [0,1] the interaction weight
#' 
#' - f_s numeric dataframe [0,1] with row names being interaction IDs, column names being sample IDs and values being abundance (as fraction) of the sending cell type
#' 
#' - f_r numeric dataframe [0,1] with row names being interaction IDs, column names being sample IDs and values being abundance (as fraction) of the receiving cell type
#' 
#' - rho_s numeric dataframe [0,1] with row names being interaction IDs, column names being sample IDs and values being normalized abundance of the sending cell type
#' 
#' - rho_r numeric dataframe [0,1] with row names being interaction IDs, column names being sample IDs and values being normalized abundance of the receiving cell type
#' 
#' - rho numeric dataframe [0,1] with row names being interaction IDs, column names being sample IDs and values being the normalized abundance component of the weight formula
#' 
#' - a_s_l numeric dataframe [0,1] with row names being interaction IDs, column names being sample IDs and values being the fraction of the sending cells expressing the ligand above the threshold
#' 
#' - a_r_r numeric dataframe [0,1] with row names being interaction IDs, column names being sample IDs and values being the fraction of the receiving cells expressing the receptor above the threshold
#' 
#' - nr_s_l_active numeric dataframe [0,Inf] with row names being interaction IDs, column names being sample IDs and values being the number of the sending cells expressing the ligand above the threshold
#' 
#' - nr_r_r_active numeric dataframe [0,Inf] with row names being interaction IDs, column names being sample IDs and values being the number of the receiving cells expressing the receptor above the threshold
#' 
#' - phi_s_l numeric dataframe [0,1] with row names being interaction IDs, column names being sample IDs and values being the normalized fraction of the sending cells expressing the ligand above the threshold
#' 
#' - phi_r_r numeric dataframe [0,1] with row names being interaction IDs, column names being sample IDs and values being the normalized fraction of the receiving cells expressing the receptor above the threshold
#' 
#' - phi numeric dataframe [0,1] with row names being interaction IDs, column names being sample IDs and values being the normalized fraction of active cells component of the weight formula
#' 
#' - e_s_l numeric dataframe [0,Inf] with row names being interaction IDs, column names being sample IDs and values being the mean expression of the ligand in the active fraction of the sending cells
#' 
#' - e_r_r numeric dataframe [0,Inf] with row names being interaction IDs, column names being sample IDs and values being the mean expression of the receptor in the active fraction of the receiving cells
#' 
#' - p_s_l numeric dataframe [0,1] with row names being interaction IDs, column names being sample IDs and values being the normalized mean expression of the ligand in the active fraction of the sending cells
#' 
#' - p_r_r numeric dataframe [0,1] with row names being interaction IDs, column names being sample IDs and values being the normalized mean expression of the receptor in the active fraction of the receiving cells
#' 
#' - p numeric dataframe [0,1] with row names being interaction IDs, column names being sample IDs and values being the normalized mean expression in the active fraction component of the weight formula
#' 
#' - weights numeric dataframe [0,1] with row names being interaction IDs, column names being sample IDs and values being the interaction weight
#' 
#' - thresholds: a list of:
#' -- threshold_expr numeric: expression threshold below which the expression of the ligand (receptor) is not considered
#' -- threshold_nr_active_cells: numeric: threshold for minimal number of active cells
#' -- threshold_celltype_size: numeric: threshold of minimum number of cells in a cell type in a sample
#' 
#' - anno_interactions: a list of:
#' -- interaction_ID: an interaction ID is constructed in the following form: "sending cell type : ligand _ receiving cell type : receptor". E.g. in an interaction "T:TGFB1_B:CXCR4", the T cells are the sending cell type expressing the TGFB1 ligand and the B cells are the receiving cell type expressing the CXCR4 receptor.
#' -- ligand_gene_name: the gene symbol of the ligand, e.g. "TGFB1"
#' -- receptor_gene_name: the gene symbol of the receptor, e.g. "CXCR4"
#' -- sending_cell_type: the sending cell type, e.g. "T"
#' -- receiving_cell_type: the receiving cell type, e.g. "B"
#' 
#' - anno_samples dataframe: concatenated data frame of the sample annotation from all samples (rows are sample IDs, columns annotation columns: must contain "sample_ID" and "case_or_control" columns).
#' 
#' - anno_cells dataframe: concatenated data frame of the cell annotation from all samples (rows are cell IDs, columns annotation columns: must contain "cell_ID", "cell_type" and "sample_ID" columns).
#' 
#' @export
#' @examples
#' counts <- data.frame(
#'       gene = c("gene1","gene2","gene3")
#'    ,cell1 = c(1,2,3)
#'  ,cell2 = c(4,5,6)
#' )
#' rownames(counts) <- counts$gene
#' counts$gene <- NULL
#' anno_cells <- data.frame(
#'      cell_ID = c("cell1","cell2")
#'   ,cell_type = c("cell_type1","cell_type2")
#' )
#' anno_samples <- data.frame(
#'      sample_ID = c("sample1")
#'  ,health_status = c("healthy")
#' ,case_or_control = c("case")
#' )
#' lrp_database <- read.csv("LRP_database.csv")
#' calculate_communication(counts = counts
#'                       ,anno_samples = anno_samples
#'                      ,anno_cells = anno_cells
#'                     ,lrp_database = lrp_database
#'                   ,verbose = TRUE
#' )
#' 
calculate_communication <- function(counts
                                    ,anno_samples
                                    ,anno_cells
                                    ,threshold_expr = 0.05
                                    ,threshold_nr_active_cells = 0
                                    ,threshold_celltype_size = 4
                                    ,lrp_database
                                    ,verbose = FALSE){
        
        # clean LRP database
        lrp_database_clean <- clean_lrp_database(counts = counts
                                                 ,lrp_database = lrp_database
                                                 ,verbose = verbose)
        
        # create annotation for the interactions for each sample
        per_sample_anno_interactions <- create_initial_anno_interactions(anno_samples = anno_samples
                                                                         ,anno_cells = anno_cells
                                                                         ,lrp_database = lrp_database_clean
                                                                         ,verbose = verbose
        )
        #print("str(anno_interactions[[1]])")
        #print(str(anno_interactions[[1]]))
        
        # Before calculating individual weights for each sample, we will calculate following parameters for all samples:
        # - f_cellType: fractions for each cell type in each samples
        # - f_cellType_max: for each cell type, max fraction over all samples
        # - TODO!!!!!!! active fraction for each cell type, each gene in each sample
        # - TODO!!!!!!! for each cell type and each gene, max active fraction over all samples
        # - e_cellType_gene: mean expression in the active fraction for each cell type, each gene in each sample
        # - e_cellType_gene_max: for each cell type and each gene, max mean expression in the active fraction over all samples
        
        # calculate f_cellTypes: 
        # - f_s
        # - f_r
        per_sample_anno_interactions <- f_ct_allSamples(anno_cells = anno_cells
                                                        ,anno_interactions_allSamples = per_sample_anno_interactions
                                                        ,threshold_celltype_size = threshold_celltype_size
                                                        ,verbose = verbose
        )
        
        # calculate f_cellTypes_max:
        # - f_s_max
        # - f_r_max
        per_sample_anno_interactions <- f_cellType_max(anno_cells = anno_cells
                                                       ,anno_interactions_allSamples = per_sample_anno_interactions
                                                       ,threshold_celltype_size = threshold_celltype_size
                                                       ,verbose = verbose
        )
        
        # calculate following parameters for the anno_interactions for each sample:
        # - a_l_s
        # - a_r_r
        # - nr_l_s_active
        # - nr_l_s_active
        per_sample_anno_interactions <- a_ct_g_allSamples(counts = counts
                                                          ,anno_cells = anno_cells
                                                          ,anno_interactions_allSamples = per_sample_anno_interactions
                                                          ,threshold_expr = threshold_expr
                                                          ,threshold_nr_active_cells = threshold_nr_active_cells
                                                          ,verbose = verbose
        )
        
        # calculate a_ct_g_max:
        # - a_s_l_max
        # - a_r_r_max
        per_sample_anno_interactions <- a_cellType_gene_max(anno_interactions_allSamples = per_sample_anno_interactions
                                                            ,verbose = verbose
        )
        
        # calculate following parameters for the anno_interactions for each sample:
        # - e_s_l
        # - e_r_r
        per_sample_anno_interactions <- e_ct_g_allSamples(counts = counts
                                                          ,anno_cells = anno_cells
                                                          ,anno_interactions_allSamples = per_sample_anno_interactions
                                                          ,threshold_expr = threshold_expr
                                                          ,verbose = verbose
        )
        
        # calculate following parmeters for the anno_interactions for each sample:
        # - e_s_l_max
        # - e_r_r_max
        per_sample_anno_interactions <- e_cellType_gene_max(anno_interactions_allSamples = per_sample_anno_interactions
                                                            ,verbose = verbose
        )
        
        
        
        # Now for each samples, individual interaction weights will be calculated
        
        # calculatefollowing parameters for the anno_interactions for each sample:
        # - rho_s -> calculate
        # - rho_r -> calculate
        # - rho -> calculate
        # - phi_s_l -> calculate
        # - phi_r_r -> calculate
        # - phi -> calculate
        # - p_l_s -> calculate
        # - p_r_r -> calculate
        # - p -> calculate
        # - w -> calculate
        samples <- names(per_sample_anno_interactions)
        per_sample_anno_interactions <- lapply(samples
                                               ,function(sample){
                                                       if(verbose)print(paste("for sample",sample,":"))
                                                       
                                                       calculate_weight(anno_interactions = per_sample_anno_interactions[[sample]]
                                                                        ,verbose = verbose
                                                       )
                                               })
        names(per_sample_anno_interactions) <- samples
        
        
        # extract f_s
        f_s <- extract_params_annointeractions(anno_interactions = per_sample_anno_interactions
                                               ,param = "f_s")
        # extract f_r
        f_r <- extract_params_annointeractions(anno_interactions = per_sample_anno_interactions
                                               ,param = "f_r")
        # extract rho_s
        rho_s <- extract_params_annointeractions(anno_interactions = per_sample_anno_interactions
                                                 ,param = "rho_s")
        # extract rho_r
        rho_r <- extract_params_annointeractions(anno_interactions = per_sample_anno_interactions
                                                 ,param = "rho_r")
        # extract rho
        rho <- extract_params_annointeractions(anno_interactions = per_sample_anno_interactions
                                               ,param = "rho")
        # extract a_s_l
        a_s_l <- extract_params_annointeractions(anno_interactions = per_sample_anno_interactions
                                                 ,param = "a_s_l")
        # extract a_r_r
        a_r_r <- extract_params_annointeractions(anno_interactions = per_sample_anno_interactions
                                                 ,param = "a_r_r")
        # extract nr_l_s_active
        nr_s_l_active <- extract_params_annointeractions(anno_interactions = per_sample_anno_interactions
                                                         ,param = "nr_s_l_active")
        # extract nr_r_r_active
        nr_r_r_active <- extract_params_annointeractions(anno_interactions = per_sample_anno_interactions
                                                         ,param = "nr_r_r_active")
        # extract phi_l_s
        phi_s_l <- extract_params_annointeractions(anno_interactions = per_sample_anno_interactions
                                                   ,param = "phi_s_l")
        # extract phi_r_r
        phi_r_r <- extract_params_annointeractions(anno_interactions = per_sample_anno_interactions
                                                   ,param = "phi_r_r")
        # extract phi
        phi <- extract_params_annointeractions(anno_interactions = per_sample_anno_interactions
                                               ,param = "phi")
        # extract e_s_l
        e_s_l <- extract_params_annointeractions(anno_interactions = per_sample_anno_interactions
                                                 ,param = "e_s_l")
        # extract e_r_r
        e_r_r <- extract_params_annointeractions(anno_interactions = per_sample_anno_interactions
                                                 ,param = "e_r_r")
        # extract p_l_s
        p_s_l <- extract_params_annointeractions(anno_interactions = per_sample_anno_interactions
                                                 ,param = "p_s_l")
        # extract phi_r_r
        p_r_r <- extract_params_annointeractions(anno_interactions = per_sample_anno_interactions
                                                 ,param = "p_r_r")
        # extract p
        p <- extract_params_annointeractions(anno_interactions = per_sample_anno_interactions
                                             ,param = "p")
        # extract interactions weights
        weights  <-  extract_params_annointeractions(anno_interactions = per_sample_anno_interactions
                                                     ,param = "w")
        
        # create anno
        anno <- data.frame(interaction_ID = rownames(weights)
                           ,ligand_gene_name = per_sample_anno_interactions[[1]]$ligand_gene_name
                           ,receptor_gene_name = per_sample_anno_interactions[[1]]$receptor_gene_name
                           ,sending_cell_type = per_sample_anno_interactions[[1]]$sending_cell_type
                           ,receiving_cell_type = per_sample_anno_interactions[[1]]$receiving_cell_type
        )
        
        # create threshold variable
        thresholds <- list(threshold_expr = threshold_expr
                           ,threshold_nr_active_cells = threshold_nr_active_cells
                           ,threshold_celltype_size = threshold_celltype_size
        )
        
        
        return(list(per_sample_anno_interactions = per_sample_anno_interactions
                    ,f_s = f_s
                    ,f_r = f_r
                    ,rho_s = rho_s
                    ,rho_r = rho_r
                    ,rho = rho
                    ,a_s_l = a_s_l
                    ,a_r_r = a_r_r
                    ,nr_s_l_active = nr_s_l_active
                    ,nr_r_r_active = nr_r_r_active
                    ,phi_s_l = phi_s_l
                    ,phi_r_r = phi_r_r
                    ,phi = phi
                    ,e_s_l = e_s_l
                    ,e_r_r = e_r_r
                    ,p_s_l = p_s_l
                    ,p_r_r = p_r_r
                    ,p = p
                    ,weights = weights
                    ,thresholds = thresholds
                    ,anno_interactions = anno
                    ,anno_samples = anno_samples
                    ,anno_cells = anno_cells
        )
        )
}