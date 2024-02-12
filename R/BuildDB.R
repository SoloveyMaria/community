#' Import Ligand-Receptor Interaction Databases
#'
#' This function imports ligand-receptor interaction data, with options to select from non-curated, curated, or both types of databases. It processes the data to avoid duplicates and provides unified formatting.
#'
#' @param db_type A character vector specifying the type of database to import.
#'                Options are "noncurated", "curated", or "both". Default is "both".
#'                The function uses this parameter to determine which database(s) to load
#'                and how to process them.
#'
#' @return Depending on the selected database type, this function returns either a dataframe
#'         of non-curated interactions, curated interactions, or a combined dataframe of both.
#'         If "both" is selected, it also provides a count of interaction pairs by their annotation strategy.
#' @import dplyr
#' @import tidyverse
#'
#' @examples
#' # Import only non-curated data
#' non_curated_data <- import_db("noncurated")
#'
#' # Import only curated data
#' #curated_data <- import_db("curated")
#'
#' # Import both curated and non-curated data
#' # both_data <- import_db("both")
#'
#' @export

import_db <- function(db_type = c("noncurated", "curated", "both")) {
    db_type <- match.arg(db_type)

    if (db_type %in% c("noncurated", "both")) {
        non_curated <- import_ligrecextra_interactions()
        non_curated <- non_curated %>% filter(!duplicated(non_curated[, c("source_genesymbol", "target_genesymbol")]))
        non_curated$Pair.Name <- paste(non_curated$source_genesymbol, non_curated$target_genesymbol, sep = "_")
        non_curated$annotation_strategy <- "LR"
  }

    if (db_type %in% c("curated", "both")) {
        curated <- curated_ligand_receptor_interactions()
        curated <- curated %>% filter(!duplicated(curated[, c("source_genesymbol", "target_genesymbol")]))
        curated$Pair.Name <- paste(curated$source_genesymbol, curated$target_genesymbol, sep = "_")
        curated$annotation_strategy <- "curated"
  }

    if (db_type == "both") {
    non_curated <- non_curated %>%
      mutate(annotation_strategy = ifelse(Pair.Name %in% curated$Pair.Name, "both", annotation_strategy))

    combined_db <- rbind(non_curated, curated)
#     combined_db <- combined_db[!duplicated(combined_db$Pair.Name), ]

    print(paste0("Retrieved interactions from ", db_type, " DB"))
    return(combined_db)
    } else if (db_type == "noncurated") {
    return(non_curated)
    } else if (db_type == "curated") {
    return(curated)
    }
}


#' Create Pairwise Pairs from Complex Interactions
#'
#' This function takes OmniPathDB, processes complex ligand-receptor interactions and
#' creates pairwise combinations. It is specifically designed to handle cases where
#' either the source or the target is a complex. The function filters for these complex rows, splits the gene
#' symbols, and creates all possible pairwise interactions.
#'
#' @param both_db A data frame containing the ligand-receptor interaction data.
#'                It is expected to have columns `source_genesymbol` and
#'                `target_genesymbol`, among others. Complex interactions are
#'                identified by the presence of "COMPLEX" in either source or target fields.
#'
#' @return A data frame containing the expanded pairwise interactions. Each row
#'         represents a unique ligand-receptor pair, along with the original complex
#'         information and other relevant data from the input data frame.
#' @importFrom stringr str_detect str_split
#' @importFrom utils combn
#' @export
#'
#' @examples
#' # Assuming `both_db` is your data frame with complex interactions
#' non_curated_data <- import_db("noncurated")
#' pairwise_interactions <- create_pairwise_pairs(non_curated_data)
#' head(pairwise_interactions)


create_pairwise_pairs <- function(both_db) {

    # Filter for complex rows
    complex <- both_db %>% 
               filter(str_detect(target, "COMPLEX") | str_detect(source, "COMPLEX"))
    print(paste0(nrow(complex), " Number of complex pairs detected"))
    # Remove pair column if exists
    complex$Pair.Name <- NULL

    # Initialize a list to store results
    results_list <- list()

    # Process each row
    for (i in 1:nrow(complex)) {
        values1 <- str_split(complex[i, "source_genesymbol"], "_", simplify = TRUE)
        values2 <- str_split(complex[i, "target_genesymbol"], "_", simplify = TRUE)

        original <- paste(complex[i, "source_genesymbol"], complex[i, "target_genesymbol"], sep = "_")
        
        pairs <- combn(c(values1, values2), 2)
        pairs <- t(pairs)
        pairs <- rbind(pairs, pairs[, 2:1])
        
        pairs_df <- as.data.frame(pairs)
        names(pairs_df) <- c("Ligand", "Receptor")
        pairs_df$complex_pair <- original
        
        pairs_df <- cbind(pairs_df, complex[i, ])
        results_list[[i]] <- pairs_df
    }

    # Combine all results
    results <- do.call(rbind, results_list)
    results$Pair.Name <- paste(results$Ligand, results$Receptor, sep = "_")
    print(paste0(nrow(results), " Number of non-redundant binary pairs produced"))
    return(results[, c("Pair.Name", names(results)[!names(results) %in% "Pair.Name"])])
}



#' Filter Pairwise Pairs Based on PPI Network
#'
#' Filters the pairwise pairs of ligand-receptor interactions by checking against PPI network.
#' 
#' @param pairwise_pairs A dataframe containing pairwise pairs of ligand-receptor interactions.
#' 
#' @return A dataframe of filtered pairwise pairs that exist in the PPI network.
#' 
#' @export
#' 
#' @examples
#' both_db <- import_db("both")
#' pairwise_pairs <- create_pairwise_pairs(both_db)
#' filtered_pairs <- filter_pairs_with_ppi(pairwise_pairs)

# Function to filter pairwise pairs based on PPI network
# Function to filter pairwise pairs based on PPI network
filter_pairs_with_ppi <- function(pairwise_pairs) {

    # Import all PPI
    pt <- import_post_translational_interactions()
    ppi_network <- pt %>% filter(!duplicated(.[, c("source_genesymbol", "target_genesymbol")]))
    ppi_network$Pair.Name <- paste(ppi_network$source_genesymbol, ppi_network$target_genesymbol, sep = "_")

    # Filter pairs that exist in the PPI network
    pt_interactions <- pairwise_pairs %>%
        filter(Pair.Name %in% ppi_network$Pair.Name) %>%
        distinct(Pair.Name, .keep_all = TRUE)
    print(paste0(nrow(pt_interactions), " Number of binary pairs detected through PPI"))
    return(pt_interactions)
}


#' Process Binary Pairs
#'
#' This function processes the single component interactions from a combined ligand-receptor database and merges them with PPI interactions.
#'
#' @param both_db A data frame of combined ligand-receptor interactions.
#' @param pt_interactions A data frame of post-translational interactions.
#'
#' @return A data frame of merged binary pairs and interactions filtered through PPI, without duplicates.
#'
#' @usage
#' process_binary_pairs(both_db, pt_interactions)
#'
#' @examples
#' both_db <- import_db("both")
#' pairwise_pairs <- create_pairwise_pairs(both_db)
#' pt_interactions <- filter_pairs_with_ppi(pairwise_pairs)
#' # Assuming 'both_db' and 'pt_interactions' are available data frames
#' complete_interactions <- process_binary_pairs(both_db, pt_interactions)
#'
#' @export


process_binary_pairs <- function(both_db, pt_interactions) {
    # Filter out single components
    single_components <- filter(both_db, !grepl('COMPLEX', target) & !grepl('COMPLEX', source))
    single_components$pair <- NULL  # Remove pair column

    # Process and rename columns
    single_components <- single_components %>%
        dplyr::rename(Ligand = source_genesymbol, Receptor = target_genesymbol) %>%
        dplyr::mutate(complex_pair = NA, 
                      Pair.Name = paste(Ligand, Receptor, sep = "_"))

    # Reorder columns
    single_components <- data.frame(Pair.Name = single_components$Pair.Name, 
                                    single_components[, !(names(single_components) %in% "Pair.Name")])

    # Clean pt_interactions data
    pt_interactions$target_genesymbol <- NULL
    pt_interactions$source_genesymbol <- NULL
    
    desired_cols <- c('Pair.Name','Ligand','Receptor','source','target','is_directed','is_stimulation',
                  'is_inhibition','consensus_direction','consensus_stimulation','consensus_inhibition',
                  'sources','references','curation_effort','n_references','n_resources','annotation_strategy',
                  'complex_pair')

    
    single_components <- single_components[desired_cols]
    pt_interactions <- pt_interactions[desired_cols]
    
    # Merge single components with PT interactions and drop duplicates
    complete <- rbind(single_components, pt_interactions)
    complete <- complete[!duplicated(complete$Pair.Name, fromLast = TRUE),]
    print(paste0(nrow(complete), " Non-redundant number of pairs in the DB"))
    return(complete)
}




#' Map Gene Symbols to Descriptions in Ligand-Receptor Interactions
#'
#' This function maps gene symbols to their corresponding descriptions in a complete dataset of ligand-receptor interactions. It utilizes mygene library for protein descriptions.
#'
#' @param complete A data frame containing ligand-receptor interactions with columns 'Ligand' and 'Receptor'.
#'
#' @return The input data frame with additional columns 'Ligand.Name' and 'Receptor.Name' containing the mapped protein descriptions. The function also adjusts specific gene symbols and reorders the columns.
#'
#' @usage
#' map_gene_data(complete)
#'
#' @examples
#' # Assuming 'complete' is a data frame with ligand-receptor interactions
#' complete_mapped <- map_gene_data(complete)
#'
#' @export
#' @importFrom stringr str_replace


map_gene_data <- function(complete) {
    # Get unique gene symbols
    warning("If this function fails, it may be due to internet connectivity issues. Try running it again.")
    
    gene_symbols <- unique(c(complete$Ligand, complete$Receptor))

    # Query for protein descriptions
    prot_descriptions <- queryMany(gene_symbols, scopes = "symbol", 
                                   fields = c("name"), 
                                   species = "human", 
                                   as_dataframe = TRUE)
    prot_descriptions <- as.data.frame(prot_descriptions)

    # Map protein descriptions to the complete dataset
    for (x in 1:nrow(complete)) {
        ligand_symbol = complete[x,]$Ligand
        receptor_symbol = complete[x,]$Receptor
        ligand_description = filter(prot_descriptions, query == ligand_symbol)$name
        receptor_description = filter(prot_descriptions, query == receptor_symbol)$name

        complete[x, "Ligand.Name"] = ligand_description[1]
        complete[x, "Receptor.Name"] = receptor_description[1]
    }

    # Handle specific case for "PIK3CD-AS1"
    complete$Ligand <- str_replace(complete$Ligand, "PIK3CD-AS1", "PIK3CD")
    complete$Pair.Name <- paste(complete$Ligand, complete$Receptor, sep = "_")
    complete$dup <- paste(complete$Receptor, complete$Ligand, sep = "_")

    # Reorder columns
    desired_order <- c("Pair.Name", "Ligand", "Ligand.Name", "Receptor", "Receptor.Name", "complex_pair")
    remaining_cols <- setdiff(names(complete), desired_order)
    final_order <- c(desired_order, remaining_cols)

    # Reorder the dataframe columns
    complete <- complete[, final_order]

    return(complete)
}

# complete <- map_and_organize_gene_data(complete)



#' Annotate Components 
#'
#' This function annotates components (genes) from the PPI network with their parent category. It uses data from OmniPath to determine if the components are categorized as ligands or receptors, and handles cases where the components fall into other categories.
#'
#' @param complete_data A data frame with columns 'Ligand' and 'Receptor'.
#'
#' @return A data frame with columns 'genesymbol', 'score', and 'parent', where 'parent' indicates the parent category of the gene
#'
#' @importFrom dplyr filter
#'
#' @usage
#' annotate_components(complete_data)
#'
#' @examples
#' mock <- data.frame(Ligand = c("IL10","IGF1"), Receptor = c("JAK2","ITGA8"))
#' annotated_data <- annotate_components(mock)
#'
#' @export



# This function is to annotate the components from the PPI network with their parent category
# The input of this function is a df with columns of "genesymbol", "score", "parent"
# The output of this function is a df with columns of "genesymbol", "score", "parent"
# This function is used to annotate the components from the PPI network with their parent category

annotate_components <- function(complete_data) {
    
    components <- unique(c(complete_data$Ligand, complete_data$Receptor))
    
    #create a df to store annotation
    df <- data.frame(genesymbol = character(length(components)), score = numeric(length(components)),
                     parent = character(length(components)), stringsAsFactors = FALSE)
    
    
    anno_raw <- import_omnipath_intercell()
    anno_lig <- anno_raw %>%
    dplyr::filter(category %in% c("receptor","ligand"))
    
    
    # Check if the components are categorized as ligands or receptors
    for (x in 1:length(components)) {
    #     maxvalue=max(filter(anno, uniprot==components[x])$consensus_score)
        genename <- components[x]
        parent_score <- sort(table(filter(anno_lig, genesymbol==components[x])$parent), decreasing = T, na.last = T)[1]
        parent_category <- names(parent_score)

        if (is.null(parent_category)) {
          parent_category <- "NA"
          parent_score <- 0
        }

        df[x, "genesymbol"] <- genename
        df[x, "score"] <- parent_score
        df[x, "parent"] <- parent_category

    #     df$genesymbol[x] <- genename
    #     df$score[x] <- parent_score
    #     df$parent[x] <- parent_category
    }
    
    # If a component is not classified as a ligand or receptor, we may consider other categories such as 
    # extracellular matrix, secreted, and transmembrane.# annotated others such as secreted, ecm etc

    df_na <- filter(df, parent=="NA")$genesymbol

    if (length(df_na) > 0) {
        for (x in 1:length(df_na)) {
            genesymbol <- df_na[x]
            parent_score <- sort(table(filter(anno_raw, genesymbol == df_na[x])$parent), decreasing = TRUE, na.last = TRUE)[1]
            parent_category <- names(parent_score)

            df <- df %>% mutate(parent = ifelse(genesymbol == df_na[x], parent_category, parent))
            df <- df %>% mutate(score = ifelse(genesymbol == df_na[x], parent_score, score))
        }
    }

    # replace ecm and secreted with ligand
    df$parent <- replace(df$parent, df$parent == "ecm", "ligand")
    df$parent <- replace(df$parent, df$parent == "secreted", "ligand")
    
    return(df)

}


#' Process Ligand-Receptor Database
#'
#' This function processes a ligand-receptor database, ensuring correct ligand-receptor direction based on provided annotation.
#'
#' @param complete Database.
#' @param annotation A data frame with gene annotations.
#' @param all_ligands (optional) A vector of all ligand genesymbols.
#' @param all_receptors (optional) A vector of all receptor genesymbols.
#'
#' @return A data frame of the processed ligand-receptor pairs with correct ligand-receptor directions.
#'
#' @usage 
#' process_lr_db(complete, annotation, all_ligands, all_receptors)
#'
#' @examples
#' # Assuming 'complete', 'annotation', 'all_ligands', and 'all_receptors' are available
#' processed_db <- process_lr_db(complete, annotation, all_ligands, all_receptors)
#'
#' @export



process_lr_db <- function(complete, annotation, all_ligands = NULL, all_receptors = NULL) {

    # If not provided, determine the true LR genespace from the annotation
    if (is.null(all_ligands) || is.null(all_receptors)) {
        true_LR_anno <- filter(annotation, parent == "receptor" | parent == "ligand")
        all_ligands <- filter(true_LR_anno, parent == "ligand")$genesymbol
        all_receptors <- filter(true_LR_anno, parent == "receptor")$genesymbol
    }

    # Filter pairs that need direction fixing
    LR_fix_dir <- filter(complete, Ligand %in% all_receptors & Receptor %in% all_ligands)

    # Swap values in Ligand and Receptor columns
    temp <- LR_fix_dir$Ligand
    LR_fix_dir$Ligand <- LR_fix_dir$Receptor
    LR_fix_dir$Receptor <- temp

    # Swap values in Ligand.Name and Receptor.Name columns
    temp <- LR_fix_dir$Ligand.Name
    LR_fix_dir$Ligand.Name <- LR_fix_dir$Receptor.Name
    LR_fix_dir$Receptor.Name <- temp

    # Swap values in source and target columns
    temp <- LR_fix_dir$source
    LR_fix_dir$source <- LR_fix_dir$target
    LR_fix_dir$target <- temp

    rm(temp)
    
    LR_fix_dir["Pair.Name"] <- paste(LR_fix_dir$Ligand, LR_fix_dir$Receptor, sep="_")
    LR_fix_dir["dup"] <- paste(LR_fix_dir$Receptor, LR_fix_dir$Ligand, sep="_")

    # Get the list of interactions that are strictly in LR direction
    true_LR_DB <- filter(complete, Ligand %in% all_ligands & Receptor %in% all_receptors)

    # Remove duplicates after fix
    LR_fix_dir <- LR_fix_dir[!LR_fix_dir$Pair.Name %in% true_LR_DB$Pair.Name,]

    # Add True_LR column
    true_LR_DB["True_LR"] <- TRUE
    LR_fix_dir["True_LR"] <- TRUE

    true_LR_DB <- rbind(true_LR_DB, LR_fix_dir)

    #move column to the first
    true_LR_DB <- true_LR_DB %>% dplyr::select(True_LR, everything())

    true_LR_DB$dup <- NULL

    # Combine and return the dataset
    return(true_LR_DB)
}

# true_LR_DB <- fix_lr_directions(complete_data, annotation_data)






#' Process Adhesive Ligand-Receptor Database
#'
#' This function processes an adhesive database. It involves filtering, rearranging, and annotating the database based on gene families and interaction properties. This includes handling of directional and non-directional interactions, as well as manual annotations for specific gene families.
#'
#' @param complete A data frame representing the complete set of ligand-receptor interactions.
#' @param annotation A data frame containing annotations for ligands and receptors.
#' @param ligand_list A list of gene symbols representing ligands. Optional, used for additional annotations.
#' @param receptor_list A list of gene symbols representing receptors. Optional, used for additional annotations.
#'
#' @return A processed data frame of the adhesive ligand-receptor database, annotated and filtered based on various criteria including gene family and interaction directionality.
#'
#' @usage 
#' process_adhesive_DB(complete, annotation, ligand_list = list(), receptor_list = list())
#'
#' @examples
#' # Assuming 'complete' and 'annotation' are previously prepared data frames
#' adhesive_db <- process_adhesive_DB(complete, annotation)
#'
#' # With additional ligand and receptor lists
#' adhesive_db <- process_adhesive_DB(complete, annotation, ligand_list, receptor_list)
#'
#' @export


process_adhesive_DB <- function(complete, annotation, ligand_list=list(), receptor_list=list()) {
    
    true_LR_anno <- filter(annotation, parent == "receptor" | parent == "ligand")
    anno_ligands <- filter(true_LR_anno, parent == "ligand")$genesymbol
    anno_receptors <- filter(true_LR_anno, parent == "receptor")$genesymbol
    
    LR_DB <- filter(complete, Ligand %in% anno_ligands & Receptor %in% anno_receptors)
    LR_fixed <- filter(complete, Ligand %in% anno_receptors & Receptor %in% anno_ligands)
    LR_DB <- rbind(LR_DB, LR_fixed)
    
    # Filter out rows not in true_LR_DB
    adhesive_DB <- filter(complete, !Pair.Name %in% LR_DB$Pair.Name)
    adhesive_DB["True_LR"] <- FALSE

    # Find reversed pairs
    reversed <- adhesive_DB[adhesive_DB$dup %in% adhesive_DB$Pair.Name,]
    adhesive_DB <- adhesive_DB[!adhesive_DB$Pair.Name %in% reversed$Pair.Name,]
    
    # Define gene families and annotations
    # manual annotation of genes

    plexin_family <- as.vector(reversed[grep("plexin", reversed$Receptor.Name), ]$Receptor)

    neuroligin_family <- as.vector(reversed[grep("neuroligin", reversed$Receptor.Name), ]$Receptor)

    adam_family <- as.vector(reversed[grep("ADAM", reversed$Receptor.Name), ]$Receptor)

    #extract all_receptors that has annotation of "receptor" under Ligand.Name
    receptor_anno <- as.vector(reversed[grep("receptor", reversed$Ligand.Name), ]$Ligand)
    
    
    # Combine the additional all_ligands with the plexin, neuroligin, and ADAM families into a vector called ligand
    ligand_list <- unique(c(ligand_list,plexin_family,neuroligin_family,adam_family))
    receptor_list <- unique(c(receptor_list,receptor_anno))
    
    
    
    # processing swapped duplicated
    # Subset the data frame to only include rows where the consensus_direction column is 1
    dir <- reversed %>% filter(Pair.Name %in% reversed$dup & consensus_direction == 1)
    reversed$dup = paste(reversed$Receptor, reversed$Ligand, sep="_")
    dir = reversed[reversed$Pair.Name %in% reversed$dup & reversed$consensus_direction == 1, ]
    

    # Subset the data frame to only include rows where the consensus_direction column is 0
    no_dir <- reversed %>% filter(Pair.Name %in% reversed$dup & consensus_direction == 0)
    
    
    # Remove rows from no_dir where the pair is already present in dir
    in_dir <- dir[dir$dup %in% no_dir$Pair.Name,] 
    no_dir <- no_dir[!no_dir$dup %in% in_dir$Pair.Name,  ] #removal of those in dir below is rm of nodir

    dir <- dir[!dir$Pair.Name %in% in_dir$Pair.Name,]
    

    # remove the interactions where receptor is annotated as ligand
    wrong_lig <- no_dir[no_dir$Receptor %in% ligand_list,]
    correct_lig <- no_dir[no_dir$dup %in% wrong_lig$Pair.Name,]
    no_dir <- no_dir[!no_dir$Pair.Name %in% c(wrong_lig$Pair.Name, correct_lig$Pair.Name),]

    wrong_rec <- no_dir[no_dir$Ligand %in% receptor_list,]
    wrong_rec <- wrong_rec[!wrong_rec$Pair.Name %in% wrong_rec$dup,]

    # remove the interactions where receptor is annotated as ligand
    wrong_rec <- no_dir[no_dir$Ligand %in% receptor_list,]
    wrong_rec <- wrong_rec[!wrong_rec$Pair.Name %in% wrong_rec$dup,]
    correct_rec <- no_dir[no_dir$dup %in% wrong_rec$Pair.Name,]
    no_dir <- no_dir[!no_dir$Pair.Name %in% c(wrong_rec$Pair.Name, correct_rec$Pair.Name),]
    
    
    # Function to lexographically sort the gene pairs
    sort_pairs <- function(pair) {
        parts <- strsplit(pair, "_")[[1]]
        sorted_parts <- sort(parts)
        return(paste(sorted_parts, collapse = "_"))
  }

    
    no_dir <- no_dir[order( no_dir[,2], no_dir[,4] ),]
    no_dir$sort <- sapply(no_dir$Pair.Name, sort_pairs)
    no_dir <- no_dir %>% distinct(sort, .keep_all = TRUE)
    no_dir$sort <- NULL
    
    
    c# some of the swapped duplicates are both directional, for this one, we also keep the lexograph order
    dir <- dir[order( dir[,2], dir[,4] ),]
    dir$sort <- sapply(dir$Pair.Name, sort_pairs)
    dir <- dir %>% distinct(sort, .keep_all = TRUE)
    dir$sort <- NULL
    
    # Combine dataframes and return
    subset_lr <- rbind(dir, no_dir, in_dir, correct_lig, correct_rec)
    subset_lr["True_LR"] <- FALSE

    adhesive_DB <- rbind(subset_lr, adhesive_DB)


    #move column to the first
    adhesive_DB <- adhesive_DB %>% dplyr::select(True_LR, everything())
    
    adhesive_DB$dup <- NULL

    return(adhesive_DB)
}

#' Auto Update Database
#'
#' This function automatically updates the database by processing 
#' non-curated, curated, or both datasets.
#' 
#' @param noncurated Logical, if TRUE processes the non-curated database.
#' @param curated Logical, if TRUE processes the curated database.
#' @param both Logical, if TRUE processes both the non-curated and curated databases.
#' @return Returns the updasted database.
#' @export
#' @examples
#' updated_db <- auto_update_db("noncurated")

auto_update_db <- function(db_type) {
    # manual annotation of genes


    ligand_list <- c("AGRN", "BMP2", "BMP4", "VTCN1", "CD244", "CD38", "GAS6", "GDNF", "GUCA2A", 
    "HHLA2", "IHH", "PSEN1", "NLGN", "NRTN", "RPH3A", "SHH","FLT3LG")

    receptor_list <- c("CD2", "CD27", "CD80", "CD86", "SELL", "CD44", "CD81", "CD8A", "CLEC1B", 
    "GLG1", "TYROBP", "FLT3", "ERBB2", "EGFR", "IL1R1", "IL1RAP", "KDR", "NRP1")

    if (db_type == "noncurated") {
    db <- import_db("noncurated")
    } else if (db_type == "curated") {
    db <- import_db("curated")
    } else if (db_type == "both") {
    db <- import_db("both")
    } else {
    stop("Invalid database type. Please choose 'noncurated', 'curated', or 'both'.")
    }

    pairwise_pairs <- create_pairwise_pairs(db)
    pt_interactions <- filter_pairs_with_ppi(pairwise_pairs)
    print("Number of PPI network interactions found:")
    print(nrow(pt_interactions))

    complete_data <- process_binary_pairs(db, pt_interactions)
    complete_data <- map_gene_data(complete_data)
    annotation <- annotate_components(complete_data)
    true_LR_DB <- process_lr_db(complete_data, annotation)
    adhesive_DB <- process_adhesive_DB(complete_data, annotation, ligand_list, receptor_list)

    LR_database <- rbind(true_LR_DB, adhesive_DB)

    return(LR_database)
}
