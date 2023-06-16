#' @title interaction_classes
#' @description This function adds columns to the annotation of the interactions
#' @param my_interactions list of interactions (see \code{\link{interaction_list}}) 
#' @param threshold_log2FC_components threshold for the log2FC of the components
#' @return my_interactions list of interactions (see \code{\link{interaction_list}})
#' @export
#' @examples
#' # load the data
#' data("my_interactions")
#' # add the columns
#' my_interactions <- interaction_classes(my_interactions
#'                                      ,threshold_log2FC_components = 1)
#' # show the first 10 rows
#' head(my_interactions$anno_interactions
#'    ,n = 10)
#' # show the last 10 rows
#' tail(my_interactions$anno_interactions
#'   ,n = 10)
#' 


interaction_classes <- function(my_interactions
                                ,threshold_log2FC_components = 1
                               ){
    
    # store threshol in the threshold list
    my_interactions$thresholds$threshold_log2FC_components <- threshold_log2FC_components
    
    if(!("log2FC_rho_s_direction" %in% colnames(my_interactions$anno_interactions))){
        df <- my_interactions$anno_interactions[,c("log2FC_rho_s"
                                     ,"log2FC_phi_s_l"
                                     ,"log2FC_p_s_l"
                                     ,"log2FC_rho_r"
                                     ,"log2FC_phi_r_r"
                                     ,"log2FC_p_r_r"
                                     )
                                  ]
        #rownames(df) <- my_anno_interactions$interaction_ID
        colnames(df) <- paste(colnames(df)
                             ,"direction"
                             ,sep = "_")

        # binarize threshold
        df[abs(df) < threshold_log2FC_components] <- 0 # unshanged
        df[df > threshold_log2FC_components] <- 1 # upregulated
        df[df < -threshold_log2FC_components] <- -1 # downregulated

        my_interactions$anno_interactions <- cbind(my_interactions$anno_interactions
                                                  ,df)

        rm(df)
    }
    
    
    
    ### ________________________
    #
    # Number of components affected (per interacting partner)
    #
    # values can be: 
    #    - 0,1,2,3 for sender and receiver
    #    - 0,1,2,3,4,5,6 for both
    
    my_interactions$anno_interactions$nr_comp_affected_s <- rowSums(abs(my_interactions$anno_interactions[,c("log2FC_rho_s_direction"
                                                                                                             ,"log2FC_phi_s_l_direction"
                                                                                                             ,"log2FC_p_s_l_direction")
                                                                                                           ]
                                                                        )
                                                                    )
    my_interactions$anno_interactions$nr_comp_affected_r <- rowSums(abs(my_interactions$anno_interactions[,c("log2FC_rho_r_direction"
                                                                                                             ,"log2FC_phi_r_r_direction"
                                                                                                             ,"log2FC_p_r_r_direction")
                                                                                                           ]
                                                                        )
                                                                    )
    my_interactions$anno_interactions$nr_comp_affected_b <- rowSums(abs(my_interactions$anno_interactions[,c("log2FC_rho_s_direction"
                                                                                                             ,"log2FC_phi_s_l_direction"
                                                                                                             ,"log2FC_p_s_l_direction"
                                                                                                             ,"log2FC_rho_r_direction"
                                                                                                             ,"log2FC_phi_r_r_direction"
                                                                                                             ,"log2FC_p_r_r_direction"
                                                                                                             )
                                                                                                          ]
                                                                        )
                                                                    )
    
    
    
    ### ________________________
    #
    # List of affected components (per interacting partner)
    #
    # values can be: 
    #    - none
    #    - p
    #    - phi
    #    - rho
    #    - phi_p
    #    - rho_p
    #    - rho_phi
    #    - rho_phi_p
    
    # iniciate column components_affected with an empty string
    my_interactions$anno_interactions$components_affected_s <- ""
    my_interactions$anno_interactions$components_affected_r <- ""
    
    for(type in c("s" # "sender"
                 ,"r" # "receiver"
                 )){
        which_components_affected_column <- paste0("components_affected_"
                                                   ,type
                                                   )
        
        # check individual complnents
    for(param in c("rho"
                 ,"phi"
                 ,"p")
       ){
        which_log2FC_column <- grep(paste0("log2FC_"
                                   ,param
                                   ,"_"
                                   )
                             ,colnames(my_interactions$anno_interactions)
                            )# vector of two: first index is for sending, second index is for receiving
        
        ifelse(type == "s", i <- 1, i <- 2)
        idx_param <- (abs(my_interactions$anno_interactions[[which_log2FC_column[i]]]) > my_interactions$thresholds$threshold_log2FC_components)
        
        my_interactions$anno_interactions[[which_components_affected_column]][idx_param] <- paste(my_interactions$anno_interactions[[which_components_affected_column]][idx_param]
                                                                                ,param)
        
    }
    
    # for those interactions that don't have any component changed, put " none"
    idx_none <- my_interactions$anno_interactions[[which_components_affected_column]] == ""
    my_interactions$anno_interactions[[which_components_affected_column]][idx_none] <- " none"
    
    # remove space in from of every strign
    my_interactions$anno_interactions[[which_components_affected_column]] <- sub("."
                                                                                  ,""
                                                                                  ,my_interactions$anno_interactions[[which_components_affected_column]])
    
    }
    
    
    
    ### ________________________
    #
    # Which counterpart is affected
    #
    # values can be:
    #
    # - none if nr component sender is = 0 and nr component receiver is = 0
    # - sender if nr component sender is > 0 and nr component receiver is = 0
    # - receiver if nr component sender is = 0 and nr component receiver is > 0
    # - both if nr component sender is > 0 and nr component receiver is > 0
    
    my_interactions$anno_interactions$sender_or_receiver_affected <- "none"
    
    idx_s <- (my_interactions$anno_interactions$nr_comp_affected_s > 0) & (my_interactions$anno_interactions$nr_comp_affected_r == 0)
    idx_r <- (my_interactions$anno_interactions$nr_comp_affected_s == 0) & (my_interactions$anno_interactions$nr_comp_affected_r > 0)
    idx_b <- (my_interactions$anno_interactions$nr_comp_affected_s > 0) & (my_interactions$anno_interactions$nr_comp_affected_r > 0)
    
    my_interactions$anno_interactions$sender_or_receiver_affected[idx_s] <- "sender"
    my_interactions$anno_interactions$sender_or_receiver_affected[idx_r] <- "receiver"
    my_interactions$anno_interactions$sender_or_receiver_affected[idx_b] <- "both"
    # store as ordered factor: sender > receiver > both
    my_interactions$anno_interactions$sender_or_receiver_affected <- factor(my_interactions$anno_interactions$sender_or_receiver_affected
                                                                            ,levels = c("none"
                                                                                        ,"sender"
                                                                                        ,"receiver"
                                                                                        ,"both"
                                                                                        )
                                                                            ,ordered = TRUE)
    
    
    
    
    ###__________________
    #
    # Direction of affected components (per interacting partner)
    #
    # values can be: 
    #
    # - "none" if zero components affected
    # - "up" if >0 components affected in the same difrection (upregulated)
    # - "down" if >0 components affected in the same difrection (downregulated)
    # - "both"  if >1 components affected in the opposite difrection
    
    sum_s <- rowSums(my_interactions$anno_interactions[,c("log2FC_rho_s_direction"
                                                         ,"log2FC_phi_s_l_direction"
                                                         ,"log2FC_p_s_l_direction")
                           ])
    sum_r <- rowSums(my_interactions$anno_interactions[,c("log2FC_rho_r_direction"
                                                         ,"log2FC_phi_r_r_direction"
                                                         ,"log2FC_p_r_r_direction")
                           ])
    sum_b <- rowSums(my_interactions$anno_interactions[,c("log2FC_rho_s_direction"
                                                         ,"log2FC_phi_s_l_direction"
                                                         ,"log2FC_p_s_l_direction"
                                                         ,"log2FC_rho_r_direction"
                                                         ,"log2FC_phi_r_r_direction"
                                                         ,"log2FC_p_r_r_direction"
                                                         )
                      ])
    
    for(cont in c("s" # sender
                 ,"r" # receiver
                 ,"b" # both
                 )
       ){
        ifelse(cont == "s"
              ,{my_sum <- sum_s
               nr_comp_affected <- "nr_comp_affected_s"}
               ,ifelse(cont == "r"
                     ,{my_sum <- sum_r
                       nr_comp_affected <- "nr_comp_affected_r"}
                      ,{my_sum <- sum_b
                       nr_comp_affected <- "nr_comp_affected_b"}
                      )
              )
        
          my_col_name <- paste("direction"
                               ,cont
                               ,sep = "_")
        
        if(my_col_name %in% colnames(my_interactions$anno_interactions)){
            my_interactions$anno_interactions[[my_col_name]] <- "none"
        } else {
            my_interactions$anno_interactions$direction <- "none"
            colnames(my_interactions$anno_interactions)[colnames(my_interactions$anno_interactions) == "direction"] <- my_col_name
        }
        
        idx_up <- (my_sum > 0 # at leas one component is upregulated
                  ) & (abs(my_sum) == my_interactions$anno_interactions[[nr_comp_affected]]) # .. in the same direction
        idx_down <- (my_sum < 0 # at leas one component is downregulated
                    ) & (abs(my_sum) == my_interactions$anno_interactions[[nr_comp_affected]]) # .. in the same direction
        idx_both <- abs(my_sum) != my_interactions$anno_interactions[[nr_comp_affected]] # .. in the same direction
        
        my_interactions$anno_interactions[[my_col_name]][idx_up] <- "up"
        my_interactions$anno_interactions[[my_col_name]][idx_down] <- "down"
        my_interactions$anno_interactions[[my_col_name]][idx_both] <- "both"
        # store as ordered factor: none > up > down > both
        my_interactions$anno_interactions[[my_col_name]] <- factor(my_interactions$anno_interactions[[my_col_name]]
                                                                     ,levels = c("none", "up", "down", "both")
                                                                     ,ordered = TRUE
                                                                    )
        # rename column
        
    }
    
    
    
    
    ###__________________
    #
    # Concordance/disconcordance of direction of affected components (per interacting partner)
    #
    # values can be:
    #
    # - NA if:
    #     - zero components are affected
    #     - one component is affected
    #
    # - concordant if:
    #     - several components are affected in the same direction
    #
    # - disconcordant if:
    #     for concordance_s and concordance_r if:
    #         - several components are affected in different directions
    #     for concordance_b if:
    #         - concordance_s concordant (direction_s up), concordance_r concordant (direction_r down), 
    #         - concordance_s concordant (direction_s down), concordance_r concordant (direction_r up), 
    #         - concordance_s concordant, concordance_r disconcordant,
    #         - concordance_s disconcordant, concordance_r concordant,
    #         - concordance_s disconcordant, concordance_r disconcordant
    
    for(cont in c("s" # sender
                 ,"r" # receiver
                 ,"b" # both
                 )
       ){
          my_col_name <- paste("concordance"
                               ,cont
                               ,sep = "_")
        
        ifelse(cont == "s"
              ,{
                  idx_notnone <- my_interactions$anno_interactions$nr_comp_affected_s >1
                  idx_different <- abs(sum_s) != my_interactions$anno_interactions$nr_comp_affected_s 
              }
              ,ifelse(cont == "r"
                     ,{
                         idx_notnone <- my_interactions$anno_interactions$nr_comp_affected_r >1
                         idx_different <- abs(sum_r) != my_interactions$anno_interactions$nr_comp_affected_r 
                     }
                     ,{
                         idx_notnone <- rowSums(my_interactions$anno_interactions[,c("nr_comp_affected_s"
                                                                                     ,"nr_comp_affected_r")])>1
                         idx_different <- (my_interactions$anno_interactions$direction_s == "up" & my_interactions$anno_interactions$direction_r == "down"
                                         ) | (my_interactions$anno_interactions$direction_s == "down" &  my_interactions$anno_interactions$direction_r == "up"
                                         ) | (my_interactions$anno_interactions$concordance_s == "disconcordant"
                                         ) | (my_interactions$anno_interactions$concordance_r == "disconcordant")
                     })
              )
        
        if(my_col_name %in% colnames(my_interactions$anno_interactions)){
            my_interactions$anno_interactions[[my_col_name]] <- "undefined" # some of these values will be changed in the next lines of code
        } else {
            my_interactions$anno_interactions$concordance <-"undefined" # some of these values will be changed in the next lines of code
            colnames(my_interactions$anno_interactions)[colnames(my_interactions$anno_interactions) == "concordance"] <- my_col_name
        }
        
        my_interactions$anno_interactions[[my_col_name]][idx_notnone] <- "concordant" # some of these values will be changed in the next lines of code
        my_interactions$anno_interactions[[my_col_name]][idx_different] <- "disconcordant" 
        # store as ordered factor: NA -> concordant -> disconcordant
        my_interactions$anno_interactions[[my_col_name]] <- factor(my_interactions$anno_interactions[[my_col_name]]
                                                                ,levels = c("undefined","concordant","disconcordant")
                                                                ,ordered = TRUE) 
    }
    
    
    
    
    
    ###__________________
    #
    # Interaction category
    #
    # values can be:
    # 
    # - no_change if:
    #    - absolute log2FC_weight ≤ threshold_log2FC and nr_comp_affected_b == 0
    #
    # - simple_decrease if:
    #    - log2FC_weight < -threshold_log2FC and direction_b == "down" and  nr_comp_affected_b == 1
    # 
    # - simple_increase if:
    #    - log2FC_weight > threshold_log2FC and direction_b == "up" and  nr_comp_affected_b == 1
    # 
    # - concordant_decrease if:
    #    - log2FC_weight < -threshold_log2FC and direction_b == "down" and  nr_comp_affected_b >1
    #
    # - concordant_increase if:
    #    - log2FC_weight > threshold_log2FC and direction_b == "up" and  nr_comp_affected_b >1
    #
    # - insufficient_compensation if:
    #    - absolute log2FC_weight > threshold_log2FC and direction_b == "both" and  nr_comp_affected_b >1
    #
    # - sufficient_compensation if:
    #    - absolute log2FC_weight ≤ threshold_log2FC and direction_b == "both" and  nr_comp_affected_b >1
    #
    # - undefined if: none of the above cases apply
    
    threshold_log2FC <- my_interactions$thresholds$threshold_log2FC 
    idx_up <- my_interactions$anno_interactions$log2FC_weights > threshold_log2FC
    idx_no_change <- abs(my_interactions$anno_interactions$log2FC_weights)<= threshold_log2FC
    idx_down <- my_interactions$anno_interactions$log2FC_weights < -threshold_log2FC

    
    idx_dir_b_down <- my_interactions$anno_interactions$direction_b == "down"
    idx_dir_b_up <- my_interactions$anno_interactions$direction_b == "up"
    idx_dir_b_both <- my_interactions$anno_interactions$direction_b == "both"

    
    idx_nrComp_s_none <- my_interactions$anno_interactions$nr_comp_affected_s == 0
    idx_nrComp_s_one <- my_interactions$anno_interactions$nr_comp_affected_s == 1
    idx_nrComp_s_several <- my_interactions$anno_interactions$nr_comp_affected_s >1

    idx_nrComp_r_none <- my_interactions$anno_interactions$nr_comp_affected_r == 0
    idx_nrComp_r_one <- my_interactions$anno_interactions$nr_comp_affected_r == 1
    idx_nrComp_r_several <- my_interactions$anno_interactions$nr_comp_affected_r >1
    
    idx_nrComp_b_none <- my_interactions$anno_interactions$nr_comp_affected_b == 0
    idx_nrComp_b_one <- my_interactions$anno_interactions$nr_comp_affected_b == 1
    idx_nrComp_b_several <- my_interactions$anno_interactions$nr_comp_affected_b >1
    
    my_interactions$anno_interactions$interaction_category <- "undefined"
    
    idx_unchanged <- idx_no_change & idx_nrComp_b_none
    idx_simple_down <- idx_down & idx_dir_b_down & idx_nrComp_b_one
    idx_simple_up <- idx_up & idx_dir_b_up & idx_nrComp_b_one
    idx_concordant_down <- idx_down & idx_dir_b_down & idx_nrComp_b_several
    idx_concordant_up <- idx_up & idx_dir_b_up & idx_nrComp_b_several
    idx_insufficient_comp <- !idx_no_change & idx_dir_b_both# & idx_nrComp_b_several
    idx_sufficient_comp <- idx_no_change & idx_dir_b_both# & idx_nrComp_b_several

    my_interactions$anno_interactions$interaction_category[idx_unchanged] <- "no_change"
    my_interactions$anno_interactions$interaction_category[idx_simple_down] <- "simple_decrease"
    my_interactions$anno_interactions$interaction_category[idx_simple_up] <- "simple_increase"
    my_interactions$anno_interactions$interaction_category[idx_concordant_down] <- "concordant_decrease"
    my_interactions$anno_interactions$interaction_category[idx_concordant_up] <- "concordant_increase"
    my_interactions$anno_interactions$interaction_category[idx_insufficient_comp] <- "insufficient_compensation"
    my_interactions$anno_interactions$interaction_category[idx_sufficient_comp] <- "sufficient_compensation"
    
    
    
    
    # return the my_interactions object
    return(my_interactions)
    
}