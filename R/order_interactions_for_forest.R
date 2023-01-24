#' @title order_interactions_for_forests
#' 
#' @description This function takes in a vector of indices and subsets the interactions dataframe based on the indices. Then it sorts the interactions based on log2FC_weights, log2FC_p_s_l, log2FC_phi_s_l, log2FC_rho_s, log2FC_p_r_r, log2FC_phi_r_r, log2FC_rho_r and interaction_ID columns. It also has helper functions like cluster_interactions and pick_param which are used in the sorting process.
#' 
#' @param my_idx A vector of indices for subsetting the interactions dataframe
#' 
#' @return A sorted dataframe of interactions
#' 
#' @export
#' @examples
#' # load example data
#' data("comm_result")
#' # calculate general statistics
#' comm_result <- general_stat(comm_result)
#' # order interactions
#' order_interactions_for_forests(my_idx = 1:10)
#' 
order_interactions_for_forests <- function(my_idx
){
        # subset anno_interactions
        my_anno_interactions <- interactions$anno_interactions[my_idx,c("log2FC_weights"
                                                                        ,"log2FC_rho_s"
                                                                        ,"log2FC_phi_s_l"
                                                                        ,"log2FC_p_s_l"
                                                                        ,"log2FC_rho_r"
                                                                        ,"log2FC_phi_r_r"
                                                                        ,"log2FC_p_r_r"
                                                                        ,"interaction_ID"
        )]
        rownames(my_anno_interactions) <- my_anno_interactions$interaction_ID
        #print(str(my_anno_interactions))
        
        comp_order <- c("none"
                        ,"one"
                        ,"both")
        params <- c("rho"
                    ,"phi"
                    ,"p")
        
        # help function
        # hclust on interactions
        # returns interactions ids vector
        cluster_interactions <- function(my_anno_interactions){
                #print((my_anno_interactions))
                #print(str(hclust(dist(my_anno_interactions$log2FC_weights))))
                
                # cluster
                # cluster by log2fc w
                my_anno_interactions  <-  my_anno_interactions[dendsort(hclust(dist(my_anno_interactions$log2FC_weights))
                                                                        ,isReverse = TRUE
                )$order,]
                #print(str(anno_interactions_sub))
                # cluster by log2fc p
                my_anno_interactions  <-  my_anno_interactions[dendsort(hclust(dist(my_anno_interactions[,c("log2FC_p_s_l"
                                                                                                            ,"log2FC_p_r_r")]))
                                                                        ,isReverse = TRUE
                )$order,]
                #print(str(anno_interactions_sub))
                # cluster by log2fc phi
                my_anno_interactions  <-  my_anno_interactions[dendsort(hclust(dist(my_anno_interactions[,c("log2FC_phi_s_l"
                                                                                                            ,"log2FC_phi_r_r")]))
                                                                        ,isReverse = TRUE
                )$order,]
                #print(str(anno_interactions_sub))
                # cluster by rho
                my_anno_interactions <- my_anno_interactions[dendsort(hclust(dist(my_anno_interactions[,c("log2FC_rho_s"
                                                                                                          ,"log2FC_rho_r")
                ]
                )
                )
                ,isReverse = TRUE
                )$order,]
                return(my_anno_interactions$interaction_ID)
                
                
        }
        
        # which columns to pick in anno_interactions. Returns a list with one time param for sending and one time param for receiving
        # returns list of columns names 
        pick_param <- function(my_param){
                if(my_param == "p"){
                        send <- "log2FC_p_s_l"
                        rec <- "log2FC_p_r_r"
                } else if(my_param == "phi"){
                        send <- "log2FC_phi_s_l"
                        rec <- "log2FC_phi_r_r"
                } else if(my_param == "rho"){
                        send <- "log2FC_rho_s"
                        rec <- "log2FC_rho_r"
                }
                return(list(send, rec))
        }
        
        # if parameter should be considered, first consider interactions with eiteher only s(_l) or r(_r) altered, then for both
        #returns index
        none_one_or_both <- function(x # "one" or "both"
                                     ,send # vector of log2FC values
                                     ,rec # vector of log2FC values
        ){
                ifelse(x == "none"
                       ,i <- (abs(send) <= 1) & (abs(rec) <= 1)
                       ,ifelse(x == "one"
                               ,i <- xor(abs(send) > 1
                                         ,abs(rec) > 1
                               )
                               ,i <- (abs(send) > 1) & (abs(rec) > 1)
                       )
                )
                #print(i)
                return(i)
        } 
        
        find_idx <- function(param
                             ,comp
                             ,my_anno_interactions){
                param_list <- pick_param(param)
                
                send_vector <- my_anno_interactions[[param_list[[1]]]]
                rec_vector <- my_anno_interactions[[param_list[[2]]]]
                
                none_one_or_both(x=comp # receive it from recursive over comp_order
                                 ,send= send_vector # receive it from pick_param func on recursive over params
                                 ,rec=rec_vector # receive it from pick_param func on recursive over params
                )
        }
        
        classes <- list(
                # none of the components
                c(rho = "rho",rho_comp = "none", phi = "phi", phi_comp = "none", p = "p", p_comp = "none")
                # only by rho
                ,c(rho = "rho",rho_comp = "one", phi = "phi", phi_comp = "none", p = "p", p_comp = "none")
                ,c(rho = "rho",rho_comp = "both", phi = "phi", phi_comp = "none", p = "p", p_comp = "none")
                # only by phi
                ,c(rho = "rho",rho_comp = "none", phi = "phi", phi_comp = "one", p = "p", p_comp = "none")
                ,c(rho = "rho",rho_comp = "none", phi = "phi", phi_comp = "both", p = "p", p_comp = "none")
                # only by p
                ,c(rho = "rho",rho_comp = "none", phi = "phi", phi_comp = "none", p = "p", p_comp = "one")
                ,c(rho = "rho",rho_comp = "none", phi = "phi", phi_comp = "none", p = "p", p_comp = "both")
                # by several
                ,c(rho = "rho",rho_comp = "none", phi = "phi", phi_comp = "one", p = "p", p_comp = "one")
                ,c(rho = "rho",rho_comp = "none", phi = "phi", phi_comp = "one", p = "p", p_comp = "both")
                ,c(rho = "rho",rho_comp = "none", phi = "phi", phi_comp = "both", p = "p", p_comp = "one")
                ,c(rho = "rho",rho_comp = "none", phi = "phi", phi_comp = "both", p = "p", p_comp = "both")
                ,c(rho = "rho",rho_comp = "one", phi = "phi", phi_comp = "none", p = "p", p_comp = "one")
                ,c(rho = "rho",rho_comp = "one", phi = "phi", phi_comp = "none", p = "p", p_comp = "both")
                ,c(rho = "rho",rho_comp = "one", phi = "phi", phi_comp = "one", p = "p", p_comp = "none")
                ,c(rho = "rho",rho_comp = "one", phi = "phi", phi_comp = "one", p = "p", p_comp = "one")
                ,c(rho = "rho",rho_comp = "one", phi = "phi", phi_comp = "one", p = "p", p_comp = "both")
                ,c(rho = "rho",rho_comp = "one", phi = "phi", phi_comp = "both", p = "p", p_comp = "none")
                ,c(rho = "rho",rho_comp = "one", phi = "phi", phi_comp = "both", p = "p", p_comp = "one")
                ,c(rho = "rho",rho_comp = "one", phi = "phi", phi_comp = "both", p = "p", p_comp = "both")
                ,c(rho = "rho",rho_comp = "both", phi = "phi", phi_comp = "none", p = "p", p_comp = "one")
                ,c(rho = "rho",rho_comp = "both", phi = "phi", phi_comp = "none", p = "p", p_comp = "both")
                ,c(rho = "rho",rho_comp = "both", phi = "phi", phi_comp = "one", p = "p", p_comp = "none")
                ,c(rho = "rho",rho_comp = "both", phi = "phi", phi_comp = "one", p = "p", p_comp = "one")
                ,c(rho = "rho",rho_comp = "both", phi = "phi", phi_comp = "one", p = "p", p_comp = "both")
                ,c(rho = "rho",rho_comp = "both", phi = "phi", phi_comp = "both", p = "p", p_comp = "none")
                ,c(rho = "rho",rho_comp = "both", phi = "phi", phi_comp = "both", p = "p", p_comp = "one")
                ,c(rho = "rho",rho_comp = "both", phi = "phi", phi_comp = "both", p = "p", p_comp = "both")
        )
        
        interactions_order <- unlist(lapply(classes
                                            ,function(c){
                                                    # idx rho
                                                    idx_rho <- find_idx(param = c["rho"]
                                                                        ,comp = c["rho_comp"]
                                                                        ,my_anno_interactions
                                                    )
                                                    
                                                    # idx phi
                                                    idx_phi <- find_idx(param = c["phi"]
                                                                        ,comp = c["phi_comp"]
                                                                        ,my_anno_interactions
                                                    )
                                                    
                                                    # idx p
                                                    idx_p <- find_idx(param = c["p"]
                                                                      ,comp = c["p_comp"]
                                                                      ,my_anno_interactions
                                                    )
                                                    
                                                    # idx
                                                    my_running_idx <- idx_rho & idx_phi & idx_p
                                                    
                                                    
                                                    # subset anno_interactions
                                                    my_anno_interactions_sub <- my_anno_interactions[my_running_idx,]
                                                    
                                                    # print out class and nr interactions
                                                    print(paste(c))
                                                    print(paste(nrow(my_anno_interactions_sub), "interactions"))
                                                    
                                                    # interactions_order
                                                    if(nrow(my_anno_interactions_sub) > 1){
                                                            return(cluster_interactions(my_anno_interactions_sub))
                                                    } else if(nrow(my_anno_interactions_sub) == 1){
                                                            return(my_anno_interactions_sub$interaction_ID)
                                                    }
                                            }
        )
        )
        # reverse order for forest plot
        interactions_order <- interactions_order[length(interactions_order):1]
        
        #print("(interactions_order)")
        #print((interactions_order[1:10]))
        
        # sort my_anno_interactions
        my_anno_interactions <- my_anno_interactions[interactions_order,]
        rownames(my_anno_interactions) <- my_anno_interactions$interaction_ID
        my_anno_interactions$interaction_ID <- factor(my_anno_interactions$interaction_ID
                                                       ,levels = my_anno_interactions$interaction_ID
                                                       ,ordered = TRUE)
        
        #print("(my_anno_interactions)")
        #print((my_anno_interactions))
        
        return(my_anno_interactions)
        
}