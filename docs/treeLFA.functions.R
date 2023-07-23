### Training with the gibbs-em algorithm: 
gibbs_EM_train <- function( topic.number,
                            data, tree_str,
                            Phi, I, rho, alpha, Z,
                            burn_in,
                            opt_N_1, cycle_1,
                            opt_N_2, cycle_2 ) {
  
  
  # Prepare the input data:
  S <- ncol(data)     # number of terminal disease codes on the tree
  S1 <- nrow(tree_str) - ncol(data)    # number of internal disease codes on the tree
  # S1 <- 5
  data <- as.matrix(data)
  
  
  
  
  
  ## Check the input diagnosis data: 
  # Check the input matrix is binary: 
  if ( all( names(table(data)) == c("0","1") ) == FALSE ) { 
    stop("The diagnosis data is not binary (coded using 0 and 1)!")
  }
  
  # Check there are no repetitive disease codes in the input diagnosis matrix: 
  if ( length(unique(colnames(data))) < ncol(data) ) { 
    stop("There are repetitive disease codes in the input data!")
  }
  
  
  
  
  
  ### Tree structure: 
  ## Check the column names:
  if ( all( colnames(tree_str) != c("node","parent","layer","terminal") ) ) { 
    stop("The columns of the tree structure matrix should be set as required!")  
  }

  # Check repetitive codes in the 1st column: 
  if ( length(unique(tree_str$node)) != nrow(tree_str) ) { 
    stop("There are repetitive disease codes in the 1st column!")
  }
  
  # Check root nodes in the 2nd column of the tree: 
  if ( !(any(tree_str$parent=="root")) ) { 
    stop("No root node in the 2nd column coded as \"root\"!")
  }
  
  # Check layers of disease codes
  if (setequal(  unique(tree_str$layer), seq( 2,(length(unique(tree_str$layer))+1) )  )==FALSE ) { 
    stop("Layers of disease codes should be a series of integers starting from 2 ")
  }
  
  # Check the coding of terminal codes: 
  if ( setequal( unique(tree_str$terminal), c("Y","N") ) == FALSE ) { 
    stop("The terminal/internal codes should be indicated with \"Y\" and \"N\"")
  }
  
  
  
  # Terminal codes on the tree:
  tree_ter <- tree_str[terminal=="Y",]
  
  # Check that terminal codes on the tree structure and codes in the input data match:
  if ( setequal(tree_ter$node, colnames(data))==FALSE ) { 
    stop("The terminal codes in the tree structure matrix and disease codes in the input data do not match!")
  }
  
  
  ## Reorder terminal/internal codes: 
  # internal codes:
  tree_int <- tree_str[!(node %in% colnames(data)),]

  # Reorder internal codes: according to layers of codes 
  tree_int<- tree_int[order(layer),]
  
  # Reorder termianl codes: according to the columns (order of disease codes) of the input data
  tree_ter <- tree_ter[match(colnames(data),tree_ter$node)]
  
  tree_str <- rbind(tree_int,tree_ter)
  
  
  
  
  # Check that the parent code of each code is on the layer above: 
  for ( i in 1:nrow(tree_str) ) { 
      
    if ( tree_str$layer[i]==2 ) { 
      if ( tree_str$parent[i]!="root" ) { 
        stop( "The parent code of disease codes on the 2nd layer must be the root code (coded as \"root\")!" ) 
      }
    } else { 
      
      # layer of parent node:
      pa <- unlist( tree_str[node==tree_str$parent[i],"node"] )
      
      if ( unlist(tree_str[node==pa,"layer"]) != (tree_str$layer[i]-1) ) { 
        stop( paste("The parent code of code", paste("\"",tree_str$node[i],"\"",sep=""), "is not on the above layer!") )
      }
      
    }
  
  }
  
  tree_str_code <- tree_str
    
  
  
  # Prepara the tree str: change the names of diseases codes to indexes:
 tree_str <- tree_str[,1:2]
 
  for ( i in 1:nrow(tree_str) ){

    node_name <- tree_str$node[i]
    tree_str$node[i] <- i

    tree_str[parent==node_name,"parent"] <- i

  }

  tree_str[parent=="root","parent"] <- 0
  
  tree_str$node <- as.double(tree_str$node)
  tree_str$parent <- as.double(tree_str$parent)
  tree_str <- as.matrix(tree_str)
  
  
  
  
  
  
  ## Set the hyper-parameters:
  # Beta prior for phi based on incidence of diseases:
  inc <- min(colSums(data)) / nrow(data)

  if ( inc > 0.01 ) {
    a00 <- 0.1
    a01 <- 500
    a10 <- 2
    a11 <- 4
  } else if ( inc <= 0.01 & inc > 0.005 ) {
    a00 <- 0.1
    a01 <- 1000
    a10 <- 1.5
    a11 <- 3
  } else {
    a00 <- 0.1
    a01 <- 4000
    a10 <- 1.2
    a11 <- 3
  }


  K <- topic.number    # number of topics 
  D <- nrow(data)    # number of individuals

  # Beta prior for the transition probability of the Markov process on the tree:
  b00 <- 3
  b01 <- 20
  b10 <- 3
  b11 <- 3

  
  
  
  # Initialization of hidden variables:
  if(missing(rho)) { 
    rho <- c(0.9,0.5)
  } 
  
  if(missing(Z)) { 
  Z <- matrix( nrow=D, ncol=S )
  for ( d in 1:D ){
    if ( sum(data[d,]) == 0 ) {
      Z[d,] <- rep(1,S)
    } else {
      for ( l in 1:S ) { Z[d,l] <- rcat( rep(1,K) )}
    }
  }
  #Z[1:D,1:S1] <- 0
  } 
  
  if(missing(I)) { 
  #I <- matrix( nrow=K, ncol=(S), 0 )
  I <- matrix( nrow=K, ncol=(S+S1), 0 )
  }
  
  if(missing(Phi)) { 
  Phi <- matrix( nrow=K, ncol=S )
  for ( j in 1:length(Phi) ) { Phi[j] <- rbeta(1,1,5000000) }
  } 
  
  
  
  
  
  
  ## Check the parameters provided: 
  # the alpha vector and number of topics should be of the same length
  if ( length(alpha) != K ) { 
    stop("The length of alpha should be the same as the number of topics to be inferred!")
  }
  
  if ( all(alpha>0) == FALSE | is.double(alpha) == FALSE ) { 
    stop("alpha should be a vector of positive real numbers!")
  }
  
  if ( K<=0 | ( K - floor(K) != 0 ) ) { 
    stop("The number of topics should be a positive integer!")
  }
  
  if ( burn_in < 0 |  ( burn_in - floor(burn_in) != 0 ) ) { 
    stop("A postive integer should be set for the \"burn_in\" argument!")
  }
  
  if ( opt_N_1 < 0 |  ( opt_N_1 - floor(opt_N_1) != 0 ) ) { 
    stop("A postive integer should be set for the \"opt_N_1\" argument!")
  }
  
  if ( cycle_1 < 0 |  ( cycle_1 - floor(cycle_1) != 0 ) ) { 
    stop("A postive integer should be set for the \"cycle_1\" argument!")
  }
  
  if ( opt_N_2 < 0 |  ( opt_N_2 - floor(opt_N_2) != 0 ) ) { 
    stop("A postive integer should be set for the \"opt_N_2\" argument!")
  }
  
  if ( cycle_2 < 0 |  ( cycle_2 - floor(cycle_2) != 0 ) ) { 
    stop("A postive integer should be set for the \"cycle_2\" argument!")
  }
  
  
  
  
  
  # Run the gibbs-em algorithm:
  result <- gibbs_EM( K, S, S1, D,
            a00, a01, a10, a11,
            b00, b01, b10, b11,
            alpha,
            rho,
            Z, I, Phi,
            data,
            tree_str,
            burn_in,
            opt_N_1, cycle_1,
            opt_N_2, cycle_2 )

  names(result) <- c("Phi_samples","I_samples","rho_samples","alpha_samples","Z_samples","L_all")

  return(result)

}










### Training with the gibbs sampling algorithm: alpha fixed at values given by gibbs-em algorithm:
gibbs_train <- function( data, tree_str,
                         Phi, I, rho, alpha, Z,
                         burn_in, cycle, interval ) {
  
  # Prepare the input data:
  topic.number <- nrow(Phi)
  S <- ncol(data)    # number of terminal disease codes on the tree 
  S1 <- nrow(tree_str) - ncol(data)      # number of internal disease codes on the tree 
  # S1 <- 5
  data <- as.matrix(data)
  
  
  
  
  ## Process the tree structure: 
  # Terminal codes on the tree:
  tree_ter <- tree_str[terminal=="Y",]
  
  ## Reorder terminal/internal codes: 
  # internal codes:
  tree_int <- tree_str[!(node %in% colnames(data)),]
  
  # Reorder internal codes: according to layers of codes 
  tree_int<- tree_int[order(layer),]
  
  # Reorder termianl codes: according to the columns (order of disease codes) of the input data
  tree_ter <- tree_ter[match(colnames(data),tree_ter$node)]
  
  tree_str <- rbind(tree_int,tree_ter)
  
  
  
  # Prepara the tree str: change the names of diseases codes to indexes:
  tree_str <- tree_str[,1:2]
  for ( i in 1:nrow(tree_str) ){
    
    node_name <- tree_str$node[i]
    tree_str$node[i] <- i
    
    tree_str[parent==node_name,"parent"] <- i
    
  }
  
  tree_str[parent=="root","parent"] <- 0
  
  tree_str$node <- as.double(tree_str$node)
  tree_str$parent <- as.double(tree_str$parent)
  tree_str <- as.matrix(tree_str)
  
  
  
  
  # Set the hyper-parameters:
  # Beta prior for phi based on the prevalences of diseases:
  inc <- min(colSums(data)) / nrow(data)
  
  if ( inc > 0.01 ) {
    a00 <- 0.1
    a01 <- 500
    a10 <- 2
    a11 <- 4
  } else if ( inc <= 0.01 & inc > 0.005 ) {
    a00 <- 0.1
    a01 <- 1000
    a10 <- 1.5
    a11 <- 3
  } else {
    a00 <- 0.1
    a01 <- 4000
    a10 <- 1.2
    a11 <- 3
  }
  
  K <- topic.number    # number of topics 
  D <- nrow(data)    # number of individuals
  
  # Beta prior for the transition probability of the Markov process on the tree
  b00 <- 3
  b01 <- 20
  b10 <- 20
  b11 <- 3
  

  
  
  
  # Check the parameters: 
  if ( burn_in < 0 |  ( burn_in - floor(burn_in) != 0 ) ) { 
    stop("A postive integer should be set for the \"burn_in\" argument!")
  }
  
  if ( cycle < 0 |  ( cycle - floor(cycle) != 0 ) ) { 
    stop("A postive integer should be set for the \"cycle\" argument!")
  }
  
  if ( interval < 0 |  ( interval - floor(interval) != 0 ) ) { 
    stop("A postive integer should be set for the \"interval\" argument!")
  }
  
  if ( burn_in >= cycle ) { 
    stop("\"burn_in\" should be smaller then \"cycle\"!")
  }
  
  if ( ((cycle-burn_in) %% interval) != 0 ) { 
    stop("(\"cycle\"-\"burn_in\") should be devisible by \"interval\"!")
  }
  
  
  
  
  
  # Run the gibbs-em algorithm:
  result <- gibbs( K, S, S1, D,
                   a00, a01, a10, a11,
                   b00, b01, b10, b11,
                   alpha,
                   rho,
                   Z, I, Phi,
                   data,
                   tree_str,
                   opt_N=1, burn_in, cycle, interval )
  names(result) <- c("Phi_samples","I_samples","rho_samples","Z_sum_samples","alpha","L_all")  
  return(result)
  
}










## Calculate predictive likelihood: 
predL <- function(alpha, phi, data, IS){
  
  data <- as.matrix(data)
  
  # Check the input matrix is binary: 
  if ( all( names(table(data)) == c("0","1") ) == FALSE ) { 
    stop("The diagnosis data is not binary (coded using 0 and 1)!")
  }
  
  # Check there are no repetitive disease codes in the input diagnosis matrix: 
  if ( length(unique(colnames(data))) < ncol(data) ) { 
    stop("There are repetitive disease codes in the input data!")
  }
  
  # Check terminal codes on the tree and disease codes in the input data match:
  if ( IS < 0 |  ( IS - floor(IS) != 0 ) ) { 
    stop("A postive integer should be set for the \"IS\" argument!")
  }
  
  
  pl <- pred_L_mca( alpha, phi, data, IS=100 )
  
  return(pl)
  
}









# plot topics using pheatmap: 
topics_plot <- function(phi) { 
  
  for ( i in 1:nrow(phi) ) { 
    phi[i,] <- rev(phi[i,])  
  }
  
  rownames(phi) <- paste("Topic",1:nrow(phi),sep="")
  colnames(phi) <- paste("Disease",ncol(phi):1,sep="")
    
  pheatmap( t(phi),
            cluster_cols=FALSE,
            cluster_rows=FALSE,
            color = colorRampPalette(brewer.pal(9,"Blues"))(300) )
  
}



# plot topics using pheatmap: 
tw_plot <- function(tw) { 
  
  for ( i in 1:nrow(tw) ) {
    tw[i,] <- rev(tw[i,])
  }
  
  rownames(tw) <- paste("Individual",1:nrow(tw),sep="")
  colnames(tw) <- paste("Topic",ncol(tw):1,sep="")
  
  pheatmap( tw,
            cluster_cols=FALSE,
            cluster_rows=FALSE,
            color = colorRampPalette(brewer.pal(9,"Reds"))(300) )

}



















### Cluster posterior sample sof topics using Louvain algorithm:
cls_louvain <- function(g_result,k) { 
  
  if ( k < 0 |  ( k - floor(k) != 0 ) ) { 
    stop("A postive integer should be set for \"k\"!")
  }
  
  # parameters: 
  K <- nrow(g_result$Phi_samples[[1]])
  D <- nrow( g_result$Z_sum_samples[[1]] )
  S <- ncol(g_result$Phi_samples[[1]])
  N_ps <- length(g_result$Z_sum_samples)
  
  
  # All posteior samples of topics: 
  topics_phi <- matrix(nrow=0, ncol=S)
  #alpha_all <- vector(mode="double",length=0)
  
  for ( ps in 1:N_ps ) { 
    topics_phi <- rbind( topics_phi, g_result$Phi_samples[[ps]] )
    #alpha_all <- c(alpha_all,g_result$alpha)
  }
  
  
  
  ## Louvain clustering of poterior topics: 
  # build KNN graph:  
  g <- buildSNNGraph(t(topics_phi), k=k, type="number")
  
  # louvain: 
  clust_all <- igraph::cluster_louvain(g)$membership
  
  
  
  # Averaged topics for clusters:  
  # topics averaged from multiple posterior samples: 
  topics_ave <- matrix(nrow=length(table(clust_all)),ncol=ncol(topics_phi))
  
  for ( cl in 1:length(table(clust_all)) ){ 
    topics_phi_cl <- topics_phi[clust_all==cl,]
    topic_phi_ave <- colMeans(topics_phi_cl)
    topics_ave[cl,] <- topic_phi_ave
  }


  
  # return results: 
  g_result$clust <- clust_all
  g_result$topics_ave <- topics_ave
  return(g_result)
  
}













### Cluster posterior sample sof topics using Louvain algorithm:
cls_hier <- function(g_result,k) { 
  
  if ( k < 0 |  ( k - floor(k) != 0 ) ) { 
    stop("A postive integer should be set for \"k\"!")
  }
  
  if ( k > length(table(g_result$clust)) ) { 
    stop("the number of topics set is larger than the current number of clusters!")
  }
  
  
  # parameters: 
  K <- nrow(g_result$Phi_samples[[1]])
  D <- nrow( g_result$Z_sum_samples[[1]] )
  S <- ncol(g_result$Phi_samples[[1]])
  N_ps <- length(g_result$Z_sum_samples)
  
  topics_ave <- g_result$topics_ave
  clust_all <- g_result$clust
  
  
  ## Hierarchical clustering of poterior topics: 
  dist_topics <- dist(topics_ave,method="manhattan")
  cls_h <- hclust( dist_topics )
  clust_h <- cutree( cls_h, k=k )
  
  
  # Averaged topics for new clusters:  
  topics_ave_h <- matrix(nrow=length(table(clust_h)),ncol=ncol(topics_ave),0)
  
  for ( cl in 1:length(table(clust_h)) ){
    topics_h_tmp <- topics_ave[clust_h==cl,,drop=FALSE]
    if ( nrow(topics_h_tmp) == 1 ) { 
      topic_h_tmp <- topics_h_tmp
    } else { 
      topic_h_tmp <- colMeans(topics_h_tmp)
    }
    topics_ave_h[cl,] <- topic_h_tmp
  }
  
  
  # return results: 
  g_result$clust_h <- clust_h
  g_result$topics_ave_h <- topics_ave_h
  return(g_result)
  
}










### Cluster posterior sample sof topics using Louvain algorithm:
cls_others <- function(g_result) { 
  
  # parameters: 
  K <- nrow(g_result$Phi_samples[[1]])
  D <- nrow( g_result$Z_sum_samples[[1]] )
  S <- ncol(g_result$Phi_samples[[1]])
  N_ps <- length(g_result$Z_sum_samples)
  
  clust_all <- g_result$clust
  
  if ( exists("clust_h",g_result) ) { 
    clust_h <- g_result$clust_h
  } else { 
    clust_h <- unique(g_result$clust)
    g_result$clust_h <- clust_h
  }
  
  
  ## Cluster other hidden variables: 
  # new clust membership: use hierarchical clustering membership to replace louvain clustering membership
  clust_all_new <- vector( mode="double",length=length(clust_all) )
  
  for ( cls in 1:length(table(clust_all)) ) {
    clust_all_new[clust_all==cls] <- clust_h[cls]
  }
  
  
  
  # Combine alpha according to the hierarchical clustering result: 
  alpha_hier <- vector( mode="double", length=length(table(g_result$clust_h)) )
  for ( i in 1:length(table(g_result$clust_h))  ) { 
    alpha_hier[i] <- sum(g_result$alpha[clust_h==i])
  }
  
  
  
  # Combine posterior samples of topic weights: 
  tw_hier <- matrix( nrow=D, ncol=length(table(clust_h)), 0 )

  for ( ps in 1:N_ps ) { 
    
    tw_ps <- g_result$Z_sum_samples[[ps]]
    tw_ps_hier <- matrix( nrow=D, ncol=length(table(clust_h)), 0 )
    
    for ( i in 1:length(table(g_result$clust_h))  ) { 
      tw_ps_hier[,i] <- ( rowSums( tw_ps[,which(clust_h==i),drop=FALSE] ) + alpha_hier[i] ) / (S+sum(alpha_hier))
    }
    
    tw_hier <- tw_hier + tw_ps_hier
    
  }
  
  tw_hier <- tw_hier / N_ps
  
  
  
  # return results: 
  g_result$tw_ave <- tw_hier
  g_result$alpha_ave <- alpha_hier
  
  return(g_result)
  
}



