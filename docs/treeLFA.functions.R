### Training with the gibbs-em algorithm: 
gibbs_EM_train <- function( topic.number,
                            data, ...,
                            Phi, I, rho, alpha, Z,
                            burn_in,
                            opt_N_1, cycle_1,
                            opt_N_2, cycle_2 ) {
  
  
  # Prepare the input data:
  S <- ncol(data)     # number of terminal disease codes on the tree
  
  
  # check if the tree str is provided:
  dots <- list(...)
  if (!"tree_str" %in% names(dots)) {
    message("tree structure is not provided, use a non-informative tree structure as default") 
    tree_str <- data.table( matrix( nrow=ncol(data), ncol=2, "") )
    S1 <- 0
    colnames(tree_str) <- c("node","parent")
    tree_str$parent <- rep("root",nrow(tree_str))
    tree_str$node <- colnames(data)
    print(tree_str)
  } else {
    tree_str <- dots$`tree_str`
    S1 <- nrow(tree_str) - ncol(data)    # number of internal disease codes on the tree
  }
  

  data <- as.matrix(data)
  colnames(tree_str) <- c("node","parent")

  
  
  
  
  # Check the input diagnosis data: 
  # (1) the matrix is a double matrix
  if ( !(is.integer(data)) & !(is.double(data)) ) { 
    stop("The input matrix is not a interger matrix!")  
  }
  # (2) check the matrix only contains 0/1
  check01_fun <- function(x) all(x %in% c(0, 1))
  if ( !check01_fun(data) ){ 
    stop("The input matrix contains values other than 0 and 1!")  
  }
  
  
  
  
  
  ### Tree structure: 
  # Check the input tree str matrix: 
  nodes_ch <- colnames(data)
  nodes_pa <- tree_str$node[which(!(tree_str$node %in% colnames(data)))]
  nodes_ch <- tree_str$node[which(!(tree_str$node %in% nodes_pa))]
  
  # (1) Check the tree str matrix contains the same terminal disease codes as the input data matrix: 
  if ( !setequal(colnames(data),c(nodes_ch)) ) { 
    stop("The input diagnosis data matrix and the tree structure matrix do not have the same set of disease codes")  
  }
  
  # (2) Check that there is no repetitive codes in the 1st column of the tree str matrix: 
  if ( sum(duplicated(tree_str$node))>0 ) { 
    stop("There is repetitive nodes in the tree structure matrix!")
  }
  
  
  # (3)-1 Check the 1st column in the tree str matrix doesn't contain the "root" node
  if( "root" %in% tree_str$node ) { 
    stop("The root node cannot be in the 1st column of the tree structure matrix!")
  }
  
  # (3)-2 Check the 2nd column in the tree str matrix contains the "root" node
  if ( !("root" %in% tree_str$parent) ) { 
    stop("The root node is not in the 2nd column of the tree structure matrix!")
  }
  
  # (4) Check the parent nodes in the 2nd column are not the terminal codes in the 1st column: 
  if ( sum( unique(tree_str$parent) %in% nodes_ch ) > 0  ) { 
    stop("The parent nodes in the 2nd column cannot be the terminal disease codes in the 1st column")  
  }
  
  # (5) Check the parent codes in the 2nd column are the same as the parent codes in the 1st column: 
  # which means all the parent codes in the 2nd column are also in the 1st column 
  if ( !setequal( unique(tree_str$parent), c("root",nodes_pa) ) ) { 
    stop("The total number of unique parent nodes in the 2nd column \ndoes not equal the total number of unique parent nodes in the 1st column")  
  }
    
  
  
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
  # Check the length of alpha is the same as the total number of topics: 
  if ( length(alpha) != topic.number ) { 
    stop("The length of alpha should be the same of the number of topics to be inferred!")
  }

  
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
gibbs_train <- function( data, ...,
                         Phi, I, rho, alpha, Z,
                         burn_in, cycle, interval ) {
  
  # Prepare the input data:
  topic.number <- nrow(Phi)
  S <- ncol(data)    # number of terminal disease codes on the tree 
  
  # check if the tree str is provided:
  dots <- list(...)
  if (!"tree_str" %in% names(dots)) {
    message("tree structure is not provided, use a non-informative tree structure as default") 
    tree_str <- data.table( matrix( nrow=ncol(data), ncol=2, "") )
    S1 <- 0
    colnames(tree_str) <- c("node","parent")
    tree_str$parent <- rep("root",nrow(tree_str))
    tree_str$node <- colnames(data)
    print(tree_str)
  } else {
    tree_str <- dots$`tree_str`
    S1 <- nrow(tree_str) - ncol(data)    # number of internal disease codes on the tree
  }

  
  S1 <- nrow(tree_str) - ncol(data)      # number of internal disease codes on the tree 
  # S1 <- 5
  data <- as.matrix(data)
  
  colnames(tree_str) <- c("node","parent")
  
  
  # ## Process the tree structure: 
  # # Terminal codes on the tree:
  # tree_ter <- tree_str[terminal=="Y",]
  # 
  # ## Reorder terminal/internal codes: 
  # # internal codes:
  # tree_int <- tree_str[!(node %in% colnames(data)),]
  # 
  # # Reorder internal codes: according to layers of codes 
  # tree_int<- tree_int[order(layer),]
  # 
  # # Reorder termianl codes: according to the columns (order of disease codes) of the input data
  # tree_ter <- tree_ter[match(colnames(data),tree_ter$node)]
  # 
  # tree_str <- rbind(tree_int,tree_ter)
  
  
  
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
  b10 <- 3
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
  
  # if ( ((cycle-burn_in) %% interval) != 0 ) { 
  #  stop("(\"cycle\"-\"burn_in\") should be devisible by \"interval\"!")
  # }
  
  
  
  
  
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


# plot distributions of ess: 
ess_plot <- function( ess, bin_N, N_ps ) { 
  
  ess_plot <- as.vector(ess)
  ess_plot_large <- ess_plot[which(ess_plot>N_ps)]
  ess_plot_small <- ess_plot[which(ess_plot<=N_ps)]
  ess_plot_small_cut <- cut(ess_plot_small,breaks=bin_N)
  
  
  df_small <- data.table( cbind( ess_plot_small, as.character(ess_plot_small_cut) ) )
  colnames(df_small) <- c("value","category")
  df_small$category <- factor(df_small$category,levels=levels(ess_plot_small_cut))
  
  df_large <- data.table( cbind(ess_plot_large,paste(">",N_ps,sep="")) )
  colnames(df_large) <- c("value","category")
  
  df_plot <- data.table( rbind(df_small,df_large) )
  colnames(df_plot) <- c("value","category")

  ggplot(df_plot, aes(x=category)) +
    geom_bar(fill="steelblue", color="white") + 
    labs(x="Range of effective sample size", y="Count of variables")

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
    dist_topics <- 1 - cosine(t(topics_ave))
  dist_topics <- dist(dist_topics)
  
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










### Cluster posterior samples of topics using Louvain algorithm:
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
  
  
  
  ## Store and combine posterior samples of topic weights: 
  # averaged topic weights for individuals and store topic weights for individuals: 
  tw_hier <- matrix( nrow=D, ncol=length(table(clust_h)), 0 )

  tw_all <- vector( "list", length=length(table(clust_h)) )
  for ( j in 1:length(table(clust_h)) ){ 
    tw_all[[j]] <- matrix(nrow=D,ncol=N_ps)
  }
  
  topics_all <- vector( "list", length=length(table(clust_h)) )
  for ( j in 1:length(table(clust_h)) ){ 
    topics_all[[j]] <- matrix(nrow=S,ncol=N_ps)
  }
  
  
  for ( ps in 1:N_ps ) { 
    
    tw_ps <- g_result$Z_sum_samples[[ps]]
    tw_ps_hier <- matrix( nrow=D, ncol=length(table(clust_h)), 0 )
    
    topics_ps <- g_result$Phi_samples[[ps]]
    
    for ( i in 1:length(table(g_result$clust_h))  ) { 
      ps_range <- c((ps-1)*10+1):c((ps-1)*10+10)
      clust_ps <- clust_all_new[ps_range]
      tw_ps_hier[,i] <- ( rowSums( tw_ps[,which(clust_ps==i),drop=FALSE] ) + alpha_hier[i] ) / (S+sum(alpha_hier))
    }
    
    tw_hier <- tw_hier + tw_ps_hier
    
    # collect posterior samples of averaged topic weights: 
    for ( j in 1:length(table(clust_h)) ){ 
      tw_all[[j]][,ps] <- tw_ps_hier[,j]
      topics_all[[j]][,ps] <- colMeans( topics_ps[which(clust_ps==j),,drop=FALSE] )
    }
  
  }
  
  tw_hier <- tw_hier / N_ps
  
  
  # calculate ess for all tw: D*K
  print("calculate effective sample size for topic weights")
  ess_tw <- matrix(nrow=D,ncol=length(table(clust_h)))
  for ( k in 1:length(table(clust_h)) ) {
    tw_ess <- tw_all[[k]]
    tw_ess <- mcmc(tw_ess)
    for ( d in 1:D )  { 
      ess_tw[d,k] <- effectiveSize(tw_ess[d,]) 
    }
  }
  
  # calculate ess for all topics: K*S
  print("calculate effective sample size for topics")
  ess_topics <- matrix(nrow=length(table(clust_h)),ncol=S)
  for ( k in 1:length(table(clust_h)) ) {
    topics_ess <- topics_all[[k]]
    topics_ess <- mcmc(topics_ess)
    for ( s in 1:S )  { 
      ess_topics[k,s] <- effectiveSize(topics_ess[s,]) 
    }
  }
  
  
  # return results: 
  g_result$tw_ave <- tw_hier
  g_result$alpha_ave <- alpha_hier
  g_result$ess_tw <- ess_tw
  g_result$ess_topics <- ess_topics
  g_result$topics_all <- topics_all
  g_result$tw_all <- tw_all
  
  return(g_result)
  
}





### evaluation of topics: 
## topic coherence: 
topics_coherence <- function(topics, data) { 
  
  TC <- 0 
  
  for ( k in 1:nrow(topics) ) { 
    
    TC_topic <- 0 
    
    for ( s1 in 1:(ncol(topics)-1) ) { 
      
      for ( s2 in (s1+1):ncol(topics) ) { 
        
        P1 <- sum(data[,s1])/nrow(data)
        P2 <- sum(data[,s2])/nrow(data)
        P12 <- sum( data[,s1] * data[,s2] ) / nrow(data)
        f_tmp <- ( log( P12/(P1*P2) ) / -log(P12) )
        
        TC_topic <- TC_topic + f_tmp
      }
      
    }
    
    TC <- TC + TC_topic/( ncol(data)*(ncol(data)-1) )
    
  }
  
  TC <- TC/ncol(data)
  
  return(TC)
  
} 





## topic diversity: 
topics_diversity <- function(topics, data, N) { 
  
  TD <- 0 
  
  codes_top <- vector(length=0,mode="character")  
  
  for ( k in 1:nrow(topics) ) { 
    codes_top <- c(codes_top,colnames(data)[order( topics[k,], decreasing=TRUE )[1:min(N,ncol(topics))]])
  }
  codes_top_unique <- unique(codes_top)
  
  TD <- codes_top_unique/(N*nrow(topics))
  
  return(TD)
  
} 
