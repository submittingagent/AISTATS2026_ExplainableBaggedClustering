# packages and functions

require(dplyr)
require(datasets)
require(MASS)
require(fossil) # Rand
require(dendextend) # FM
require(crank) # permute
require(data.table)
require(clValid) # silhouette
require(MCMCpack) # rdirichlet

generate_synthetic_dataset <- function(prp = c(3,3,3),size=33,
                                       means = rbind(c(0,0,1),c(3,3,3),c(3,6,5)), sg = 1, 
                                       cst = 0, mat = rbind(diag(c(0.5,0.1,0.09)))) {
    dim = ncol(means)
    k = length(prp)
    df <- data.frame()
    for (i in 1:k) {
      df_k <- as.data.frame(mvrnorm(n = size*prp[i], 
                                    mu = means[i,],  
                                    Sigma = (mat*sg) ))
      df_k$cluster <- i
      df <- rbind(df, df_k)
    }
    df_normalized <- df
    cls <- c()
    for (dmn in 1:dim) {cls <- append(cls, paste("X", dmn, sep = "_"))}
    colnames(df_normalized) <- c(cls,"cluster")
    cluster <- df_normalized$cluster
    df_normalized$X_4 <- runif(nrow(df_normalized), min = min(c(df_normalized$X_1,df_normalized$X_2)), max = max(c(df_normalized$X_1,df_normalized$X_2)))
    df_normalized$X_5 <- runif(nrow(df_normalized), min = min(c(df_normalized$X_1)), max = max(c(df_normalized$X_1)))
    df_normalized$X_6 <- runif(nrow(df_normalized), min = min(c(df_normalized$X_2)), max = max(c(df_normalized$X_2)))
    df <- as.data.frame(scale(df_normalized[,c(1,2,3,5,6,7,4)]))
    df$cluster <- cluster
    df
}

clust2 <- function(df, k=3){
  
  df <- as.data.table(df)
  
  km <- kmeans(subset(df, select = -index),k, nstart = 20, iter.max = 15)
  
  df$cluster <- km$cluster
  
  df
  
}

mypermute <- function(dataf){
  
  vec <- as.numeric(dataf$cluster)
  
  df <- numeric()
  
  val <- unique(vec)
  
  perm_tot <- permute(val)
  
  for (i in 1:nrow(perm_tot)) {
    
    perm <- perm_tot[i,]
    
    v <- vec
    
    for (j in 1:ncol(perm_tot)) {
      
      v[which(vec == val[j])] <- perm[j]
      
    }
    
    df <- cbind(df,v)
    
  }
  
  as.data.frame(df)
  
}

clu_labels_adj <- function(df1,df2){
  
  temp <- df2
  temp$original_cluster <- 0
  
  for (i in 1:nrow(temp)) {
    temp$original_cluster[i] <- df1$cluster[which(df1$index == temp$index[i])]
  }
  
  clu_ind <- mypermute(temp)
  
  cnt <- 0
  
  ind <- 1
  
  for (l in 1:ncol(clu_ind)) {
    if (length(which((temp$original_cluster-clu_ind[[l]]) == 0)) > cnt) {
      cnt <- length(which((temp$original_cluster-clu_ind[[l]]) == 0))
      ind <- l
    }
    
  }
  
  df2$cluster <- clu_ind[[ind]]
  
  df2
  
}

fin_res <- function(my_data){
  # Initialize an empty dataframe to store the result
  # result <- data.frame(index = unique(my_data$index))
  result <- data.frame(index = sort(unique(my_data$index[my_data$index != 0])))
  additional_columns <- matrix(0, nrow = nrow(result), ncol = max(my_data$cluster))
  colnames(additional_columns) <- paste0("partecipation_to_cluster_", 1:max(my_data$cluster))
  
  # Concatenate additional columns to my_data
  result <- cbind(result, additional_columns)
  
  # Loop through each unique value of "index"
  for (i in 1:nrow(result)) {
    # Subset the data for the current index value
    subset_data <- my_data[my_data$index == result$index[i], ]
    
    # Count the occurrences of each cluster label
    counts <- table(subset_data$cluster)
    
    counts <- counts[order(as.integer(names(counts)))]
    
    result[i, paste0("partecipation_to_cluster_", names(counts))] <- counts
    
  }
  
  # Replace NA values with 0
  result[is.na(result)] <- 0
  
  result <- result[order(result$index),]
  
  result_fin_temp <- result[,-1]
  
  for (i in 1:nrow(result)) {
    result$cluster[i] <- which(result_fin_temp[i,] == max(result_fin_temp[i,]))[1]
  }
  
  result
  
}

kde_entropy <- function(x, n = 512) {
  x <- na.omit(x)
  dens <- density(x, kernel = "gaussian", n = n)
  p <- dens$y
  p[p <= 0] <- .Machine$double.eps
  dx <- diff(dens$x)[1]
  return(-sum(p * log(p)) * dx) 
}

conditional_entropy <- function(x, cluster_labels) {
  x <- na.omit(data.frame(x = x, cluster = cluster_labels))
  clusters <- unique(x$cluster)
  H_cond <- 0
  
  for (cl in clusters) {
    x_c <- x$x[x$cluster == cl]
    p_c <- length(x_c) / nrow(x)
    H_xc <- kde_entropy(x_c)
    H_cond <- H_cond + p_c * H_xc
  }
  return(H_cond)
}

mutual_information <- function(ds) {
  cluster_lab<-ds$cluster
  ds_clean <- ds[, !(names(ds) %in% "cluster"), drop = FALSE]
  ds_clean<-data.frame(ds_clean)
  MI_values<-numeric()
  for(i in 1:ncol(ds_clean)){
    H <- kde_entropy(ds_clean[,i])
    H_given_C <- conditional_entropy(ds_clean[,i], cluster_lab)
    MI_values <- c(MI_values,H - H_given_C)
  }
  names(MI_values)<-names(ds_clean)
  return(MI_values)
}

explainable_clustering <- function(data, ground_truth = NA, nvar = ncol(data)-1, ncl, B = 100) {
  # Dataset must have an "index" column as the first column
  
  cl0 <- clust2(data, ncl)  # initial k-means for label alignment
  
  res_final <- data.frame()
  dunn_ind <- numeric()
  rand_ind <- numeric()
  FM_ind <- numeric()
  sil_ind <- numeric()
  Mut_info <- numeric()
  var_selected <- list()
  
  for (i in 1:B) {
    print(i)
    n <- nrow(data)
    
    # bootstrap resampling (currently just copy)
    bdts <- data
    index <- bdts$index
    
    # select random subset of variables (excluding index)
    col_ind <- sample(which(!names(bdts) %in% c("index")), nvar)
    df <- bdts[, sort(col_ind), drop = FALSE]
    
    # k-means clustering
    km <- kmeans(df, centers = ncl, iter.max = 50)
    df$cluster <- km$cluster
    bdts$cluster <- km$cluster
    df$index <- index
    
    # adjust labels to match initial k-means
    df <- clu_labels_adj(cl0, df)
    bdts$cluster <- km$cluster
    bdts <- clu_labels_adj(cl0, bdts)
    
    # --- Dunn index (base R subsetting) ---
    dunn_ind <- c(dunn_ind, clValid::dunn(
      clusters = df$cluster,
      Data = df[, !(names(df) %in% c("index", "cluster")), drop = FALSE],
      method = "euclidean"
    ))
    
    # Silhouette index
    sil <- silhouette(df$cluster, dist(as.matrix(df[, !(names(df) %in% c("index", "cluster"))])))
    sil_ind <- c(sil_ind, mean(sil[, "sil_width"]))
    
    # Rand index
    if (length(ground_truth) == nrow(df)) {
      rand_ind <- c(rand_ind, rand.index(df$cluster, ground_truth))
      FM_ind   <- c(FM_ind, FM_index(df$cluster, ground_truth))
    } else {
      rand_ind <- NA
      FM_ind   <- NA
    }
    
    # Mutual information (drop index)
    Mut_info <- rbind(Mut_info, mutual_information(bdts[, !(names(bdts) %in% "index"), drop = FALSE]))
    
    var_selected[[i]] <- names(df)[1:length(col_ind)]
    
    # store results
    final <- bdts
    final$cluster <- bdts$cluster
    final$copy <- i
    res_final <- rbind(res_final, final)
  }
  
  Mut_info <- data.frame(Mut_info)
  df_fin <- fin_res(res_final)
  
  # final normalization (drop index and cluster safely)
  res_df_fin <- df_fin[, !(names(df_fin) %in% c("index", "cluster")), drop = FALSE]
  res_df_fin <- res_df_fin / rowSums(res_df_fin)
  dati_final <- data.frame(index = df_fin$index, res_df_fin, cluster = df_fin$cluster)
  
  return(list(
    final_results = dati_final,
    dunn_index = dunn_ind,
    sil_index = sil_ind,
    rand_index = rand_ind,
    FM_index = FM_ind,
    MI = Mut_info,
    variable_sel = var_selected
  ))
}

explainable_efron_clustering <-function(data, ground_truth=NA, nvar=ncol(data)-1,ncl,B=100){
  #Dataset deve avere come prima colonna una colonna chiamata "index"
  
  cl0 <- clust2(data,ncl) #kmean
  # table(cl0$cluster, ground_truth)
  
  res_final <- data.frame()
  dunn_ind <- numeric()
  rand_ind <- numeric()
  FM_ind <- numeric()
  sil_ind <- numeric()
  
  Mut_info <- numeric()
  var_selected<-list()
  
  for (i in 1:B) {
    print(i)
    n <- nrow(data)
    #### bootstrap resampling if implemented
    bdts <- data[sample.int(nrow(data), nrow(data), replace = TRUE),]
    ####
    index <- bdts$index
    col_ind<-sample(which(!names(bdts)%in%c("index")),nvar)
    df <- bdts[,sort(col_ind)]
    km <- kmeans(df, ncl, iter.max = 50)
    df$cluster <- km$cluster
    df$index <- index
    df <- clu_labels_adj(cl0,df)
    bdts$cluster <- km$cluster
    bdts <- clu_labels_adj(cl0,bdts)
    # Dunn index su tutte le variabili
    dunn_ind<-c(dunn_ind,clValid::dunn(clusters=df$cluster,Data=df%>%dplyr::select(-index,-cluster),method = "euclidean"))
    # Silhouette index su tutte le variabili
    sil <- silhouette(df$cluster, dist(as.matrix(df)))
    sil_ind <- c(sil_ind,mean(sil[, "sil_width"]))
    # rand ind, with the sorted truelabels
    if (length(ground_truth)==nrow(df)){rand_ind <- c(rand_ind,rand.index(df$cluster,ground_truth[index]))}
    else {rand_ind <- NA}
    if (length(ground_truth)==nrow(df)){FM_ind <- c(FM_ind,FM_index(df$cluster,ground_truth[index]))}
    else {FM_ind <- NA}
    # mutual information, tolgo la colonna index
    Mut_info<-rbind(Mut_info,mutual_information(bdts%>%dplyr::select(-index)))
    var_selected[[i]] <- names(df)[1:length(col_ind)]
    
    final<-bdts
    final$cluster <- bdts$cluster
    final$copy <- i
    res_final <- rbind(res_final, final)
  }
  
  Mut_info <- data.frame(Mut_info)
  df_fin <- fin_res(res_final)
  res_df_fin <- subset(df_fin, select = -c(index, cluster))
  res_df_fin <- res_df_fin/rowSums(res_df_fin)
  dati_final<-data.frame(index=df_fin$index,res_df_fin,cluster=df_fin$cluster)
  
  return(list(final_results=dati_final,dunn_index=dunn_ind,sil_index=sil_ind, rand_index=rand_ind,FM_index=FM_ind,MI=Mut_info,variable_sel=var_selected))
}

getcenters <- function(dtst) {
  centers <- data.frame()
  dtst <- as.data.frame(dtst)
  for (i in 1:length(unique(dtst$cluster))) {
    centers <- rbind(centers, colMeans(dtst[dtst$cluster == i, !(names(dtst) %in% c("index", "cluster"))]))
  }
  centers
}

getcovariances <- function(dtst) {
  dtst <- as.data.frame(dtst)
  cls <- length(unique(dtst$cluster))
  
  covariances <- list()
  
  for (i in 1:cls) {
    
    # Store the covariance matrix in the list
    covariances[[i]] <- cov(dtst[dtst$cluster == i, !(names(dtst) %in% c("index", "cluster"))])
    
  }
  
  covariances
  
}

explainable_BBC <-function(data, ground_truth=NA, nvar=ncol(data)-1,ncl,wg=0.1,l=1,B=100){
  #Dataset deve avere come prima colonna una colonna chiamata "index" di interi
  
  cl0 <- clust2(data,ncl) #kmean
  cl_ds <- data.frame()
  # table(cl0$cluster, ground_truth)
  
  for (i in 1:20){ #stima su 20 tentativi pr trovare i centroidi e proporzioni stabili
    cl1 <- clust2(data,ncl)
    cl1 <- clu_labels_adj(cl0,cl1)
    cl_ds <- rbind(cl_ds, cl1)
  }
  
  centers <- getcenters(cl_ds)
  names(centers)<-names(cl_ds%>%dplyr::select(-index,-cluster))
  
  cluster_proportions <- sapply(1:ncl, function(i) {
    length(which(cl_ds$cluster == i)) / nrow(cl_ds)
  })
  
  sig <- getcovariances(cl0) # evaluated on the first cluster
  
  res_final <- data.frame()
  dunn_ind<-numeric()
  rand_ind<-numeric()
  FM_ind<-numeric()
  sil_ind<-numeric()
  
  Mut_info<-numeric()
  var_selected<-list()
  
  for (i in 1:B) {
    
    print(i)
    ds_clusters<-cl0
    w=wg
    mu=centers
    sigma=sig
    dts <- as.data.table(subset(ds_clusters, select = -c(cluster)))
    # Calculate the number of clusters
    num_clusters <- nrow(mu)
    
    # Parameters bayesian resampling  
    n <- nrow(dts)
    k <- n * w / (1 - w)
    m_ <- n
    
    weights <- rdirichlet(1, rep((n + k) / m_, m_))
    
    xstar <- data.table()
    index <- numeric()
    
    for (j in 1:m_) {
      soglia <- runif(1, 0, 1)
      
      if (soglia > w) {
        ind <- sample.int(n, size = 1)
        xstar <- rbind(xstar, dts[ind, ])
      } else {
        soglia2 <- runif(1, 0, 1)
        cumulative_proportion <- 0
        
        for (p in 1:num_clusters) {
          cumulative_proportion <- cumulative_proportion + cluster_proportions[p]
          if (soglia2 < cumulative_proportion) {
            media <- mu[p,]
            sigma_cl <- sigma[[p]]
            #break
          }
        }
        
        v_0 <- data.table(0,t(mvrnorm(1, as.numeric(media), sigma_cl)))
        names(v_0)<-colnames(dts)
        xstar <- rbind(xstar, v_0)
      }
    }
    
    resample_cumulative <- cumsum(weights)
    rnd <- runif(m_)
    
    ds <- data.table()
    
    for (boh in 1:m_) {
      
      if(rnd[boh]<min(resample_cumulative)){
        ds <- rbind(ds,xstar[1,])
      }
      else {
        ds <- rbind(ds,xstar[max(which(resample_cumulative <= rnd[boh]))+1,])
      }
      
    }
    bsdts_mixture<-ds
    
    df <- as.data.table(bsdts_mixture)
    col_ind<-sample(which(!names(df)%in%c("index","cluster")),nvar)
    df_cluster<-df %>% dplyr::select(index,names(df)[sort(col_ind)])
    centers_i <- centers %>% dplyr::select(names(df)[sort(col_ind)])
    
    # nstart per evitare cluster vuoti
    km <- kmeans(subset(df_cluster, select = -index), centers_i, nstart = 20, iter.max = 50)
    
    df$cluster <- km$cluster
    #Dunn index su tutte le variabili
    dunn_ind<-c(dunn_ind,clValid::dunn(clusters=df$cluster,Data=df%>%dplyr::select(-index,-cluster),method = "euclidean"))
    # Silhouette index su tutte le variabili
    sil <- silhouette(df$cluster, dist(as.matrix(df)))
    sil_ind <- c(sil_ind,mean(sil[, "sil_width"]))
    # rand ind, with the sorted truelabels
    if (length(ground_truth)==nrow(df)){rand_ind <- c(rand_ind,rand.index(df$cluster[which(bsdts_mixture$index !=0)], ground_truth[bsdts_mixture$index[which(bsdts_mixture$index !=0)]]))}
    else {rand_ind <- NA}
    if (length(ground_truth)==nrow(df)){FM_ind <- c(FM_ind,FM_index(df$cluster[which(bsdts_mixture$index !=0)], ground_truth[bsdts_mixture$index[which(bsdts_mixture$index !=0)]]))}
    else {FM_ind <- NA}
    df <- as.data.frame(df)
    Mut_info<-rbind(Mut_info,mutual_information(df%>%dplyr::select(-index)))
    var_selected[[i]]<-names(df)[col_ind]
    
    final<-df
    final$copy <- i
    
    res_final <- rbind(res_final, final)
    
  }
  
  Mut_info<-data.frame(Mut_info)
  df_fin <- fin_res(res_final)
  res_df_fin <- subset(df_fin, select = -c(index, cluster))
  res_df_fin <- res_df_fin/rowSums(res_df_fin)
  dati_final<-data.frame(index=df_fin$index,res_df_fin,cluster=df_fin$cluster)
  
  return(list(final_results=dati_final,dunn_index=dunn_ind,sil_index=sil_ind,rand_index=rand_ind,FM_index=FM_ind,MI=Mut_info,variable_sel=var_selected))
}

MI_read <- function(data_results) {
  MI_selected<-data_results$MI
  for(i in 1:nrow(MI_selected)){
    MI_selected[i,!names(MI_selected)%in%data_results$variable_sel[[i]]]<-NA
  }
  MI_ <-colMeans(MI_selected,na.rm=T)
  MI_ <- MI_/max(MI_)
  MI_weighted <- colMeans(MI_selected*data_results$dunn_ind*100, na.rm=T)
  MI_weighted <- MI_weighted/max(MI_weighted)
  return(list(MI_, MI_weighted))
}

format_latex_mean_sd <- function(mean_mat, sd_mat, digits = 2) {
  # Ensure dimensions match
  stopifnot(all(dim(mean_mat) == dim(sd_mat)))
  
  n_rows <- nrow(mean_mat)
  n_cols <- ncol(mean_mat)
  
  # Start LaTeX table
  latex_output <- "\\begin{tabular}{|"
  latex_output <- paste0(latex_output, paste(rep("c|", n_cols), collapse = ""), "}\n\\hline\n")
  
  for (i in 1:n_rows) {
    row_entries <- sapply(1:n_cols, function(j) {
      mean_val <- round(mean_mat[i, j], digits)
      sd_val <- round(sd_mat[i, j], digits)
      paste0(mean_val, " $\\pm$ ", sd_val)
    })
    latex_output <- paste0(latex_output, paste(row_entries, collapse = " & "), " \\\\ \\hline\n")
  }
  
  latex_output <- paste0(latex_output, "\\end{tabular}")
  
  return(latex_output)
}