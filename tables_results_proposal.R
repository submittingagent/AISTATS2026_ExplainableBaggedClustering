source("./functions.R")
## select the dataset type
dataset <- "synthetic"
# dataset <- "iris"

rand_matrix_list <- list()
FM_matrix_list <- list()
MI_matrix_selected_list <- list()
kmeans_rand <- numeric()
h_rand <- numeric()
kmeans_FM <- numeric()
h_FM <- numeric()
res_basic_list <- list()
res_efron_list <- list()
res_BBC_list <- list()
  
## proposal: computation of results
  for (n in 1:20) {
  # dataset read
  if (dataset == "iris") {ds1 <- read.csv(paste0("./", dataset, ".csv"))}
  if (dataset == "synthetic") {ds1 <- read.csv(paste0("./", dataset, n, ".csv"))}
  
  sample_size = nrow(ds1)
  n_cluster = length(unique(ds1$cluster))
  trueLabels <- ds1$cluster
  ds1 <- subset(ds1, select = -cluster)
  total_features = ncol(ds1)
  
  # feature dropout selection: m number of features involved in clustering
  m = ncol(ds1)-1
  
  # K-means
  cl2 <- ds1
  cl2$cluster <- kmeans(ds1, n_cluster)$cluster
  kmeans_rand = append(kmeans_rand, rand.index(cl2$cluster, trueLabels)) 
  kmeans_FM = append(kmeans_FM, FM_index(cl2$cluster, trueLabels)) 
  cl2$index <- 1:nrow(cl2)
  
  # hierarchical
  cl3 <- ds1
  dist_matrix <- dist(cl3, method = "euclidean")
  hc <- hclust(dist_matrix, method = "complete")
  cl3$cluster <- cutree(hc, k = n_cluster)
  cl3$index <- 1:nrow(cl3)
  # permutation of labels for fair assessment of the validation indices
  cl3 <- clu_labels_adj(cl2,cl3)
  h_rand <- append(h_rand, rand.index(cl3$cluster, trueLabels))
  h_FM <- append(h_FM, FM_index(cl3$cluster, trueLabels))
  
  # proposal basic
  res_basic <- explainable_clustering(data=data.frame(index=1:nrow(ds1),ds1), ground_truth = trueLabels, nvar=m,ncl=n_cluster,B=100)
  res_basic$rand_index_aggregated <- rand.index(res_basic$final_results$cluster,trueLabels)
  res_basic$FM_index_aggregated <- FM_index(res_basic$final_results$cluster,trueLabels)
  
  # proposal Efron
  res_efron <- explainable_efron_clustering(data=data.frame(index=1:nrow(ds1),ds1), ground_truth = trueLabels, nvar=m, ncl=n_cluster, B=100)
  res_efron$rand_index_aggregated <- rand.index(res_efron$final_results$cluster,trueLabels)
  res_efron$FM_index_aggregated <- FM_index(res_efron$final_results$cluster,trueLabels)
  
  # proposal BBC
  res_BBC <- explainable_BBC(data=data.frame(index=1:nrow(ds1),ds1), ground_truth = trueLabels, nvar=m, ncl=n_cluster, wg=0.3, l=1, B=100)
  res_BBC$rand_index_aggregated <- rand.index(res_BBC$final_results$cluster,trueLabels)
  res_BBC$FM_index_aggregated <- FM_index(res_BBC$final_results$cluster,trueLabels)
  
  res_basic_list[[n]] <- res_basic
  res_efron_list[[n]] <- res_efron
  res_BBC_list[[n]] <- res_BBC
  
  #### mean on the bootstrap copies of rand, FM index
  rand_matrix <- matrix(c(mean(res_basic$rand_index), mean(res_efron$rand_index), mean(res_BBC$rand_index)),
                        nrow = 1,
                        dimnames = list(NULL, c("Standard", "Efron", "BBC")))
  
  rand_matrix_list[[n]] <- rand_matrix
  
  FM_matrix <- matrix(c(mean(res_basic$FM_index), mean(res_efron$FM_index), mean(res_BBC$FM_index)),
                      nrow = 1,
                      dimnames = list(NULL, c("Standard", "Efron", "BBC")))
  
  FM_matrix_list[[n]] <- FM_matrix
  
  mutual_info_basic <- MI_read(res_basic)
  mutual_info_efron <- MI_read(res_efron)
  mutual_info_BBC <- MI_read(res_BBC)
  
  MI_matrix_selected <- rbind(mutual_info_basic[[1]], mutual_info_efron[[1]], mutual_info_BBC[[1]], mutual_info_basic[[2]], mutual_info_efron[[2]], mutual_info_BBC[[2]])
  rownames(MI_matrix_selected) = c("MI_basic", "MI_efron", "MI_BBC", "MI_basic_weighted", "MI_efron_weighted", "MI_BBC_weighted")
  MI_matrix_selected_list[[n]] <- MI_matrix_selected
  
}

## Table generation: Table 2, clustering results
methods <- c("Kmeans", "Hierarchical", "Basic", "Efron", "BBC")
# Rand
vals_basic <- sapply(res_basic_list, function(x) x$rand_index_aggregated)
vals_efron <- sapply(res_efron_list, function(x) x$rand_index_aggregated)
vals_BBC   <- sapply(res_BBC_list,   function(x) x$rand_index_aggregated)
mean_kmeans = mean(kmeans_rand)
sd_kmeans = sd(kmeans_rand)
mean_h = mean(h_rand)
sd_h = sd(h_rand)
mean_basic <- mean(vals_basic)
sd_basic   <- sd(vals_basic)
mean_efron <- mean(vals_efron)
sd_efron   <- sd(vals_efron)
mean_BBC <- mean(vals_BBC)
sd_BBC   <- sd(vals_BBC)
means <- c(mean_kmeans, mean_h, mean_basic, mean_efron, mean_BBC)
sds   <- c(sd_kmeans,   sd_h,   sd_basic,   sd_efron,   sd_BBC)
table_Rand <- matrix(
  mapply(function(m, s) paste0(round(m, 2), " ± ", round(s, 2)),
         m = means,
         s = sds),
  nrow = 1,
  byrow = TRUE
)

# FM
vals_basic <- sapply(res_basic_list, function(x) x$FM_index_aggregated)
vals_efron <- sapply(res_efron_list, function(x) x$FM_index_aggregated)
vals_BBC   <- sapply(res_BBC_list,   function(x) x$FM_index_aggregated)
mean_kmeans = mean(kmeans_FM)
sd_kmeans = sd(kmeans_FM)
mean_h = mean(h_FM)
sd_h = sd(h_FM)
mean_basic <- mean(vals_basic)
sd_basic   <- sd(vals_basic)
mean_efron <- mean(vals_efron)
sd_efron   <- sd(vals_efron)
mean_BBC <- mean(vals_BBC)
sd_BBC   <- sd(vals_BBC)
means <- c(mean_kmeans, mean_h, mean_basic, mean_efron, mean_BBC)
sds   <- c(sd_kmeans,   sd_h,   sd_basic,   sd_efron,   sd_BBC)
table_FM <- matrix(
  mapply(function(m, s) paste0(round(m, 2), " ± ", round(s, 2)),
         m = means,
         s = sds),
  nrow = 1,
  byrow = TRUE
)

table_clustering_results <- rbind(table_Rand,table_FM)
colnames(table_clustering_results) <- methods
row.names(table_clustering_results) <- c("Rand", "FM")
write.table(table_clustering_results, "./table_clustering_results.txt", sep = "\t", quote = FALSE, fileEncoding = "UTF-8")

## Table generation: Table 3 and 5, MI results
mean_MI_matrix_selected <- Reduce("+", MI_matrix_selected_list) / length(MI_matrix_selected_list)
sd_MI_matrix_selected <- sqrt(Reduce("+", lapply(MI_matrix_selected_list, function(mat) {
  (mat - mean_MI_matrix_selected)^2})) / length(MI_matrix_selected_list))

table_unweighted <- matrix(
  mapply(function(m, s) paste0(round(m, 2), " ± ", round(s, 2)),
         m = as.vector(mean_MI_matrix_selected[c(1,2,3), ]),
         s = as.vector(sd_MI_matrix_selected[c(1,2,3), ])),
  nrow = nrow(mean_MI_matrix_selected[c(1,2,3), ]),
  byrow = F
)

table_weighted <- matrix(
  mapply(function(m, s) paste0(round(m, 2), " ± ", round(s, 2)),
         m = as.vector(mean_MI_matrix_selected[c(4,5,6), ]),
         s = as.vector(sd_MI_matrix_selected[c(4,5,6), ])),
  nrow = nrow(mean_MI_matrix_selected[c(1,2,3), ]),
  byrow = F
)

colnames(table_unweighted) <- paste0("X_", seq(1:total_features))
row.names(table_unweighted) <- c("Basic MI", "Efron MI", "BBC MI")
print(table_unweighted)
write.table(table_unweighted, "./table_MI.txt", sep = "\t", quote = FALSE, fileEncoding = "UTF-8")
colnames(table_weighted) <- paste0("X_", seq(1:total_features))
row.names(table_weighted) <- c("Basic WMI", "Efron WMI", "BBC WMI")
print(table_weighted)
write.table(table_weighted, "./table_WMI.txt", sep = "\t", quote = FALSE, fileEncoding = "UTF-8")

## comparisons with literature
source("packages_literaturecomparison.R")
result_scores = NULL

X <- ds1
# URF_gini
# Step 1: Generate synthetic data from marginals
X_synth <- apply(X, 2, sample)
# Step 2: Combine real + synthetic
X_combined <- rbind(X, X_synth)
y_combined <- factor(c(rep(1, nrow(X)), rep(0, nrow(X_synth))))
# Step 3: Train RF classifier
rf_unsup <- randomForest(x = X_combined, y = y_combined, ntree = 500,
  nodesize = 1, mtry = 2, proximity = TRUE, oob.prox = FALSE)
URF_gini = importance(rf_unsup, type = 2)
URF_gini = URF_gini/max(URF_gini)  
print(URF_gini)
result_scores = rbind(result_scores,t(as.matrix(round(URF_gini,2))))

# URF_perm
base_pred <- predict(rf_unsup, type = "response")
base_acc <- mean(base_pred == y_combined)
perm_import <- numeric(ncol(X))
for (j in 1:ncol(X)) {
  X_perm <- X_combined
  X_perm[, j] <- sample(X_perm[, j])  # permute column j
  pred <- predict(rf_unsup, newdata = X_perm, type = "response")
  perm_import[j] <- base_acc - mean(pred == y_combined)
}
names(perm_import) <- colnames(X)
perm_import = perm_import/max(perm_import)
result_scores = rbind(result_scores,t(as.matrix(round(perm_import, 2))))

############ CUBT 
X <- as.matrix(X)
tree_cubt = cubt(X)
cubt_scores = variable.importance(tree_cubt, dat = X)
cubt_scores = cubt_scores/max(cubt_scores)
names(cubt_scores) = names(result_scores)
result_scores = rbind(result_scores,t(as.matrix(round(cubt_scores, 2))))

########### TWKM
# The parameters are (representing penalty terms in the two-level entropy weighted KMeans):
# the group subdivision of features, set as null to make single feature grouping
# lambda: parameter for individual feature weights
# eta: parameter for group weights

TWKM_matscores <- list()
for (i in 1:20) {
  TWKM_clustering = fgkm(X, centers = n_cluster, groups = rep(1, ncol(X)), 
                         lambda = 1, eta = 1, maxiter=100, delta=0.000001,
                         maxrestart=10,seed=-1)
  TWKM_weigths = TWKM_clustering$featureWeight
  # the result is given as a matrix, ranking the importance of each feature for each cluster. Max is taken
  max_scores = apply(as.data.frame(TWKM_weigths), 2, max)
  TWKM_matscores[[i]] = max_scores/max(max_scores)
}

mean_TWKM <- Reduce("+", TWKM_matscores) / length(TWKM_matscores)
sd_TWKM <- sqrt(Reduce("+", lapply(TWKM_matscores, function(mat) {
  (mat - mean_TWKM)^2})) / length(TWKM_matscores))

results_TWKM <- matrix(
  mapply(function(m, s) paste0(round(m, 2), " ± ", round(s, 2)),
         m = as.vector(mean_TWKM),
         s = as.vector(sd_TWKM)),
  nrow = 1,
  byrow = F
)
result_scores = rbind(result_scores,(as.matrix(results_TWKM)))


# LS for all features
X = as.matrix(X)
# proportion: fraction of the features selected
out = do.lscore(X, ndim = 4, type = c("proportion", 1), 
                preprocess = c("null"), t = 10)
ls_scores = out$lscore
ls_scores = ls_scores/max(ls_scores)
names(ls_scores) = names(result_scores)
print(ls_scores) 
result_scores = rbind(result_scores,t(as.matrix(round(ls_scores, 2))))

## Table generation: Table 4 and 6, literature methods
row.names(result_scores) <- c("URF_gini", "URF_perm", "CUBT", "TWKM", "LS")
print(result_scores)
write.table(result_scores, "./table_comparison.txt", sep = "\t", quote = FALSE, fileEncoding = "UTF-8")