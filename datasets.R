# dataset generation
source("./functions.R")

## iris ##
ds1 <- iris
colnames(ds1) <- c("X_1", "X_2", "X_3", "X_4", "cluster")
cluster <- as.numeric(ds1$cluster)
ds1 <- as.data.frame(scale(ds1[,c(1:4)]))
ds1$cluster <- cluster
write.csv(ds1, "./iris.csv", row.names = F)

## synthetic ##
for (k in 1:20) {
ds1 <- generate_synthetic_dataset()
write.csv(ds1, paste0("./synthetic", k, ".csv"), row.names = F)
}