# kde
filenames <- c("./iris.csv","./synthetic1.csv")
for (filename in filenames) {
  # read each file
  ds1 <- read.csv(filename)
  # plot with cluster labels of the dataset
  cluster_colors <- c("red", "blue", "green")
  # for each variable, plot the KDE and colored rug
  for (col_name in colnames(subset(ds1, select = -cluster))) {
    png(filename = paste0(sub(".*/(.*)\\.csv$", "\\1", filename), col_name, ".png"), width = 800, height = 600)
    # density computation
    d <- density(ds1[[col_name]])
    # manual plot
    plot(d$x, d$y, type = "l",
         main = col_name,
         xlab = "", ylab = "",
         col = "black", lwd = 2)
    # rug for each point
    for (i in 1:nrow(ds1)) {
      rug(ds1[[col_name]][i], col = cluster_colors[  ds1[["cluster"]][i] ], lwd = 2)
    }
    dev.off()
  }
}