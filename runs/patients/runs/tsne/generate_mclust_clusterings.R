library(mclust)

# Read data.
y = {
  prefix = "img/y"
  suffix = "-prepped.csv"
  lapply(as.list(paste0(prefix, 1:3, suffix)), read.csv)
}

# Mclust-selected number of clusters (according to BIC for VVV)
# For VII, the marginal increase in BIC starts decreasing at G=7.
num_clusters = 7

# Indices
N = sapply(y, nrow)
idx_lower = cumsum(c(1, idx[-length(N)]))
idx_upper = cumsum(N)
# idx_upper - idx_lower + 1

# Stack into one matrix 
Y = Reduce(rbind, y)

# Do Mclust on all samples
mclust.result = Mclust(Y, G=num_clusters, model="VII")

for (i in 1:length(y)) {
  write.table(mclust.result$class[idx_lower[i]:idx_upper[i]],
              file=paste0("img/mclust-", i, ".csv"),
              quote=F, row.names=F, col.names=F)
}
