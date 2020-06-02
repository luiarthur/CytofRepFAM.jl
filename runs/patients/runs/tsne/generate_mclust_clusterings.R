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

# Number of cells per sample
N = sapply(y, nrow)

# Get sample index
I = length(N)
sample_idx = unlist(sapply(1:I, function(i) rep(i, N[i])))

# Stack into one matrix 
Y = Reduce(rbind, y)

# Do Mclust on all samples
mclust.result = Mclust(Y, G=num_clusters, model="VII")

# Write results
write.table(cbind(mclust.result$class, sample_idx),
            file=paste0("img/mclust-clusterings.csv"),
            quote=F, row.names=F, col.names=F)
