library(FlowSOM)

flowsanitize = function(Y) {
  colnames(Y) = 1:NCOL(Y)
  return(flowCore::flowFrame(Y))
}
  

# Read data.
y = {
  prefix = "img/y"
  suffix = "-prepped.csv"
  lapply(as.list(paste0(prefix, 1:3, suffix)), read.csv)
}

# Indices
N = sapply(y, nrow)
idx_lower = cumsum(c(1, N[-length(N)]))
idx_upper = cumsum(N)
# idx_upper - idx_lower + 1

# Get sample index
I = length(N)
sample_idx = unlist(sapply(1:I, function(i) rep(i, N[i])))

# Stack into one matrix 
Y = Reduce(rbind, y)
J = NCOL(Y)

# Prep data to be in FlowSOM data format
ff_Y = flowsanitize(as.matrix(Y))

# Flowsom clusterings
flowsom.result = FlowSOM(ff_Y,
                         colsToUse=1:J,
                         maxMeta=20, # same as MCMC
                         seed=0)
# Get cluster labels
fsmeta = flowsom.result$meta
fsclus = fsmeta[flowsom.result$FlowSOM$map$mapping[, 1]]

# Print clusterings for each sample
for (i in 1:length(y)) {
  idx_i = idx_lower[i]:idx_upper[i]
  cluster_i = as.numeric(fsclus)[idx_i]
  write.table(cbind(cluster_i, i),
              file=paste0("img/flowsom-", i, ".csv"),
              quote=F, row.names=F, col.names=F)
}

# Print clusterings for all samples
write.table(cbind(fsclus, sample_idx),
            file=paste0("img/flowsom-clusterings.csv"),
            quote=F, row.names=F, col.names=F)
