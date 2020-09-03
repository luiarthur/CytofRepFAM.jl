library(FlowSOM)

flowsanitize = function(Y) {
  colnames(Y) = 1:NCOL(Y)
  return(flowCore::flowFrame(Y))
}
  

# Read data.
y = {
  prefix = 'results/phi25-pi0.2/img/txt/y'
  suffix = '_mean.csv'
  lapply(as.list(paste0(prefix, 1:2, suffix)), read.csv, header=F)
}

# Number of cells in each sample
N = sapply(y, nrow)

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

# Print clusterings for all samples
write.table(cbind(fsclus, sample_idx),
            file=paste0("results/tsne/flowsom.csv"),
            quote=F, row.names=F, col.names=F)
