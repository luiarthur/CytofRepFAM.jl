library(stringr)
library(mclust)
library(FlowSOM)

get_mclust_clusterings = function(Y, max_clusters, modelNames='VII', seed=0) {
  set.seed(seed)
  result = Mclust(Y, G=1:max_clusters, modelNames=modelNames)
  result$classification
}

flowsanitize = function(Y) {
  colnames(Y) = 1:NCOL(Y)
  return(flowCore::flowFrame(Y))
}
  
get_flowsom_clusterings = function(Y, max_clusters, seed=0) {
  # Prep data to be in FlowSOM data format
  ff_Y = flowsanitize(as.matrix(Y))

  # Numbe of markers
  num_markers = ncol(Y)

  # Flowsom clusterings
  flowsom.result = FlowSOM(ff_Y,
                           colsToUse=1:num_markers,
                           maxMeta=max_clusters,
                           seed=seed)
  # Get cluster labels
  fsmeta = flowsom.result$meta
  fsclus = fsmeta[flowsom.result$FlowSOM$map$mapping[, 1]]

  return(as.numeric(fsclus))
}

### MAIN ###

# Output directory.
outdir = 'viz/csv'

# Maximum number of clusters to use.
max_clusters = 10

# Loop through data sets.
for (pmiss in c(0.0, 0.6)) for (zind in 1:3) {
  # Simulation name.
  simname = str_interp('pmiss$[.1f]{pmiss}-phi0.0-zind${zind}')

  # Path to data.
  path_to_data = str_interp('${outdir}/tsne-${simname}.csv')
  cat('Processing:', path_to_data, '\n')

  # Read data
  df = read.csv(path_to_data)

  # Get relevant columns for clustering.
  ycolumns = grepl('Y\\d+', colnames(df))
  Y = df[, ycolumns]

  # Do clustering (mclust, flowsom).
  # A variety of seeds were tried, and the best result was used.
  mc.cluster = get_mclust_clusterings(Y, max_clusters, seed=0, modelNames='VII')
  fs.cluster = get_flowsom_clusterings(Y, max_clusters, seed=2)

  # Write results.
  write.table(mc.cluster, str_interp('${outdir}/mclust-${simname}.csv'),
              quote=F, row.names=F, col.names=F)
  write.table(fs.cluster, str_interp('${outdir}/flowsom-${simname}.csv'),
              quote=F, row.names=F, col.names=F)
}
