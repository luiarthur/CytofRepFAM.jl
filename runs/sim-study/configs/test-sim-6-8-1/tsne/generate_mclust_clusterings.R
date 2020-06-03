library(mclust)

get_mclust_clusterings = function(Y, num_clusters, modelNames='VII') {
  result = Mclust(Y, G=num_clusters, modelNames=modelNames)
  result$classification
}



# path_to_data = 'viz/csv/tsne-pmiss0.0-phi0.0-zind1.csv'
# df = read.csv(path_to_data)
# ycolumns = grepl('Y\\d+', colnames(df))
# Y = df[, ycolumns]
# 
# result = get_mclust_clusterings(Y, 7)
# result
