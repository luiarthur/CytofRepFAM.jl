library(mclust)

# Directory with patients data  
data_dir = '../../data/patients/transformed-data/'  

set.seed(1)

# Filenames of files around day 30  
fnames_day30 = c('001_d31_clean.csv', 
                 '007_d35_clean.csv', 
                 '010_d35_clean.csv') 

# Samples into a list 
samples = lapply(as.list(fnames_day30), 
                 function(fname) read.csv(paste0(data_dir, fname))) 


# Samples into a matrix
Y = do.call(rbind, samples)
idx = {
  num_samples = length(samples)
  .idx = sapply(samples, nrow)
  cumsum_idx = cumsum(.idx)
  rbind(c(1, cumsum_idx[-num_samples] + 1),
        c(cumsum_idx))
}

# Fill missing values
yneg = Y[Y < -0.5 & !is.na(Y)]
Y[is.na(Y)] <- sample(yneg, size=sum(is.na(Y)), replace=TRUE)
colnames(Y) <- colnames(samples[[1]])
head(Y)

# Mclust
Y.orig = do.call(rbind, samples)
pc.Y = princomp(Y)

for (G in 5:7) {
  for (modelName in c("VVI", "VVV")) {
    # G = 6
    # modelName = "VVI"
    model = mclust::Mclust(Y, G=G, modelNames=modelName)
    # Plot heatmap for each sample
    for (i in 1:length(samples)) {
      png(paste0('img/heatmap/heatmap', '-G', G, '-model', modelName,
                 '-sample', i, '.png'), h=700, w=700)
      idx.i = idx[1, i]:idx[2, i]
      Yi = Y.orig[idx.i, ]
      Zi = Yi[order(model$class[idx.i]), ]
      cytof3::my.image(as.matrix(Zi), col=cytof3::blueToRed(9),
                       na.color='black', zlim=c(-4, 4), addL=TRUE,
                       f=function(Z) cytof3::add.cut(clus=model$class[idx.i]))
      dev.off()

      pdf(paste0('img/pc/pc', '-G', G, '-model', modelName,
                 '-sample', i, '.pdf'))
      plot(pc.Y$scores[idx.i, 1:2], col=model$class[idx.i], cex=.5, pch=20)
      dev.off()
    }
  }
}
