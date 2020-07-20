R_plots = function(countspath, savepath=NULL, xlim=NULL) {
  R = read.table(countspath)
  I = NROW(R)
  rownames(R) = paste0('S', 1:I, ' counts')

  if (!is.null(savepath)) {
    pdf(savepath)
  }

  plot.ts(t(R), xlab='number of cell subpopulations', cex=2,
          type='b', pch=20, main='', cex.lab=1.5, cex.axis=1.3,
          mar.multi=c(2, 5.1, 0, 2.1))

  if (!is.null(savepath)) {
    dev.off()
  }
}

make_R_plots = function(path, xlim=NULL) {
  countspath = paste0(path, '/img/txt/Rcounts.txt')
  savepath = paste0(path, '/img/Rcounts.pdf')
  R_plots(countspath, savepath, xlim=xlim)
}

# NOTE: edit this.
Z_dist_plots = function(z1_path, z2_path, w1_path, w2_path, savepath) {
  Z1 = read.table(z1_path)
  Z2 = read.table(z2_path)
  W1 = read.table(w1_path)
  W2 = read.table(w2_path)

  z1_pairwise_dist = c(dist(t(Z1[, W1>0]), method='manhattan'))
  z2_pairwise_dist = c(dist(t(Z2[, W2>0]), method='manhattan'))
  npairs = length(z1_pairwise_dist)
  J = dim(Z1)[1]

  pdf(savepath)
  par(mfrow=c(2, 1))
  plot(table(z1_pairwise_dist)/npairs, xlim=c(0, J), 
       ylab='proportion of column pairs', xlab='Manhattan distance',
       main='Z estimate for Sample 1')
  # abline(v=mean(z1_pairwise_dist), lwd=3, col='red')
  plot(table(z2_pairwise_dist)/npairs, xlim=c(0, J), 
       ylab='proportion of column pairs', xlab='Manhattan distance',
       main='Z estimate for Sample 2')
  # abline(v=mean(z2_pairwise_dist), lwd=3, col='red')
  par(mfrow=c(1, 1))
  dev.off()
}

# NOTE: edit this.
make_Z_dist_plots = function(path) {
  z1_path = paste0(path, '/img/txt/Z1_best.txt')
  z2_path = paste0(path, '/img/txt/Z2_best.txt')
  w1_path = paste0(path, '/img/txt/W1_best.txt')
  w2_path = paste0(path, '/img/txt/W2_best.txt')
  savepath = paste0(path, '/img/Z_dist.pdf')
  Z_dist_plots(z1_path, z2_path, w1_path, w2_path, savepath)
}


### MAIN ###
for (phi in c(0, 1, 100, 10000)) {
  path = paste0('results/phi', phi)
  make_R_plots(path, xlim=c(10, 25))  # NOTE: mind xlim!
  make_Z_dist_plots(path)
}

