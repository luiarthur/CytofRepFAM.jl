R_plots = function(countspath, savepath=NULL, ...) {
  R = read.table(countspath)
  I = NROW(R)
  K = NCOL(R)

  if (!is.null(savepath)) {
    pdf(savepath)
  }

  prop = apply(R, 1, function(r) r/sum(r))
  par(mfrow=c(I, 1), mar=c(4.1, 5.1, .2, 2.1))
  for (i in 1:I) {
    xlab = ifelse(i == I, 'number of selected subpopulations', '')
    plot(prop[,i], xlab=xlab, cex=2,
         type='b', pch=20, main='', cex.lab=1.5, cex.axis=1.3,
         ylab=paste0('Proportion in Sample ', i), xaxt='n', ...)
    axis(1, 1:K, 1:K, cex.axis=1.3, las=2)
  }
  par(mfrow=c(1, 1), mar=c(5.1, 4.1, 4.1, 2.1))

  if (!is.null(savepath)) {
    dev.off()
  }
}

make_R_plots = function(path, ...) {
  countspath = paste0(path, '/img/txt/Rcounts.txt')
  savepath = paste0(path, '/img/Rcounts.pdf')
  R_plots(countspath, savepath, ...)
}

# NOTE: edit this.
Z_dist_plots = function(z1_path, z2_path, w1_path, w2_path, savepath) {
  Z1 = read.table(z1_path)
  Z2 = read.table(z2_path)
  W1 = read.table(w1_path)
  W2 = read.table(w2_path)

  z1_pairwise_dist = c(dist(t(Z1[, W1>0]), method='manhattan'))
  z2_pairwise_dist = c(dist(t(Z2[, W2>0]), method='manhattan'))
  npairs1 = length(z1_pairwise_dist)
  npairs2 = length(z2_pairwise_dist)
  J = dim(Z1)[1]

  pdf(savepath)
  par(mfrow=c(2, 1))

  plot(table(z1_pairwise_dist)/npairs1, xlim=c(0, J), 
       ylab='proportion of column pairs', xlab='Manhattan distance',
       main='Z estimate for Sample 1', lwd=3, xaxt='n')
       # type='b', pch=20, cex.lab=1.5, cex.axis=1.3, cex=2)
  # abline(v=mean(z1_pairwise_dist), lwd=3, col='red')
  axis(1, 1:J, 1:J, las=2)

  plot(table(z2_pairwise_dist)/npairs2, xlim=c(0, J), 
       ylab='proportion of column pairs', xlab='Manhattan distance',
       main='Z estimate for Sample 2', lwd=3, xaxt='n')
       # type='b', pch=20, cex.lab=1.5, cex.axis=1.3, cex=2)
  # abline(v=mean(z2_pairwise_dist), lwd=3, col='red')
  axis(1, 1:J, 1:J, las=2)

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
for (phi in c(0, 1, 10, 25, 100)) {  # NOTE: mind phi
  for (pmiss in c(0.0, 0.2)) {
    for (z in 1:3) {
      path = paste0('../results/pmiss', format(pmiss, nsmall=1), '-phi', phi, '-zind', z)
      make_R_plots(path, xlim=c(3, 15), ylim=c(0, 1))  # NOTE: mind xlim!
      make_Z_dist_plots(path)
    }
  }
}
