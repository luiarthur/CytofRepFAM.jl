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
  K = dim(Z1)[2]

  pdf(savepath)
  par(mfrow=c(2, 1), mar=c(4.1, 5.1, .2, 2.1))
  plot(table(z1_pairwise_dist)/npairs1, xlim=c(0, J), 
       xlab='',
       main='', cex.lab=1.5, cex.axis=1.3, lwd=3,
       ylab='Proportion in Sample 1', xaxt='n')
  # abline(v=mean(z1_pairwise_dist), lwd=3, col='red')
  axis(1, 1:J, 1:J, cex.axis=1.3, las=2)

  plot(table(z2_pairwise_dist)/npairs2, xlim=c(0, J), 
       xlab='Pairwise column distance',
       main='', cex.lab=1.5, cex.axis=1.3, lwd=3,
       ylab='Proportion in Sample 2', xaxt='n')
  # abline(v=mean(z2_pairwise_dist), lwd=3, col='red')
  axis(1, 1:J, 1:J, cex.axis=1.3, las=2)

  par(mfrow=c(1, 1), mar=c(5.1, 4.1, 4.1, 2.1))
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
# NOTE: Settings:
phis = 25
p_is = c(0.0, 0.1, 0.2, 0.3)
results_dir = '../results/phi'
xlim = c(10, 25)
ylim = c(0, 0.8)

# Run
for (phi in phis) for (p_i in p_is) {
  path = paste0(results_dir, phi, '-pi', format(p_i, nsmall=1))
  make_R_plots(path, xlim=xlim, ylim=ylim)
  make_Z_dist_plots(path)
}
