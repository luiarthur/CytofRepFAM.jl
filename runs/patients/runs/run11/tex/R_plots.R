R_plots = function(countspath, savepath=NULL) {
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

make_R_plots = function(path) {
  countspath = paste0(path, '/img/txt/Rcounts.txt')
  savepath = paste0(path, '/img/Rcounts.pdf')
  R_plots(countspath, savepath)
}

# for (phi in c(0, 1, 10, 100)) {
for (phi in c(0, 1, 100, 10000)) {
  path = paste0('results/phi', phi)
  make_R_plots(path)
}

