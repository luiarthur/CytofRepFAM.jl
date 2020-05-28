beta = matrix(c(-14.01933, -14.78794, -15.41919,
                 -9.51325,  -9.16730,  -9.68004,
                 -1.13317,  -1.06930,  -1.12893), nrow=3)

sigmoid = function(x) 1 / (1 + exp(-x))
pm = function(y, beta) {
  sigmoid(beta[1] + beta[2] * y + beta[3] * y^2)
}
ys = seq(-10, 6, length=100)

pdf('img/missmech.pdf')
plot(ys, pm(ys, beta[1, ]), type='n', xlab='y', ylab='prob miss')
for (i in 1:3) {
  lines(ys, pm(ys, beta[i, ]), type='l', col=i + 1, lwd=3)
}
legend('topright', legend=paste0('sample: ', 1:3), lwd=3, col=2:4)
dev.off()
