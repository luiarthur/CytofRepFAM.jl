rtnorm = function(n) qnorm(runif(n, .5, 1))

rsn = function(n, loc, scale, skew) {
  z = if (skew == 0) 0 else rtnorm(n)
  loc + scale * skew * z + scale * sqrt(1 - skew^2) * rnorm(n)
}

x = rsn(1e5, 1, .5, -.97) 
plot(density(x))
abline(v=0, col='red', lwd=3)
cat('Pr(x < 0):', mean(x < 0), '\n')
