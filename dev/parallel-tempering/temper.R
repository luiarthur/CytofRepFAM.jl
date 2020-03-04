x = seq(-5, 5, len=1e4)
ld1 = dnorm(x, -2, .1, log=TRUE)
ld2 = dnorm(x, 2, .1, log=TRUE)
logsumexp = function(logx, logy) {
  lxy = cbind(logx, logy)
  mx = apply(lxy, 1, max)
  # log(sum(exp(logx - logy))) + logy
  log(rowSums(exp(lxy - mx))) + mx
}

lf1 = function(t) {
  # (.7 * d1 + .3 * d2) ^ (1/t)
  logsumexp(log(.7) + ld1, log(.3) + ld2) / t
}

lf2 = function(t) {
  # .7 * d1 ^ (1/t) + .3 * d2 ^ (1/t)
  logsumexp(log(.7) + ld1/t, log(.3) + ld2/t)
}

f1 = function(t) exp(lf1(t))
f2 = function(t) exp(lf2(t))

ntemp = 20
maxtemp = 10000
degree = 4
ts = maxtemp ^ (((1:ntemp) ^ degree - 1) / (ntemp ^ degree - 1))
tempstep = maxtemp ^ (1/(ntemp - 1))

# Plots
plot(x, f1(1), type='l', lwd=3)
for (ti in ts[-1]) {
  # lines(x, f1(ti))
  lwd=ifelse(ti == maxtemp, 3, 1)
  lines(x, f1(ti), col='grey', lwd=lwd)
  lines(x, f2(ti), col='blue', lwd=lwd)
}

# ntemp=100; maxtemp=100
# for (d in 1:4) {
#   ts = maxtemp ^ (((1:ntemp) ^ d - 1) / (ntemp ^ d - 1))
#   if (d == 1) plot(ts, type='l', lwd=3) else lines(ts)
# }
