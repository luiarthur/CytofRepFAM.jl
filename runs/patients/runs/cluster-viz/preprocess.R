is.bad.row = function(x, lower, upper) {
  any(x < lower, na.rm=TRUE) || any(x > upper, na.rm=TRUE)
}

is.good.row = function(x, lower, upper) !is.bad.row(x, lower, upper)

# Remove rows that have very extreme values (i.e. any cell > upper or any cell < lower)
preprocess = function(yi, lower, upper) {
  idx.keep = which(apply(yi, 1, is.good.row, lower, upper))
  yi[idx.keep, ]
}
