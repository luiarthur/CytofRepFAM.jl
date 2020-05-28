data_dir = "../../data/patients/transformed-data/"

fnames = c("001_d31_clean.csv",
           "007_d35_clean.csv",
           "010_d35_clean.csv")

y = lapply(as.list(paste0(data_dir, fnames)), read.csv)

for (thresh in -4:-7) for (i in 1:3){
  yi = as.matrix(y[[i]])
  prop_bad = apply(yi, 1, function(x) any(x < thresh, na.rm=TRUE))
  cat(paste0('prop. cells where at least 1 marker < ', thresh), paste0('sample: ', i),
      mean(prop_bad), '\n', sep=' | ')
}

cat('\n')

for (thresh in 3:5) for (i in 1:3){
  yi = as.matrix(y[[i]])
  prop_bad = apply(yi, 1, function(x) any(x > thresh, na.rm=TRUE))
  cat(paste0('prop. cells where at least 1 marker > ', thresh), paste0('sample: ', i),
      mean(prop_bad), '\n', sep=' | ')
}

is.bad.row = function(x, lower=-7, upper=4) any(x < -7, na.rm=T) || any(x > 4, na.rm=T)
lower = -7
upper = 4
for (i in 1:3){
  yi = as.matrix(y[[i]])
  prop_bad = apply(yi, 1, function(x) is.bad.row(x, lower, upper))
  cat('prop. cells where at least 1 marker outside of (-7, 4)', paste0('sample: ', i),
      mean(prop_bad), '\n', sep=' | ')
}


