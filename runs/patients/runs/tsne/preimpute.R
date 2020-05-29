preimpute = function(y, seed=0) {
  # set seed for reproducibility
  set.seed(0)

  # Get negative values of y
  yneg = yi[which(yi < 0)]

  # Get indices of all missing
  idx_na = which(is.na(yi))

  # Get number of missing data
  num_na = length(idx_na)

  # Replace missing data with samples of the negatives
  y[idx_na] = sample(yneg, num_na, replace=TRUE)

  return(y)
}

