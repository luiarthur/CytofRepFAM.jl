filter.data.files = function(files) {
  out = c()

  for (f in files) {
    is.data.file = grepl('clean.csv', f, ignore.case=TRUE)
    if (is.data.file) {
      out = c(out, f)
    }
  }

  return(out)
}

# List files in current directory
files = list.files('data')

# Filter out the expression-level files 
data.files = filter.data.files(files)

# Get the dimensions of each data set
data.sizes = sapply(data.files, function(f) dim(read.csv(f)))

# Plot histogram of the number of observations (Ni)
hist(t(data.sizes)[,1], breaks=20)
