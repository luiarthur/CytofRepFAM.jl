# NOTE: Skip patient 005. The cutoff files are in Excel 2007. I cannot
# open them.

source('read.tcsv.R')

# Sanitize marker name
sanitize.marker = function(marker) {
  out = toupper(marker)
  out = gsub(x=out, pattern='\\.|-|\\s', '')
  out = gsub(x=out, pattern='KIR', '')
  out = gsub(x=out, pattern='X', '')
  out = gsub(x=out, pattern='\\w+SYK', 'SYK')
  out = gsub(x=out, pattern='2D4', '2DS4')
  out
}

# Sanitize marker names
sanitize.markers = function(markers) {
  as.character(sapply(markers, sanitize.marker))
}


# Keep data files that end with variants of 'clean.csv'.
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

# Function to get dimensions, percent of missing for each file.
# Make plots of observed data?
get.info = function(f, imgdir='img/', datadir=data.dir,
                    cutoffdir=paste0(datadir, 'cutoff/')) {
  f.contents = read.csv(f)
  colnames(f.contents) = sanitize.markers(colnames(f.contents))
  data.size = dim(f.contents)
  prop.miss = apply(f.contents, 2, function(x) mean(x==0))

  # Create names of output if needed
  outfile = gsub(x=f, pattern=datadir, replacement='')
  outfile = gsub(x=outfile, pattern='.csv', replacement='.pdf')
  outfile = paste0(imgdir, outfile)
  cat('Processing: ', outfile, '\n')

  # Cutoff filepath
  cutoff.path = gsub(x=f, pattern=datadir, replacement=cutoffdir)
  cutoff = read.tcsv(cutoff.path)
  colnames(cutoff) = sanitize.markers(colnames(cutoff))

  # Get marker names
  data.markernames = colnames(f.contents)
  cutoff.markernames = colnames(cutoff)

  print(sort(data.markernames))
  print(sort(cutoff.markernames))
  stopifnot(all(sort(data.markernames) == sort(cutoff.markernames)))

  # Plot data density for each marker
  pdf(outfile)
  par(mfrow=c(4, 2), mar=c(5, 4, 1, 2))
  for (j in 1:data.size[2]) {
    markername = data.markernames[j]
    cutoff.index = which(cutoff.markernames == markername)
    marker.trans = log(as.numeric(f.contents[, j])) / as.numeric(cutoff[cutoff.index])
    marker.trans = ifelse(is.infinite(marker.trans), NA, marker.trans)
    hist(marker.trans, breaks=20,
         prob=TRUE, xlab=markername, main='')
    legend('topleft', legend=c(paste0('missing: ', round(prop.miss[j] * 100, 1), '%'),
                               paste('cells:', data.size[1])),
           text.font=2, text.col='red', cex=.9, bty='n')
  }
  par(mfrow=c(1, 1), mar=rcommon::mar.default(), oma=rcommon::oma.default())
  dev.off()

  c(nrows=data.size[1], ncols=data.size[2], prop.miss=prop.miss)
}

# Data directory
data.dir = 'data/'

# List files in current directory
files = paste0(data.dir, list.files(data.dir))

# Filter out the expression-level files 
data.files = filter.data.files(files)

# FIXME: Removing patient 005 because I can't open the cutoff files.
data.files = data.files[!grepl('005', data.files)]

# Get the dimensions of each data set
data.sizes = sapply(data.files, get.info)

# Plot histogram of the number of observations (Ni)
pdf('img/hist_nobs.pdf')
hist(data.sizes[1,], breaks=20, xlab='data size', main='Histogram of sample size')
dev.off()
