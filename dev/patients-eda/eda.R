# NOTE: Skip patient 005. The cutoff files are in Excel 2007. I cannot
# open them.

source('read.tcsv.R')
library(stringr)

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
                    rawdatadir=paste0(datadir, 'raw/'),
                    cutoffdir=paste0(datadir, 'cutoff/'),
                    transdatadir=paste0(datadir, 'transformed-data/'),
                    verbose=0) {
  f.contents = read.csv(f)
  colnames(f.contents) = sanitize.markers(colnames(f.contents))
  data.size = dim(f.contents)
  prop.miss = apply(f.contents, 2, function(x) mean(x==0))

  # Create names of output if needed
  outfile = gsub(x=f, pattern=rawdatadir, replacement='')
  outfile = gsub(x=outfile, pattern='.csv', replacement='.pdf')
  outfile = paste0(imgdir, outfile)
  cat('Processing: ', outfile, '\n')

  # Cutoff filepath
  cutoff.path = gsub(x=f, pattern=rawdatadir, replacement=cutoffdir)
  cutoff = read.tcsv(cutoff.path)
  colnames(cutoff) = sanitize.markers(colnames(cutoff))

  # Get marker names
  data.markernames = colnames(f.contents)
  cutoff.markernames = colnames(cutoff)

  if (verbose > 0) {
    print(sort(data.markernames))
    print(sort(cutoff.markernames))
  }
  stopifnot(all(sort(data.markernames) == sort(cutoff.markernames)))

  markers.sorted = sort(data.markernames)

  # Plot data density for each marker
  trans.data = matrix(NA, data.size[1], data.size[2])
  colnames(trans.data) = markers.sorted
  pdf(outfile)
  par(mfrow=c(4, 2), mar=c(5, 4, 1, 2))

  # Column index for transformed data matrix
  j = 0
  for (marker in markers.sorted) {
    j = j + 1
    data.index = which(data.markernames == marker)
    cutoff.index = which(cutoff.markernames == marker)
    numer = as.numeric(f.contents[, data.index])
    denom = as.numeric(cutoff[cutoff.index])
    marker.trans =  log(numer) - log(denom)
    marker.trans = ifelse(is.infinite(marker.trans), NA, marker.trans)
    trans.data[, j] = marker.trans
    hist(marker.trans, breaks=30,
         prob=TRUE, xlab=marker, main='',
         xlim=c(-1, 1) * 6)
    legend('topleft', legend=c(paste0('missing: ',
                                      round(prop.miss[data.index] * 100, 1), '%'),
                               paste('cells:', data.size[1]),
                               paste('sd: ', 
                                     round(sd(marker.trans, na.rm=TRUE), 2))),
           text.font=2, text.col='red', cex=.9, bty='n')
  }
  par(mfrow=c(1, 1), mar=rcommon::mar.default(), oma=rcommon::oma.default())
  dev.off()

  # Write transformed data
  dir.create(transdatadir, showWarn=FALSE)
  fout = gsub(x=f, pattern=rawdatadir, transdatadir)
  write.csv(trans.data, file=fout, quote=FALSE, row.names=FALSE)
  
  c(nrows=data.size[1], ncols=data.size[2], prop.miss=prop.miss)
}

# Data directory
data.dir = 'data/'
raw.data.dir = paste0(data.dir, 'raw/')

# List files in current directory
files = paste0(raw.data.dir, list.files(raw.data.dir))

# Filter out the expression-level files 
data.files = filter.data.files(files)

# FIXME: Removing patient 005 because I can't open the cutoff files.
data.files = data.files[!grepl('005', data.files)]

# Get the dimensions of each data set
data.sizes = sapply(data.files, get.info)

# Save data info
write.csv(data.sizes, file='img/data-info.csv', quote=FALSE, row.names=TRUE)


# Plot histogram of the number of observations (Ni)
pdf('img/hist_nobs.pdf')
hist(data.sizes[1,], breaks=20, xlab='data size', main='Histogram of sample size')
dev.off()

# Number of samples for each patient
days = as.numeric(str_match(data.files, '(?<=_d)\\d+(?=_)'))
pdf('img/days.pdf')
hist(days, breaks=10)
dev.off()

# Number of patients
patients = as.numeric(str_match(data.files, '(?<=/)\\d+(?=_)'))

# TODO: Order patients by days
cells.days.patients = cbind(num_cells=data.sizes[1, ],
                            day=days, patient=patients)
write.csv(cells.days.patients[order(days), ],
          file='img/samples-ordered-by-day.csv',
          quote=FALSE, row.names=TRUE)
