# Directory with patients data
data_dir = '../../data/transformed-patients-data/'

# Filenames of files around day 30
fnames_day30 = c('001_d31_clean.csv',
                 '007_d35_clean.csv',
                 '010_d35_clean.csv')

# Samples into a list
samples = lapply(as.list(fnames_day30),
                 function(fname) read.csv(paste0(data_dir, fname)))

# SD of the markers for each sample
sd_samples = rbind(apply(samples[[1]], 2, sd, na.rm=T),
                   apply(samples[[2]], 2, sd, na.rm=T),
                   apply(samples[[3]], 2, sd, na.rm=T))

# Sample names
sample_names = gsub(x=fnames_day30, pattern='_clean\\.csv', replace='')

# Give the samples names
rownames(sd_samples) = sample_names

# Print and save result
write.csv(t(sd_samples),
          file='out/sd_samples.csv',
          quote=FALSE, row.names=TRUE, col.names=TRUE)
