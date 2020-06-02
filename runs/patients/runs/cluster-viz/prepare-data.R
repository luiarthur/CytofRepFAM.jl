source('preimpute.R')
source('preprocess.R')

# Set random number generator seed for reproducibility.
set.seed(0)

# Data directory.
data_dir = "../../../data/patients/transformed-data/"

# File names>
fnames = c("001_d31_clean.csv",
           "007_d35_clean.csv",
           "010_d35_clean.csv")

# Read data.
y.orig = lapply(as.list(paste0(data_dir, fnames)), read.csv)

# Preprocess data by removing rows with cells < -7 or > 4.
y = lapply(y.orig, preprocess, lower=-7, upper=4)

for (i in 1:length(y)) {
  yi = round(as.matrix(y[[i]]), 5)
  # Preimpute missing values
  yi = preimpute(yi)
  # Write data.
  write.csv(yi, file=paste0("img/y", i, "-prepped.csv"), quote=F, row.names=F)
}

for (i in 1:length(y)) {
  yi = round(as.matrix(y[[i]]), 5)
  write.csv(yi, file=paste0("img/y", i, "-prepped-with-missing.csv"),
            quote=F, row.names=F)
}
