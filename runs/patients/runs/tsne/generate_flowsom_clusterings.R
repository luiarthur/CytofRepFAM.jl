data_dir = "../../../data/patients/transformed-data/"

fnames = c("001_d31_clean.csv",
           "007_d35_clean.csv",
           "010_d35_clean.csv")

y = lapply(as.list(paste0(data_dir, fnames)), read.csv)

i = 1
yi = y[[1]]

