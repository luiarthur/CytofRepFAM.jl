# NOTE:
# In this file, I preimpute missing expression levels for each marker and
# sample  by subsampling observed expression levels in the same marker and
# sample. I subsample only observed values between a certain range (e.g.
# between -6 and 0), so that I don't sample values that are too extreme on
# either side. The limitation is I may not get a very smooth distribution.

Y = read.csv('../../data/cb_transformed.csv')

# Get sample 1
y1 = Y[Y$sample_id == 1, 1:32]

# Resample observed values, giving more weight to values smaller values
sample_obs = function(n, yij_obs, thresh=c(-6, 0)) {
  lower = thresh[1]
  upper = thresh[2]
  population = yij_obs[yij_obs > lower & yij_obs < upper]
  sample(population, n, replace=TRUE)
}

j = 8
nmiss_1j = sum(is.na(y1[, j]))

y1j_new = y1[, j]
mis_idx = which(is.na(y1[, j]))
y1j_new[mis_idx] = sample_obs(nmiss_1j, y1[-mis_idx, j])

hist(y1j_new, border=NA, col=rgb(1, 0, 0, .5), nclass=50, xlim=range(y1[-mis_idx, j]))
hist(y1[, j], border=NA, col=rgb(0, 0, 1, .5), add=TRUE, nclass=50)
mean(is.na(y1[, j]))
