library(mclust)

# Directory with patients data  
data_dir = '../../data/patients/transformed-data/'  

set.seed(1)

# Filenames of files around day 30  
fnames_day30 = c('001_d31_clean.csv', 
                 '007_d35_clean.csv', 
                 '010_d35_clean.csv') 

# Samples into a list 
samples = lapply(as.list(fnames_day30), 
                 function(fname) read.csv(paste0(data_dir, fname))) 


# Samples into a matrix
Y = do.call(rbind, samples)

# Fill missing values
yneg = Y[Y < -0.5 & !is.na(Y)]
Y[is.na(Y)] <- sample(yneg, size=sum(is.na(Y)), replace=TRUE)
colnames(Y) <- colnames(samples[[1]])
head(Y)

# Kmeans
# Ysub = Y[sample(1:nrow(Y), 2000), ]
km = kmeans(Y, centers=30, iter=1000)
clus = km$clus

pdf("img/clus.pdf")
cytof3::my.image(as.matrix(Y[order(clus), ]),
                 col=cytof3::blueToRed(9), zlim=c(-4, 4))
cytof3::add.cut(clus=clus)

cytof3::my.image(km$centers > 0)
abline(v=0:32 + .5, h=0:nrow(km$centers) + .5)

ucenters = unique(km$centers > 0, MARGIN=1)
cytof3::my.image(ucenters)
abline(v=0:32 + .5, h=0:nrow(ucenters) + .5)
dev.off()

# mclust
# BIC = mclust::mclustBIC(Ysub)
# pdf('img/bic.pdf')
# plot(BIC)
# dev.off()
# 
# BIC_big = mclust::mclustBIC(Ysub, 1:30)
# pdf('img/bic-big.pdf')
# plot(BIC_big)
# dev.off()

# Requries n <= d
# BIC_all = mclust::mclustBIC(Y, c(1, seq(5, 50, by=5)), modelNames="VII")

BIC_all = mclust::mclustBIC(Y, c(1, seq(5, 50, by=5)), modelNames="EEE")
pdf('img/bic-all.pdf')
plot(BIC_all)
dev.off()

BIC_vvv = mclust::mclustBIC(Y, seq(1, 20, by=2),
                            modelNames=c("VVV", "VVI"))
pdf('img/bic-vvv.pdf')
plot(BIC_vvv)
plot(diff(BIC_vvv), legendArgs=list(x="topright"), ylab="BIC differences")
abline(h=0, lty=2, col='grey')
dev.off()

mclust_vvv = Mclust(Y, G=5, modelNames=c("VVV"))
pdf('img/mclust-vvv.pdf')
mclus = mclust_vvv$class
cytof3::my.image(as.matrix(Y[order(mclus), ]),
                 col=cytof3::blueToRed(9), zlim=c(-4, 4))
cytof3::add.cut(clus=mclus)
dev.off()
