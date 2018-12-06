#!/usr/bin/Rscript --vanilla
library(pvclust)
library(parallel)

cargs = commandArgs(trailingOnly=T)


orthomatrad = cargs[1]
nboot = cargs[2]

load(sprintf('%s_gene_abspres.mat.RData', orthomatrad))

jacc.dist = dist(t(genocount), method='binary')
jacc.dist.ward.clust = hclust(jacc.dist, method='ward.D2')
jacc.dist.complete.clust = hclust(jacc.dist, method='average')
save(jacc.dist, jacc.dist.ward.clust, jacc.dist.complete.clust, file=sprintf('%s_gene_abspres.clustering.RData', orthomatrad))
pdf(file=sprintf('%s_gene_abspres.clustering.pdf', orthomatrad), width=16, height=9)
plot(jacc.dist.ward.clust)
plot(jacc.dist.complete.clust)
dev.off()

jacc.dist.ward.pvclust = pvclust(genocount, method.dist='binary', method.hclust='ward.D2', parallel=T, nboot=nboot)
save(jacc.dist.ward.pvclust, file=sprintf('%s_gene_abspres.pvclustering.RData', orthomatrad))
pdf(file=sprintf('%s_gene_abspres.clustering.pdf', orthomatrad), width=16, height=9)
plot(jacc.dist.ward.pvclust)
dev.off()
