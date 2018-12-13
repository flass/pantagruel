#!/usr/bin/Rscript --vanilla
library(pvclust)
#~ library(parallel)
library(ape)

cargs = commandArgs(trailingOnly=T)


orthomatrad = cargs[1]
nboot = as.numeric(cargs[2])

load(sprintf('%s_gene_abspres.mat.RData', orthomatrad))

jacc.dist = dist(t(genocount), method='binary')
jacc.dist.ward.clust = hclust(jacc.dist, method='ward.D2')
save(jacc.dist, jacc.dist.ward.clust, file=sprintf('%s_gene_abspres.clustering.RData', orthomatrad))
pdf(file=sprintf('%s_gene_abspres.clustering.pdf', orthomatrad), width=16, height=9)
plot(jacc.dist.ward.clust)
dev.off()
write.tree(as.phylo(jacc.dist.ward.clust), file=sprintf('%s_gene_abspres.clustering.nwk', orthomatrad))

jacc.dist.ward.pvclust = pvclust(genocount, method.dist='binary', method.hclust='ward.D2', parallel=T, nboot=nboot)
save(jacc.dist.ward.pvclust, file=sprintf('%s_gene_abspres.pvclustering.RData', orthomatrad))
pdf(file=sprintf('%s_gene_abspres.pvclustering.pdf', orthomatrad), width=16, height=9)
plot(jacc.dist.ward.pvclust)
dev.off()
jacc.dist.ward.pvclust.tree = as.phylo(jacc.dist.ward.pvclust$hclust)
for (support in c('bp', 'au')){
	jacc.dist.ward.pvclust.tree$node.label = as.character(format(jacc.dist.ward.pvclust$edges[,support], digits=2))
	write.tree(jacc.dist.ward.pvclust.tree, file=sprintf('%s_gene_abspres.pvclustering.%s.nwk', orthomatrad, support))
}
jacc.dist.ward.pvclust.tree$node.label = paste(format(jacc.dist.ward.pvclust$edges[,'au'], digits=2),
                                               format(jacc.dist.ward.pvclust$edges[,'bp'], digits=2), sep='/')
write.tree(jacc.dist.ward.pvclust.tree, file=sprintf('%s_gene_abspres.pvclustering.au-bp.nwk', orthomatrad, support))
write.nexus(jacc.dist.ward.pvclust.tree, file=sprintf('%s_gene_abspres.pvclustering.au-bp.nex', orthomatrad, support))
