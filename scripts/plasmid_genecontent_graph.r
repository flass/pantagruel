#!/usr/bin/Rscript --vanilla

library(parallel)
library(vegan)
library(igraph)

writeCommunityContent = function(refp.comm, refplasmid, matreplicontent, nfinfammat){
	nfoutrefpcomm = paste(nfinfammat, sprintf("%s-centred_community_size%d.replicons", refplasmid, length(refp.comm)), sep='.')
	write(refp.comm, file=nfoutrefpcomm)
	allfam.refp.group = names(which(apply(matreplicontent[refp.comm,], 2, sum)>0))
	nfoutrefpcommfam = paste(nfinfammat, sprintf("%s-centred_community_size%d.families", refplasmid, length(refp.comm)), sep='.')
	write(sort(allfam.refp.group), file=nfoutrefpcommfam)
}


maxsizerefcomm = 25
minmediansimrefcomm = 0.7
#~ distmethods = c("jaccard", "bray")
distmethods = c("bray")

nfinfammat = file.path(Sys.getenv()['plasmids'], 'plasmid2genefamily.mat_excl-transposases')
refplasmids = c('NC_022651.1', 'NZ_CP008933.1')
#~ nfinfammat = cargs[1]
#~ refplasmids = cargs[2:length(cargs)]

matreplicontent = t(read.table(nfinfammat, h=T, row.names=1))
emptyrepli = which(apply(matreplicontent, 1, sum)==0)
if (length(emptyrepli) > 0){
	matreplicontent = matreplicontent[-emptyrepli, ]
}
for (distmet in distmethods){
	print(distmet)
#~ distreplicontent = vegdist(matreplicontent, method=distmet, binary=T)
	distreplicontent = vegdist(matreplicontent, method=distmet)
	nfoutdistreplicontent = paste(nfinfammat, sprintf("%s_dist.RData", distmet), sep='.')
	save(distreplicontent, file=nfoutdistreplicontent)

	sim.repli = as.matrix(1 - distreplicontent) # diagonal stay zero so no self adjacency
	graph.repli = graph_from_adjacency_matrix(sim.repli, weighted=T, mode="undirected")
	graph.repli.communities = cluster_fast_greedy(graph.repli)
#~ 	graph.repli.communities = cluster_louvain(graph.repli)
	
	nfoutrepligraph = paste(nfinfammat, sprintf("%s_graph.pdf", distmet), sep='.')
	pdf(nfoutrepligraph, width=15, height=15)
	layout(matrix(1:6, 3, 2, byrow=T))
	plot(graph.repli.communities, graph.repli)
	hist(as.dist(sim.repli), main=sprintf("'%s' similarity within complete graph", distmet))
	dev.off()
	
	# try to further split the group containing the reference until it reaches a decent size
	for (refplasmid in refplasmids){
		ncomm = length(graph.repli.communities)
		nfoutrefpcommpdf = paste(nfinfammat, sprintf("%s-centred_communities.pdf", refplasmid), sep='.')
		pdf(nfoutrefpcommpdf, width=15, height=15)
		layout(matrix(1:6, 3, 2, byrow=T))
#~ 		refp.comm = graph.repli.communities[[which(sapply(1:ncomm, function(i){refplasmid %in% graph.repli.communities[[i]]}))]]
		refp.comm = names(V(graph.repli))
		subgraph.p.communities = graph.repli.communities
		n = 0
		
		while (length(refp.comm) > maxsizerefcomm){
#~ 		while ((median(as.dist(sim.repli[refp.comm, refp.comm])) < minmediansimrefcomm) & (ncomm > 1)){
			writeCommunityContent(refp.comm, refplasmid, matreplicontent, nfinfammat)
			n = n + 1
			print(sprintf("round #%d", n), quote=F)
			whereis.refp = which(sapply(1:ncomm, function(i){refplasmid %in% subgraph.p.communities[[i]]}))
			print(sprintf("%s in community #%d", refplasmid, whereis.refp), quote=F)
			refp.comm = subgraph.p.communities[[whereis.refp]]
			subgraph.p = induced_subgraph(graph.repli, refp.comm)
			subgraph.p.communities = cluster_fast_greedy(subgraph.p)
			plot(subgraph.p.communities, subgraph.p, main=sprintf("subgraph focused on %s, round#%d", refplasmid, n))
			print(median(as.dist(sim.repli[refp.comm, refp.comm])))
			hist(as.dist(sim.repli[refp.comm, refp.comm]), main=sprintf("'%s' similarity within subgraph, round#%d", distmet, n))
			ncomm = length(subgraph.p.communities)
			commsizes = sapply(1:ncomm, function(i){length(subgraph.p.communities[[i]])}) ; names(commsizes) = 1:ncomm
			print(commsizes)
		}
		
		dev.off()
		writeCommunityContent(refp.comm, refplasmid, matreplicontent, nfinfammat)
	}
}
