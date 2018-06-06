#!/usr/bin/Rscript --vanilla
library(gplots)
library(RColorBrewer)
library(igraph)

plotHMevents = function(mat, matname, cap=99){
#~ 	qmat = quantile(mat, (1:100)/100, na.rm=T)
#~ 	widthbin = qmat[cap]/8
#~ 	brk = c(seq(from=0, to=qmat[cap], by=widthbin), max(mat, na.rm=T))
	brk = quantile(mat, (0:9)/9, na.rm=T)
	print(brk)
	heatmap.2(mat, Rowv=F, Colv=F, dendrogram='none', scale='none', trace='none', breaks=brk, col=brewer.pal(9, 'YlOrRd'), na.col='white', main=matname)
#~ 		heatmap.2(mat, Rowv=F, Colv=F, dendrogram='none', scale='none', trace='none', col=brewer.pal(9, 'YlOrRd'), na.col='white', main=matname)
	for (i in 1:(ncol(mat)-1)){
	 for (j in 2:ncol(mat)){
	  mat[j,i] = mat[i,j]
	}}
	layout(matrix(1:2, 2,1))
	barplot(apply(mat, 1, sum)/nsample, las=2, main=sprintf("%s\nsummed matched event freq.", matname))
	n = 5
	barplot(sapply(1:nrow(mat), function(i){
		js = sapply(seq(i-n, i+n), function(k){ ifelse(k>0, ifelse(k<=ncol(mat), k, k-nrow(mat)), k+ncol(mat)) })
		sum(mat[i,-js]) 
	})/nsample, las=2, main=sprintf("%s\nsummed matched event freq. excluding %d next genes", matname, n))
}


#~ nfmatchevents = '/enterobac/PanteroDB_v0.3/07.compare_scenarios/bs_stem70_withinmedian35/replaceCCinGasinS-collapsePOPinSnotinG/ale_collapsed_undat/gene_event_matches/matching_events_10253genes_06062018.tab'
#~ nfmatchevents = '/enterobac/PanteroDB_v0.3/07.compare_scenarios/bs_stem70_withinmedian35/replaceCCinGasinS-collapsePOPinSnotinG/ale_collapsed_undat/gene_event_matches/plasmid2genefamily.mat_excl-transposases.NC_022651.1-centred_community_size13/matching_events.tab'
#~ refrepli = 'NC_022651.1'
#~ nfrefrepliordlineages = '/enterobac/PanteroDB_v0.3/07.compare_scenarios/bs_stem70_withinmedian35/replaceCCinGasinS-collapsePOPinSnotinG/ale_collapsed_undat/gene_event_matches/plasmid2genefamily.mat_excl-transposases.NZ_CP008933.1-centred_community_size24/NC_022651.1_ordered_gene_lineage_products.tab'
cargs = commandArgs(trailingOnly=T)
nfmatchevents = cargs[1]
if (length(cargs)>1){
	refrepli = cargs[2]
	nfrefrepliordlineages = cargs[3]
}else{
	nfrefrepliordlineages = NULL
	refrepli = NULL
}

matchev = unique(read.table(nfmatchevents)) ; colnames(matchev) = c('rlocdsid1', 'rlocdsid2', 'freq')
qmatchev = quantile(matchev$freq, p=c(0:100/100, .999, .9999, .99999, 1.0))
# top 1.0% have co-evol score in [11.0; 23.3] (2M pairwise associations)
# top 0.1% have co-evol score in [13.5; 23.3] (200k pairwise associations)
# top 0.01% have co-evol score in [15.7; 23.3] (20k pairwise associations)
# top 0.001% have co-evol score in [17.7; 23.3] (2k pairwise associations)
topmatchev = matchev[matchev$freq>=qmatchev['99.999%'],]
u1 = unique(topmatchev$rlocdsid1)
u2 = unique(topmatchev$rlocdsid2)
u = intersect(u1, u2)
# ~100 genes tested both ways, with strong association
graph_topmatchev = graph_from_data_frame(topmatchev)

nfgraphout = paste(nfmatchevents, 'top_associations.RData', sep='.')
save(topmatchev, graph_topmatchev, file=nfgraphout)


if (!is.null(nfrefrepliordlineages)){
	refrepli = cargs[2]
	nfrefrepliordlineages = cargs[3]
	refrepliordlineages = read.table(nfrefrepliordlineages, h=T, sep='\t')

	matmatchev = data.matrix(sapply(refrepliordlineages$rlocds_id, function(rlocdsid1){
		sapply(refrepliordlineages$rlocds_id, function(rlocdsid2){
			i = which((matchev$rlocdsid1==rlocdsid1 & matchev$rlocdsid2==rlocdsid2) | (matchev$rlocdsid1==rlocdsid2 & matchev$rlocdsid2==rlocdsid1))
			if (length(i)==1){ return(matchev$freq[i])
			}else{ if (length(i)==0){ return(NA) 
			}else{ print(matchev[i,]) ; stop(sprintf('multiple co-evolution score values for lineage pair %d, %d', rlocdsid1, rlocdsid2)) }}
		})
	}))



	nfpdf = paste(nfmatchevents, sprintf('co-evol_scores_projection_%s.pdf', refrepli) sep='.')
	pdf(nfpdf, height=20, width=20)
	plotHMevents(matmatchev, sprintf('co-evolution scores projected on %s map', refrepli))
	dev.off()
	print(nfpdf)
}
