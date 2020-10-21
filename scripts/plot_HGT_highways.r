#!/usr/bin/env Rscript
library(ggplot2)
library(ggtree)
library(igraph)

ptgdb = '/Users/fl4/pantagruel_databases/choleraclub315plus7'

parsedrecdir=file.path(ptgdb, '07.reconciliations/fullgenetree_GeneRax_recs/nocollapse/noreplace/generax_fullgenetree_recsampling_1/reconciliations/')
#nfspetree = file.path(ptgdb, '05.core_genome/core-genome-based_reference_tree_choleraclub315plus7.full_clade_defs.nwk')
nfspetree = file.path(ptgdb, '05.core_genome/core-genome-based_reference_tree_choleraclub315plus7.only314_clade_defs.nwk')
spetree = read.tree(nfspetree)

presentoutput = '~/Documents/presentations/group_and_collab_meetings/cholera_club_pangenome_2020-10-19'

read.parsed.rec = function(nfparsedrec, event.filter=c('T'), freq.thresh=10, max.freq=1000, scale.freq=TRUE, recdon.filter=NULL){
	parsed.rec = read.table(nfparsedrec, sep='\t', header=F, col.names=c('evtype', 'rec', 'don', 'freq'))
	parsed.rec = parsed.rec[parsed.rec[,'evtype'] %in% event.filter & parsed.rec[,'freq'] > freq.thresh,]
	if (!is.null(recdon.filter)){
		parsed.rec = parsed.rec[as.character(parsed.rec[,'don']) %in% recdon.filter & as.character(parsed.rec[,'rec']) %in% recdon.filter,]
	}
	if (scale.freq){
		parsed.rec[,'freq'] = as.double(parsed.rec[,'freq'] / max.freq)
	}
	return(parsed.rec[, c('rec', 'don', 'freq')])
}

all.parsed.recs.fams = sapply(strsplit(list.files(parsedrecdir), '_'), `[`, 1)
all.parsed.recs = lapply(list.files(parsedrecdir, full.names=T), read.parsed.rec, recdon.filter=c(spetree$tip.label, spetree$node.label))
names(all.parsed.recs) = all.parsed.recs.fams

combined.recs = list()
for (parsed.rec in all.parsed.recs){
	if (nrow(parsed.rec)>0){
		recdon = paste(parsed.rec[,'rec'], parsed.rec[,'don'], sep='_')
		for (i in 1:nrow(parsed.rec)){
			rd = recdon[i]
			rdf = combined.recs[[rd]]
			if (is.null(rdf)){
				combined.recs[[rd]] = parsed.rec[i,'freq']
			}else{
				combined.recs[[rd]] = parsed.rec[i,'freq'] + rdf
			}
		}
	}
}
combined.recs.df = as.data.frame(t(simplify2array(strsplit(names(combined.recs), split='_'))), col.names=c('rec', 'don'))
combined.recs.df[['freq']] = unlist(combined.recs)
top.combined.recs = combined.recs.df[combined.recs.df[['freq']]>=1,] ; colnames(top.combined.recs) = c('rec', 'don', 'freq')
save(top.combined.recs, file=file.path(presentoutput, 'HGTfreq_combinedfams_topevents.RData'))

top.combined.recs = combined.recs.df[combined.recs.df[['freq']]>=3,] ; colnames(top.combined.recs) = c('rec', 'don', 'freq')

#p4 <- ggtree(spetree, layout="inward_circular", xlim=c(150, 0)) +
p <- ggtree(spetree, layout="inward_circular")

#p1 = p + geom_taxalink(data=all.parsed.recs[[15]], 
p1 = p + geom_taxalink(data=top.combined.recs, 
                         mapping=aes(taxa1=don, taxa2=rec, color=freq, size=freq)) + 
          scale_colour_viridis_c()

ggsave(file=file.path(presentoutput, 'HGTfreq_inwardcirculartree.pdf'), device='pdf', plot=p1)

pdf(file=file.path(presentoutput, 'HGTfreq_graphtree.pdf'), width=13, height=8)
spetreegraph = graph_from_edgelist(spetree$edge)
spetreegraph = set_vertex_attr(spetreegraph, 'name', V(spetreegraph), c(spetree$tip.label, spetree$node.label))
spetreegraph = set_edge_attr(spetreegraph, 'length', E(spetreegraph), spetree$edge.length)
spetreegraph = set_edge_attr(spetreegraph, 'weight', E(spetreegraph), 1)
plot(spetreegraph, vertex.size=.1, vertex.label=NA)

hgtgraph = graph_from_edgelist(as.matrix(top.combined.recs[,1:2]))
hgtgraph = set_edge_attr(hgtgraph, 'weight', E(hgtgraph), top.combined.recs[,'freq'])
hgt.freq.cols = c('grey', rainbow(15))
plot(hgtgraph, vertex.size=.1,
	edge.color=hgtgraph[as.integer(round(edge_attr(hgtgraph, 'weight', E(hgtgraph))))])

spetree.hgt.graph = union(spetreegraph, hgtgraph, byname=T)
w1 = edge_attr(spetree.hgt.graph, 'weight_1', E(spetree.hgt.graph))
w2 = edge_attr(spetree.hgt.graph, 'weight_2', E(spetree.hgt.graph))
w = ifelse(is.na(w1), w2, w1)
spetree.hgt.graph = set_edge_attr(spetree.hgt.graph, 'weight', E(spetree.hgt.graph), w)

plot(spetree.hgt.graph, vertex.size=.1, vertex.label=NA,
	 edge.color=hgt.freq.cols[as.integer(round(edge_attr(spetree.hgt.graph, 'weight', E(spetree.hgt.graph))))])

dev.off()

library(phytools)


asscode = read.table('/Users/fl4/pantagruel_databases/choleraclub315plus7/03.database/genome_codes.tab', header=F, comment.char='', sep='\t') 
codes = asscode[,2] ; names(codes) = gsub('#', '.', asscode[,1])

spetree$tip.label = gsub('#', '.', spetree$tip.label)

pangenome.mat = read.table('/Users/fl4/pantagruel_databases/choleraclub315plus7/02.gene_alignments/full_families_genome_counts-noORFans.mat', header=T, comment.char='', row.names=1)

ass = gsub('^X', '', colnames(pangenome.mat))
ass[ass=="RKI.ZBS2.CH129_TACAGC_L002.contigs_spades.1"] = 'RKI-ZBS2-CH129_TACAGC_L002.contigs_spades.1'
colnames(pangenome.mat) = codes[ass]
pangenome.mat[pangenome.mat>1] = 1

tpangenome.mat = t(pangenome.mat)[spetree$tip.label, (apply(pangenome.mat, 1, sum) %in% 16:(322-16))]
pangenome.clust = hclust(dist(t(tpangenome.mat)))
pdf(file=file.path(presentoutput, 'pangenome_tree.pdf'), width=13, height=8)
phylo.heatmap(spetree, tpangenome.mat[,pangenome.clust$order], labels=F)
dev.off()
