#/usr/bin/R
library(ape)

#~ cargs = commandArgs(trailing=T)
#~ nflineevt = cargs[1]
#~ nfspetree = cargs[2]
nflineevt = '/home/flassall/PanteroDB_v0.3/07.compare_scenarios/bs_stem70_withinmedian35/replaceCCinGasinS-collapsePOPinSnotinG/ale_collapsed_undat/summary_gene_tree_events'
nfspetree = '/home/flassall/PanteroDB_v0.3/04.core_genome/raxml_tree/RAxML_bipartitions.pseudo-core-878-unicopy_concat_cds_880entero.166paramBS_MADrooted.full.lsd_internalPopulations.nwk'
lineevt = read.table(nflineevt, h=T, sep='\t')
spetree = read.tree(file=nfspetree)

event_text= c('Speciations', 'Transfers') ; names(event_text) = c('S', 'T')
metrics = c('nb_lineages', 'cum_freq')
metric_texts = c('- # gene lineages featuring event', '- cummulative event observation frequency') ; names(metric_texts) = metrics
recordons = c('rec_branch_name', 'don_branch_name')
recordon_texts = c('by recipient branch', 'by donor branch') ; names(recordon_texts) = recordons

plot_eventcounts = function(event_counts, spetree, main, log10cex=T){
	plot.spetree = plot.phylo(spetree, show.tip.label=F, show.node.label=T, main=main)
	ntips=length(spetree$tip.label)
	event_counts = event_counts[spetree$node.label]
	names(event_counts) = spetree$node.label
	event_counts[is.na(event_counts)] = 0
	log_event_counts = ceiling(log10(event_counts+1))
	if (log10cex) event_counts = log_event_counts
	nodelabels(node=ntips+(1:spetree$Nnode), pch=1, cex=event_counts, col=cols[log_event_counts])
	legend('topleft', legend=paste('1e', 0:(length(cols)-1), ' - 1e', 1:(length(cols)), sep=''), pch=1, col=cols, bg='white')	
}

cols = c('black', 'blue', 'green', 'gold', 'orange', 'red', 'purple')
pdf(paste(nflineevt, 'species_tree_density.pdf', sep='.'), height=18, width=32)
for (e in levels(lineevt$event_type)){
	print(e)
	for (metric in metrics){
		print(metric)
		if (e=='T'){
			for (recordon in recordons){
				print(recordon)
				events_brname = levels(lineevt[lineevt$event_type==e, recordon])
				event_counts = sapply(events_brname, function(x){
				 sum(lineevt[(lineevt$event_type==e & lineevt[[recordon]]==x), metric]) 
				})
				names(event_counts) = events_brname
				plot_eventcounts(event_counts, spetree, main=paste(event_text[e], metric_texts[metric], recordon_texts[recordon]))
			}
		}else{
			recordon='rec_branch_name'
			events_brname = lineevt[lineevt$event_type==e, recordon]
			event_counts = lineevt[lineevt$event_type==e, metric]
			names(event_counts) = events_brname
			plot_eventcounts(event_counts, spetree, main=paste(event_text[e], metric_texts[metric]))
		}
	}
}
dev.off()

# compute estimated time to process these event data
time1lineagepair = as.double(470.0)/79360180 # emprical estimate (in seconds) from parsing all lineage pairs for evtid=1795683

esttimeperevt = sapply(lineevt[['nb_lineages']], function(x){ x*(x-1)*time1lineagepair })

