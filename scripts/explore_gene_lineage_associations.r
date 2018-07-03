#!/usr/bin/Rscript --vanilla
library(gplots)
library(RColorBrewer)
library(igraph)
library(getopt)
library(parallel)
library(RPostgreSQL)

plotCoEvolQuantiles = function(qmatchev, detail.quant, sample.name=''){
	layout(matrix(c(1,1,1,2), 1, 4, byrow=T))
	plot(qmatchev[as.character(0:100/100)], xlab='quantiles', ylab='co-evolution score',
	 main=sprintf("distribution of co-evolution score\nin sample %s", sample.name))
	plot(qmatchev[c(detail.quant, 1.0)], xlab='quantiles')
}

parseMatchEventFile = function(nfmatchevents, quant.only=F, nmaxrlocsid=NULL, refrepliordlineages=NULL,
                                detail.quant=c(.999, .9999, .99999), top.quant.cutoff=.99999, top.val.cutoff=NULL,
                                comp.quant=T, plot.quant=F, verbose=F){
	# load (partial) gene lineage match report
	if (verbose) print(nfmatchevents, quote=F)
#~ 	matchev = unique(read.table(nfmatchevents))
	matchev = read.table(nfmatchevents, sep='\t', header=F, colClasses="numeric")
	colnames(matchev) = c('rlocdsid1', 'rlocdsid2', 'freq')
	p = unique(sort(c(0:99/100, detail.quant, top.quant.cutoff, 1.0)))
	ncomp = nrow(matchev)
	# every ~1.1Gb table has ~46.7M rows of reported pairwise associations, each using between 6 and 8 Gb of memory when loaded
	## empirical distrubution of co-evolution scores
	if (comp.quant){
		# compute
		qmatchev = quantile(matchev$freq, probs=p, names=F) ; names(qmatchev) = as.character(p)
		if (plot.quant){
			plotCoEvolQuantiles(qmatchev, detail.quant, sample.name=basename(nfmatchevents))
		}
	}else{
		qmatchev = NULL
	}
	qn = c(qmatchev, ncomp) ; names(qn) = c(names(qmatchev), 'nb.reported.comp')
	if (verbose) print(c(sprintf("%s ... done. quantiles:", nfmatchevents), qn) quote=F)
	# can stop here
	if (quant.only){ return( qn ) }
	## get set of unique query lineages and can derive the total count of comparisons that could have been attempted
	urlocdsid1 = unique(matchev$rlocdsid1)
	if (!is.null(nmaxrlocsid)){ totcomp = sum(sapply(urlocdsid1, function(k){ nmaxrlocsid - k }))
	}else{ totcomp = NULL }
	## reduce the table to top hits
	if (!is.null(top.val.cutoff)){
		topmatchev = matchev[matchev$freq>=top.val.cutoff,]
	}else{
		topmatchev = matchev[matchev$freq>=qmatchev[as.character(top.quant.cutoff)],]
	}
	## can project results on a reference replicon map (the resulting matrix can be summed with those from other partial results)
	if (!is.null(refrepliordlineages)){
		matmatchev = data.matrix(sapply(refrepliordlineages$rlocds_id, function(rlocdsid1){
			sapply(refrepliordlineages$rlocds_id, function(rlocdsid2){
				i = which((matchev$rlocdsid1==rlocdsid1 & matchev$rlocdsid2==rlocdsid2) | (matchev$rlocdsid1==rlocdsid2 & matchev$rlocdsid2==rlocdsid1))
				if (length(i)==1){ return(matchev$freq[i])
				}else{ if (length(i)==0){ return(NA) 
				}else{ print(matchev[i,]) ; stop(sprintf('multiple co-evolution score values for lineage pair %d, %d', rlocdsid1, rlocdsid2)) }}
			})
		}))
	}else{ matmatchev = NULL }
	if (verbose){ print(sprintf("%s ... done. Top co-evolution scores:", nfmatchevents), quote=F) ; print(head(topmatchev)) ; print(sprintf("... (%d lines)", nrow(topmatchev))) }
	return( list(quantiles=qmatchev, unique.query.rlocdsid=urlocdsid1, nb.reported.comp=ncomp, tot.comp=totcomp, top.matches=topmatchev, ref.repli.proj.mat=matmatchev) )
}

plotHMevents = function(mat, matname, cap=99){
#~ 	qmat = quantile(mat, (1:100)/100, na.rm=T)
#~ 	widthbin = qmat[cap]/8
#~ 	brk = c(seq(from=0, to=qmat[cap], by=widthbin), max(mat, na.rm=T))
	brk = quantile(mat, (0:9)/9, na.rm=T)
	print(brk)
	print(summary(mat))
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


#~ dirmatchevents = '/enterobac/PanteroDB_v0.3/07.compare_scenarios/bs_stem70_withinmedian35/replaceCCinGasinS-collapsePOPinSnotinG/ale_collapsed_undat/gene_event_matches'
#~ dirmatchevents = '/enterobac/PanteroDB_v0.3/07.compare_scenarios/bs_stem70_withinmedian35/replaceCCinGasinS-collapsePOPinSnotinG/ale_collapsed_undat/gene_event_matches_10253genes_06062018/plasmid2genefamily.mat_excl-transposases.NC_022651.1-centred_community_size13'
#~ refrepli = 'NC_022651.1'
#~ nfrefrepliordlineages = '/enterobac/PanteroDB_v0.3/07.compare_scenarios/bs_stem70_withinmedian35/replaceCCinGasinS-collapsePOPinSnotinG/ale_collapsed_undat/gene_event_matches_10253genes_06062018/plasmid2genefamily.mat_excl-transposases.NZ_CP008933.1-centred_community_size24/NC_022651.1_ordered_gene_lineage_products.tab'
#~ nftopmatchlinedetails = '/enterobac/PanteroDB_v0.3/07.compare_scenarios/bs_stem70_withinmedian35/replaceCCinGasinS-collapsePOPinSnotinG/ale_collapsed_undat/gene_event_matches_10253genes_06062018/covered_10253genes_lineages.details'
#~ nmaxrlocsid = 640228
#~ opt = list()

### argument (option) parsing
spec = matrix(c(
  'dir_match_events',   'i', 1, "character", "[mandatory] (preferably absolute) path of folder with matched events tables (output of compare_collapsedALE_scenarios.py, which names are expected to be starting with 'matching_events.tab')",
  'ref_replicons',      'r', 2, "character", "comma-separated list of replicons on which to project the LD-like matrix of co-evolution scores",
  'replicon_annot_map', 'm', 2, "character", "path to table file of gene lineage ordered as in reference replicons, with functional annotation of gene products; table must be tab-delimited and contain at least columns named: 'rlocds_id', 'cds_code', 'product'",
  'max_lineage_id',     'l', 2, "integer",   "maximum value if serial id of lineages to be compared; allows to compute the number of comparisons that were evaluated (usually much more than those reported due to filters on e.g. minimum co-evol scolution score number and correct quantiles accordingly)",
  'output_dir',         'o', 2, "character", "path to folder where to write output files; default to the parent folder of input match event folder",
  'output_prefix',      'p', 2, "character", "prefix for output files; default to the name of input match event folder",
  'quant_cutoff',       'q', 2, "double",    "cut-off fraction of co-evolution score distribution over which to to report top gene lineage associations",
  'score_cutoff',       's', 2, "double",    "cut-off co-evolution score over which to to report top gene lineage associations (bypass quantile computation)",
  'load_quantiles',     'Q', 0, "character", "will attempt to load quantile information from R object archive file saved in a previous run; file must be named like `dir_match_events` argument value with prefix '.co-evolution_quantiles.RData', and contain a data.frame bearing the name `matchevdf` (bypass quantile computation)",
  'db_name',            'D', 2, "character", "name of database containing detailed gene annotation information, to be collated for top asocciated gene lineages",
  'db_type',            'T', 2, "character", "SQL database type; default to PostgreSQL",
  'db_host',            'H', 2, "character", "SQL database host; default to 'localhost' (only relevant for DB with client/server system, e.g. PostgreSQL, but not for file-based systems like SQLite)",
  'db_pw',              'W', 2, "character", "SQL database password; default to none",
  'db_user',            'U', 2, "character", "SQL database user; default to system user executing this script",
  'db_port',            'P', 2, "character", "SQL database port; default to 5432",
  'threads',         't', 2, "integere",  "number of parrallel threads for input data processing (each process may use from 6 to 10 Gb when loading a 1 Gb table file)"
), byrow=TRUE, ncol=5);
opt = getopt(spec, opt=commandArgs(trailingOnly=T))

if ( is.null(opt$dir_match_events  ) ){ 
	stop("path to input file folder required! please specify it through '--dir_match_events' option")  
}else{ 
	dirmatchevents = gsub('/$', '', opt$dir_match_events) 
}
patnfmatchevents = file.path(dirmatchevents, 'matching_events.tab*')
#~ print(patnfmatchevents)
lnfmatchevents = Sys.glob(patnfmatchevents)
#~ print(head(lnfmatchevents))
if (length(lnfmatchevents)>1){
	lnfmatchevents = lnfmatchevents[order(as.numeric((sapply(strsplit(lnfmatchevents, split='\\.'), function(x){ x[length(x)] }))))]
	print(sprintf("will load lineage co-evolution scores from %d files matchin '%s'", length(lnfmatchevents), patnfmatchevents), quote=F)
}
if ( is.null(opt$output_dir) ){ 
	dirout = dirname(dirmatchevents)  
}else{
	dirout = opt$output_dir
}
if ( is.null(opt$output_repfix) ){ 
	prefixout = basename(dirmatchevents)  
}else{
	prefixout = opt$output_repfix
}
qcutoff = scutoff = NULL
if (is.null(opt$score_cutoff)){
	if (is.null(opt$quant_cutoff)){
		qcutoff = .99999
	}else{
		qcutoff = opt$quant_cutoff
	}
}else{
	scutoff = opt$score_cutoff
}
refrepli = opt$ref_replicons
nmaxrlocsid = opt$max_lineage_id
nfrefrepliordlineages = opt$replicon_annot_map
if (!is.null(opt$replicon_annot_map)){
	refrepliordlineages = read.table(opt$replicon_annot_map, h=T, sep='\t')
}else{
	refrepliordlineages = NULL
}
if (is.null(opt$threads)){
	ncores = 6
}else{
	ncores = opt$threads
}

if (!is.null(qcutoff)){
	nfsavematchev = paste(file.path(dirout, prefixout), 'co-evolution_quantiles.RData', sep='.')
	if (opt$load_quantiles & file.exists(nfsavematchev)){
		load(nfsavematchev)
		print(sprintf("loaded quantiles of co-evolution score from pre-exesting file: '%s", nfsavematchev), quote=F)
	}else{
		### first pass: parse data to extract co-evolution score quantiles
		print("will load lineage co-evolution score tables a first time to compute quantiles", quote=F)
		lparsematchev = mclapply(lnfmatchevents, parseMatchEventFile, quant.only=T, top.quant.cutoff=qcutoff, mc.cores=ncores)
		names(lparsematchev) = basename(lnfmatchevents)

		matchevdf = as.data.frame(t(simplify2array(lparsematchev)))

		save(matchevdf, file=nfsavematchev)
		print(sprintf("saved quantiles of co-evolution score to file: '%s", nfsavematchev), quote=F)
	}
	pdf(paste(file.path(dirout, prefixout), 'associations_perfile.pdf', sep='.'), height=10, width=10)
	layout(matrix(c(1,1,1,2), 1, 4, byrow=T))
	q1 = as.character(0:100/100)
	boxplot(lapply(q1, function(j){ matchevdf[,j] }), names=q1, xlab='quantiles', 'co-evolution score',
	 main="distribution of co-evolution score\nacross samples")
	q2 = colnames(matchevdf)[!(colnames(matchevdf) %in% c(q1, 'nb.reported.comp'))]
	boxplot(lapply(q2, function(j){ matchevdf[,j] }), names=q2, xlab='quantiles')
	dev.off()

	totrepcomp = sum(matchevdf[,'nb.reported.comp'])
	avg.cutoff.val = sum(matchevdf[,as.character(qcutoff)]*matchevdf[,'nb.reported.comp'])/totrepcomp

#~ 	save(matchevdf, totrepcomp, avg.cutoff.val, file=paste(file.path(dirout, prefixout), 'co-evolution_quantiles.RData', sep='.'))

	print(sprintf("found an average value of %f for quantiles at p=%f of empircal distribution of co-evolution scores over %g reported gene lineage comparisons (stored within %d separate files)",
				   avg.cutoff.val, qcutoff, totrepcomp, length()), quote=F)
	scutoff = avg.cutoff.val
}

### second pass: parse data to extract top co-evolution scores
lparsematchev = mclapply(lnfmatchevents, parseMatchEventFile, comp.quant=F, top.val.cutoff=avg.cutoff.val, 
                         refrepliordlineages=refrepliordlineages, mc.cores=ncores)
save(lparsematchev, file=paste(file.path(dirout, prefixout), 'top_association.RData', sep='.'))
topmatchev = lparsematchev[[1]]$top.matches
for (i in 2:length(lparsematchev)){
	topmatchev = rbind(topmatchev, lparsematchev[[i]]$top.matches)
}
write.table(topmatchev, file=paste(file.path(dirout, prefixout), 'top_associations.tab', sep='.'), sep='\t', row.names=F)

# create network of gene lineage association
graph_topmatchev = graph_from_data_frame(topmatchev, directed=F)
save(graph_topmatchev, file=paste(file.path(dirout, prefixout), 'top_association_network.RData', sep='.'))
pdf(file=paste(file.path(dirout, prefixout), 'top_association_network.pdf', sep='.'), height=30, width=30)
plot(graph_topmatchev)
dev.off()

# collect detailed annotation data on top associated lineages
if (!is.null(opt$db_name)){
	if (is.null(opt$db_type){ dbtype = "PostgreSQL" }else{ dbtype = opt$db_type }
	if (grepl('postgres', dbtype, ignore.case=T){
		library("RPostgreSQL")
		dbdrv = dbDriver("PostgreSQL")
		dbtype = "PostgreSQL"
	}else{ if (grepl('sqlite', dbtype, ignore.case=T)){
		library("RSQLite")
		dbdrv = dbDriver("SQLite")
		dbtype = "SQLite"
	}else{ stop(sprintf("database system not supported: %s; please use either 'PostgreSQL' or 'SQLite'", dbtype)) }}
	dbcon = dbConnect(dbdrv, 
					  dbname   = opt$db_name,
					  host     = ifelse(!is.null(opt$db_host), opt$db_host, "localhost"), 
					  port     = ifelse(!is.null(opt$db_port), opt$db_port, 5432), 
					  password = ifelse(!is.null(opt$db_password), opt$db_password, ''), 
					  user     = ifelse(!is.null(opt$db_user), opt$db_user, Sys.info()["user"]))
	u1 = unique(topmatchev$rlocdsid1)
	u2 = unique(topmatchev$rlocdsid2)
	u = intersect(u1, u2)
	if (dbtype=="PostgreSQL"){ dbExecute(dbcon, "set search_path = genome , phylogeny;") }
	
	dbExecute(dbcon, "create temp table toplineages (rlocds_id int); ")
	dbWriteTable(dbcon, "toplineages", value = data.frame(rlocds_id=u), append=T, row.names=F)
	
	topmatchlinedetails = dbGetQuery(dbcon, "
    select distinct rlocds_id, replacement_label_or_cds_code, max(product) as product, max(cds_code) as repr_cds_code
      from toplineages
      inner join replacement_label_or_cds_code2gene_families using (rlocds_id)
      inner join gene_tree_label2cds_code using (replacement_label_or_cds_code)
      inner join coding_sequences using (cds_code)
      inner join proteins using (genbank_nr_protein_id)
    group by rlocds_id, replacement_label_or_cds_code order by repr_cds_code ;")
    
    dbExecute(dbcon, "drop table toplineages; ")
	
	topmatchevdetail = merge(merge(topmatchev, topmatchlinedetails, by.x='rlocdsid1', by.y='rlocds_id'), topmatchlinedetails, by.x='rlocdsid2', by.y='rlocds_id', suffixes=1:2)
	topmatchevdetail = topmatchevdetail[order(topmatchevdetail$repr_cds_code1, topmatchevdetail$repr_cds_code2),c('rlocdsid1', 'rlocdsid2', 'freq', 'replacement_label_or_cds_code1', 'replacement_label_or_cds_code2', 'product1', 'product2', 'repr_cds_code1', 'repr_cds_code2')]
	topmatchevdetail$species_pop1 = as.factor(sapply(strsplit(as.character(topmatchevdetail[,'replacement_label_or_cds_code1']), split='_'), `[`, 1))
	topmatchevdetail$species_pop2 = as.factor(sapply(strsplit(as.character(topmatchevdetail[,'replacement_label_or_cds_code2']), split='_'), `[`, 1))
	genome_cdspos1 = as.data.frame(t(simplify2array(strsplit(as.character(topmatchevdetail[,'repr_cds_code1']), split='_')))) ; colnames(genome_cdspos1) = c('genome1', 'cds_pos1')
	genome_cdspos2 = as.data.frame(t(simplify2array(strsplit(as.character(topmatchevdetail[,'repr_cds_code2']), split='_')))) ; colnames(genome_cdspos2) = c('genome2', 'cds_pos2')
	topmatchevdetail = cbind(topmatchevdetail, genome_cdspos1, genome_cdspos2)
	write.table(topmatchevdetail, file=paste(file.path(dirout, prefixout), 'top_associations.detail.tab', sep='.'), sep='\t', row.names=F)
	
	samegenometopmatches = topmatchevdetail[as.character(topmatchevdetail$genome1)==as.character(topmatchevdetail$genome2),]
	pdf(paste(file.path(dirout, prefixout), 'same_genome_top_associations.jointfreq_vs_gene_dist.pdf', sep='.'), height=10, width=15)
	plot(abs(as.numeric(samegenometopmatches$cds_pos1) - as.numeric(samegenometopmatches$cds_pos2)), samegenometopmatches$freq, xlab='gene-gene distance', ylab='Co-evolution score', main='Co-evolution score and inter-gene distance\n(lineages leading to genes residing in the same genome)')
	dev.off()
	
	graph_topmatchevdet = graph_from_data_frame(topmatchevdetail[,c('repr_cds_code1', 'repr_cds_code2')], directed=F)
	graph_topmatchevdet$color = rainbow(nlevels(topmatchevdetail$species_pop1))[as.integer(topmatchevdetail$species_pop1)]
	pdf(paste(file.path(dirout, prefixout), 'top_associations_representative_host_colours.pdf', sep='.'), height=10, width=10)
	plot(graph_topmatchevdet, vertex.color=graph_topmatchevdet$color, vertex.label=NA)
	dev.off()
}


if (!is.null(refrepliordlineages)){
	matmatchev = lparsematchev[[1]]$ref.repli.proj.mat
	for (i in 2:length(lparsematchev)){
		matmatchev = matmatchev + lparsematchev[[i]]$ref.repli.proj.mat
	}
	nfpdf = paste(file.path(dirout, prefixout), sprintf('co-evol_scores_projection_%s.pdf', refrepli), sep='.')
	pdf(nfpdf, height=20, width=20)
	plotHMevents(matmatchev, sprintf('co-evolution scores projected on %s map', refrepli))
	dev.off()
	print(nfpdf)
}
