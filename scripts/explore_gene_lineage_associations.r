#!/usr/bin/Rscript --vanilla
library(gplots)
library(RColorBrewer)
library(igraph)
library(getopt)
library(parallel)
library(RPostgreSQL)

eventstablim = 50
plotlabellim = 100

repli.invar = c("assembly_id", "genomic_accession", "replicon_type", "replicon_size")
comm.algos =  c("louvain", "fast_greedy", "infomap", "label_prop")

outannotcols = c('freq', sapply(1:2, function(i){ 
	paste0( c("rlocdsid","replacement_label_or_cds_code","product","repr_cds_code","gene_family_id","tot_lineage_freq","chromorplas","species_pop","singmulti","coreoracc"), i)
}))

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
	if (verbose) print(sprintf("loading %s ...", nfmatchevents), quote=F)
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
	if (verbose) print(c(sprintf("%s ... done. quantiles:", nfmatchevents), qn), quote=F)
	# can stop here
	if (quant.only){ return( qn ) }
	## get set of unique query lineages and can derive the total count of comparisons that could have been attempted
	urlocdsid1 = unique(matchev$rlocdsid1)
	if (!is.null(nmaxrlocsid)){ totcomp = sum(sapply(urlocdsid1, function(k){ nmaxrlocsid - k }))
	}else{ totcomp = NULL }
	## reduce the table to top hits
	if (!is.null(top.val.cutoff)){
		topmatchev = matchev[matchev$freq>=top.val.cutoff,]
		topstr = sprintf("scores >= %.2g", top.val.cutoff)
	}else{
		topmatchev = matchev[matchev$freq>=qmatchev[as.character(top.quant.cutoff)],]
		topstr = sprintf("top %.g", 1-top.quant.cutoff)
	}
	if (verbose){ print(sprintf("%s ... done. Top co-evolution scores (%s):", nfmatchevents, topstr), quote=F) ; print(head(topmatchev)) ; print(sprintf("... (%d lines)", nrow(topmatchev))) }
	## can project results on a reference replicon map (the resulting matrix can be summed with those from other partial results)
	if (!is.null(refrepliordlineages)){
		if (verbose){ print("build matrix of gene-to-gene coevolution scoore based on selected replicon map", quote=F) }
		matmatchev = data.matrix(sapply(refrepliordlineages$rlocds_id, function(rlocdsid1){
			i1 = which(matchev$rlocdsid1==rlocdsid1)
			i2 = which(matchev$rlocdsid2==rlocdsid1)
			rowmat = sapply(refrepliordlineages$rlocds_id, function(rlocdsid2){
				j1 = which(matchev$rlocdsid1[i2]==rlocdsid2)
				j2 = which(matchev$rlocdsid2[i1]==rlocdsid2)
				k = c(i1[j2], i2[j1])
				if (length(k)==1){ return(matchev$freq[k])
				}else{ if (length(k)==0){ return(NA) 
				}else{ print(matchev[k,]) ; stop(sprintf("multiple co-evolution score values for lineage pair %d, %d", rlocdsid1, rlocdsid2)) }}
			})
			if (verbose) cat(rlocdsid1, " ")
			return(rowmat) 
		}))
		if (verbose) cat("... done.\n")
	}else{ matmatchev = NULL }
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

connectionDB = function(opt){
	if (is.null(opt$db_type)){ dbtype = "PostgreSQL" }else{ dbtype = opt$db_type }
	if (grepl('postgres', dbtype, ignore.case=T)){
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
	if (dbtype=="PostgreSQL"){ dbExecute(dbcon, "set search_path = genome , phylogeny;") }	
	return(dbcon)
}

computeCommunities = function(graph, method="louvain"){
	if (method=="louvain") return( cluster_louvain(graph, weight=graph$freq) )
	if (method=="fast_greedy") return( cluster_fast_greedy(graph, weight=graph$freq, merges=T, membership=T) )
	if (method=="edge_betweenness") return( cluster_edge_betweenness(graph, weight=graph$freq, merges=T, membership=T) )
	if (method=="infomap") return( cluster_infomap(graph, e.weight=graph$freq) )
	if (method=="label_prop") return( cluster_label_prop(graph, weight=graph$freq) )
#~ 	stop("incorrect method specified")
}

#~ dirmatchevents = '/enterobac/PanteroDB_v0.3/07.compare_scenarios/bs_stem70_withinmedian35/replaceCCinGasinS-collapsePOPinSnotinG/ale_collapsed_undat/gene_event_matches'
#~ dirmatchevents = '/enterobac/PanteroDB_v0.3/07.compare_scenarios/bs_stem70_withinmedian35/replaceCCinGasinS-collapsePOPinSnotinG/ale_collapsed_undat/gene_event_matches_10253genes_06062018/plasmid2genefamily.mat_excl-transposases.NC_022651.1-centred_community_size13'
#~ refrepli = 'NC_022651.1'
#~ nfrefrepliordlineages = '/enterobac/PanteroDB_v0.3/07.compare_scenarios/bs_stem70_withinmedian35/replaceCCinGasinS-collapsePOPinSnotinG/ale_collapsed_undat/gene_event_matches_10253genes_06062018/plasmid2genefamily.mat_excl-transposases.NZ_CP008933.1-centred_community_size24/NC_022651.1_ordered_gene_lineage_products.tab'
#~ nftoplinereprannot = '/enterobac/PanteroDB_v0.3/07.compare_scenarios/bs_stem70_withinmedian35/replaceCCinGasinS-collapsePOPinSnotinG/ale_collapsed_undat/gene_event_matches_10253genes_06062018/covered_10253genes_lineages.details'
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
  'db_name',            'D', 2, "character", "name of database containing detailed gene annotation information, to be collated for top asocciated gene lineages",
  'db_type',            'T', 2, "character", "SQL database type; default to PostgreSQL",
  'db_host',            'H', 2, "character", "SQL database host; default to 'localhost' (only relevant for DB with client/server system, e.g. PostgreSQL, but not for file-based systems like SQLite)",
  'db_pw',              'W', 2, "character", "SQL database password; default to none",
  'db_user',            'U', 2, "character", "SQL database user; default to system user executing this script",
  'db_port',            'P', 2, "character", "SQL database port; default to 5432",
  'threads',            't', 2, "integer",   "number of parrallel threads for input data processing (each process may use from 6 to 10 Gb when loading a 1 Gb table file)",
  'verbose',            'v', 0, "logical",   "verbose mode",
  'parse_max_files',    'X', 2, "integer",   "maximum number of files from input list that will actually be parsed (in order from the list)",
  'load_quantiles',     'Q', 0, "logical",   "will attempt to load quantile information from R object archive file saved in a previous run; file must be named like `dir_match_events` argument value with prefix '.co-evolution_quantiles.RData', and contain a data.frame named `matchevdf` (bypass quantile computation)",
  'load_top_ass',       'A', 0, "logical",   "will attempt to load top lineage association table from R object archive file saved in a previous run; file must be named like `dir_match_events` argument value with prefix '.top_association.RData', and contain a list named `lparsematchev` (bypass quantile computation and top association retrieval)",
  'load_top_network',   'N', 0, "logical",   "will attempt to load top lineage association network from R object archive file saved in a previous run; file must be named like `dir_match_events` argument value with prefix '.top_association_network.RData', and contain a list named `lparsematchev` (bypass quantile computation, top association retrieval and graph computation)",
  'load_up_to_step',    'L', 2, "integer",   "will attempt to load data from R object archive file saved during a previous run, up to the desired step (1:quantile computation; 2:top association retrieval; 3:top association graph computation; 4:top associated lineage annotation retrival from SQL db)",
  'plot_network',       'n', 0, "logical",   "enable plotting of the networks (long)",
  'stop_after_step',    'Z', 2, "integer",   "will run until reaching the desired progress of computation, then quit (steps: quantile computation=1; top association retrieval=2; top association graph computation=3; top associated lineage annotation retrival from SQL db=4; compute association network=5; search communities in network=6; plot top association network with metadata=7; plot heatmaps of lineages association intensity as projected on maps of selected replicons [end])"
), byrow=TRUE, ncol=5);
opt = getopt(spec, opt=commandArgs(trailingOnly=T))

print("options:", quote=F)
print(opt)

if (!is.null(opt$verbose)){ verbose = TRUE }else{ verbose = FALSE }
if (!is.null(opt$plot_network)){ plotnw = TRUE }else{ plotnw = FALSE }

if ( is.null(opt$dir_match_events  ) ){ 
	stop("path to input file folder required! please specify it through '--dir_match_events' option")  
}else{ 
	dirmatchevents = gsub('/$', '', opt$dir_match_events) 
}
patnfmatchevents = file.path(dirmatchevents, 'matching_events.tab*')
#~ print(patnfmatchevents)
lnfmatchevents = Sys.glob(patnfmatchevents)
if (!is.null(opt$parse_max_files)){
	lnfmatchevents = lnfmatchevents[1:opt$parse_max_files]
}
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
if ( is.null(opt$output_prefix) ){ 
	prefixout = basename(dirmatchevents)  
}else{
	prefixout = opt$output_prefix
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
	refrepliordlineages = read.table(opt$replicon_annot_map, h=T, sep='\t', quote='')
}else{
	refrepliordlineages = NULL
}
if (is.null(opt$threads)){
	ncores = 6
}else{
	ncores = opt$threads
}
verbose = opt$verbose
if (is.null(opt$stop_after_step)){
	endscript = 10
}else{
	endscript = opt$stop_after_step
}
if (is.null(opt$load_up_to_step)){
	loadup = 0
}else{
	loadup = opt$load_up_to_step
}

dbcon = NULL
if (is.null(opt$db_type)){
	dbtype = "PostgreSQL" 
}else{ 
	dbtype = opt$db_type 
}

if (!is.null(qcutoff)){
	nfsavematchev = paste(file.path(dirout, prefixout), 'co-evolution_quantiles.RData', sep='.')
	if ((!is.null(opt$load_quantiles) | loadup>=1) & file.exists(nfsavematchev)){
		load(nfsavematchev)
		if (verbose) print(sprintf("loaded quantiles of co-evolution score from pre-exesting file: '%s", nfsavematchev), quote=F)
	}else{
		### step 1: parse data to extract co-evolution score quantiles
		if (verbose) print("will read lineage co-evolution score tables a first time to compute quantiles", quote=F)
		lparsematchev = mclapply(lnfmatchevents, parseMatchEventFile, quant.only=T, top.quant.cutoff=qcutoff, verbose=verbose, mc.cores=ncores)
		names(lparsematchev) = basename(lnfmatchevents)

		matchevdf = as.data.frame(t(simplify2array(lparsematchev)))

		save(matchevdf, file=nfsavematchev)
		if (verbose) print(sprintf("saved quantiles of co-evolution score to file: '%s", nfsavematchev), quote=F)
	}
	q1 = as.character(0:100/100)
	q2 = colnames(matchevdf)[!(colnames(matchevdf) %in% c(q1, 'nb.reported.comp'))]
	pdf(paste(file.path(dirout, prefixout), 'co-evolution_quantiles.pdf', sep='.'), height=6, width=10)
	layout(matrix(c(1,1,1,2), 1, 4, byrow=T))
	plot(y=quantile(unlist(matchevdf[,q1]), probs=as.numeric(q1), names=F), x=q1, xlab='quantiles', ylab='co-evolution score',
	 main="distribution of co-evolution scores")
	plot(y=apply(matchevdf[,q2], 2, mean), x=q2, xlab='quantiles', ylab='co-evolution score')
	dev.off()
	pdf(paste(file.path(dirout, prefixout), 'co-evolution_quantiles_perfile.pdf', sep='.'), height=6, width=10)
	layout(matrix(c(1,1,1,2), 1, 4, byrow=T))
	boxplot(lapply(q1, function(j){ matchevdf[,j] }), names=q1, xlab='quantiles', ylab='co-evolution score',
	 main="distribution of co-evolution scores\nacross data partitions")
	boxplot(lapply(q2, function(j){ matchevdf[,j] }), names=q2, xlab='quantiles')
	dev.off()
	

	totrepcomp = sum(matchevdf[,'nb.reported.comp'])
	avg.cutoff.val = sum(matchevdf[,as.character(qcutoff)]*matchevdf[,'nb.reported.comp'])/totrepcomp

#~ 	save(matchevdf, totrepcomp, avg.cutoff.val, file=paste(file.path(dirout, prefixout), 'co-evolution_quantiles.RData', sep='.'))

	print(sprintf("found an average value of %f for quantiles at p=%f of empircal distribution of co-evolution scores over %g reported gene lineage comparisons (stored within %d separate files)",
				   avg.cutoff.val, qcutoff, totrepcomp, length(lnfmatchevents)), quote=F)
	scutoff = avg.cutoff.val
}
if (endscript<=1){ quit(save='no') }

### step 2: parse data to extract top co-evolution scores
nftopass = paste(file.path(dirout, prefixout), 'top_associations.RData', sep='.')
nftopasstab = paste(file.path(dirout, prefixout), 'top_associations.tab', sep='.')
if ((!is.null(opt$load_top_ass) | loadup>=2) & file.exists(nftopasstab) & is.null(refrepliordlineages)){
	topmatchev = read.table(nftopasstab, sep='\t', header=T)
	if (verbose) print(sprintf("loaded top association data from file: '%s'", nftopasstab))
}else{
	if ((!is.null(opt$load_top_ass) | loadup>=2) & file.exists(nftopass)){
	load(nftopass)
	if (verbose) print(sprintf("loaded top association data from file: '%s'", nftopass))
	}else{
		if (verbose) print("extract top associations from bulk data")
		lparsematchev = mclapply(lnfmatchevents, parseMatchEventFile, comp.quant=F, top.val.cutoff=scutoff, 
								 refrepliordlineages=refrepliordlineages, verbose=verbose, mc.cores=ncores)
		save(lparsematchev, file=nftopass)
	}
	topmatchev = lparsematchev[[1]]$top.matches
	if (!is.null(refrepliordlineages)){ matmatchev = lparsematchev[[1]]$ref.repli.proj.mat }
	if (length(lparsematchev)>1){
		for (i in 2:length(lparsematchev)){
			topmatchev = rbind(topmatchev, lparsematchev[[i]]$top.matches)
			if (!is.null(refrepliordlineages)){ matmatchev = matmatchev + lparsematchev[[i]]$ref.repli.proj.mat }
		}
	}
	rm(lparsematchev) ; gc()
	write.table(topmatchev, file=nftopasstab, sep='\t', row.names=F)
}
if (endscript<=2){ quit(save='no') }

### step 3: create network of gene lineage association
nftopassgraph = paste(file.path(dirout, prefixout), 'top_association_network.RData', sep='.')
if ((!is.null(opt$load_top_network) | loadup>=3) & file.exists(nftopassgraph)){
	load(nftopassgraph)
	if (verbose) print(sprintf("loaded network of top associated gene lineages from file: '%s'", nftopassgraph))
}else{
	if (verbose) print("compute network of top associated gene lineages")
	graph_topmatchev = graph_from_data_frame(topmatchev, directed=F)
	save(graph_topmatchev, file=nftopassgraph)
}
#~ pdf(file=paste(file.path(dirout, prefixout), 'top_association_network.pdf', sep='.'), height=30, width=30)
#~ plot(graph_topmatchev)
#~ dev.off()
if (endscript<=3){ quit(save='no') }

### step 4: collect detailed annotation data on top associated lineages
if (!is.null(opt$db_name)){
	nftopmatchevdetail = paste(file.path(dirout, prefixout), 'top_associations_annot.RData', sep='.')
	if ((!is.null(opt$load_top_annot) | loadup>=4) & file.exists(nftopmatchevdetail)){
		load(nftopmatchevdetail)
	if (verbose) print(sprintf("loaded annotation of (genes in) top associated gene lineages from file: '%s'", nftopmatchevdetail))
	}else{
		if (verbose) print(sprintf("retrieve annotation of (genes in) top associated gene lineages from %s database: '%s'", dbtype, opt$db_name))
		
		if (is.null(dbcon)) dbcon = connectionDB(opt)
		
		u1 = unique(topmatchev$rlocdsid1)
		u2 = unique(topmatchev$rlocdsid2)
		u = union(u1, u2)
		
		dbExecute(dbcon, "create temp table toplineages (rlocds_id int); ")
		dbWriteTable(dbcon, "toplineages", value = data.frame(rlocds_id=u), append=T, row.names=F)
		
		# query annotations - 1 row / lineage, with a representative CDS
		annotreprquery = "
		select distinct rlocds_id, replacement_label_or_cds_code, max(product) as product, max(cds_code) as repr_cds_code,
		 count(cds_code) as tot_lineage_freq, count(distinct assembly_id) as lineage_pres_freq, gene_family_id, size as fam_size, genome_present as fam_occurence,
		 sum(case when replicon_type='chromosome' then 1 else 0 end) as count_chromosome, sum(case when replicon_type='plasmid' then 1 else 0 end) as count_plasmid
		  from toplineages
		  inner join replacement_label_or_cds_code2gene_families using (rlocds_id)
		  inner join gene_tree_label2cds_code using (replacement_label_or_cds_code)
		  inner join coding_sequences using (cds_code, gene_family_id)
		  inner join proteins using (genbank_nr_protein_id)
		  inner join replicons using (genomic_accession)
		  inner join gene_family_sizes using (gene_family_id)
		group by rlocds_id, replacement_label_or_cds_code, gene_family_id, size, genome_present 
		order by repr_cds_code ;
		"
		if (verbose) cat(annotreprquery)
		toplinereprannot = dbGetQuery(dbcon, annotreprquery)
		toplinereprannot$gene_family_id = as.factor(trimws(as.character(toplinereprannot$gene_family_id)))
		nb_genome = max(toplinereprannot$fam_occurence)
		toplinereprannot$chromorplas = as.factor(ifelse(toplinereprannot$count_chromosome>0, ifelse(toplinereprannot$count_plasmid==0, 'chromosome', 'both'), ifelse(toplinereprannot$count_plasmid>0, 'plasmid', 'none')))
		toplinereprannot$species_pop = as.factor(sapply(strsplit(as.character(toplinereprannot[,'replacement_label_or_cds_code']), split='_'), `[`, 1))
		toplinereprannot$coreoracc = ifelse(toplinereprannot$fam_occurence >= nb_genome*0.98, "core", "accessory")
		toplinereprannot$singmulti = ifelse((toplinereprannot$fam_size - toplinereprannot$fam_occurence)/toplinereprannot$fam_occurence <= 0.10, "single", "multi")
		write.table(toplinereprannot, file=paste(file.path(dirout, prefixout), 'top_lineage_annot.tab', sep='.'), sep='\t', row.names=F)

		topmatchevdetail = merge(merge(topmatchev, toplinereprannot, by.x='rlocdsid1', by.y='rlocds_id'), toplinereprannot, by.x='rlocdsid2', by.y='rlocds_id', suffixes=1:2)
		topmatchevdetail = topmatchevdetail[order(topmatchevdetail$repr_cds_code1, topmatchevdetail$repr_cds_code2), outannotcols]
		write.table(topmatchevdetail, file=paste(file.path(dirout, prefixout), 'top_associations_repr_gene_annot.tab', sep='.'), sep='\t', row.names=F)
		
		# query annotations - for all CDS in the lineages
		
		multiannotquery = "
		select distinct rlocds_id, replacement_label_or_cds_code, cds_code, product,
		  assembly_id, genomic_accession, replicon_type, replicon_size, gene_family_id, cds_begin
		  from toplineages
		  inner join replacement_label_or_cds_code2gene_families using (rlocds_id)
		  inner join gene_tree_label2cds_code using (replacement_label_or_cds_code)
		  inner join coding_sequences using (gene_family_id, cds_code)
		  inner join proteins using (genbank_nr_protein_id)
		  inner join replicons using (genomic_accession);
		  "
		if (verbose) cat(multiannotquery)
		toplinemultiannot = dbGetQuery(dbcon, multiannotquery)
		dbExecute(dbcon, "drop table toplineages; ")
		toplinemultiannot$assembly_id = as.factor(trimws(as.character(toplinemultiannot$assembly_id)))
		toplinemultiannot$replicon_type = as.factor(trimws(as.character(toplinemultiannot$replicon_type)))
		toplinemultiannot$gene_family_id = as.factor(trimws(as.character(toplinemultiannot$gene_family_id)))
		
		# add the details of all genes represented in the lineages, but keep only rows where the associated genes are on the same replicon.
		topmatchevmultidetail = merge(merge(topmatchev, toplinemultiannot, by.x='rlocdsid1', by.y='rlocds_id'),
		                              toplinemultiannot, by.x=c('rlocdsid2', repli.invar), by.y=c('rlocds_id', repli.invar), suffixes=1:2)
		topmatchevmultidetail = topmatchevmultidetail[order(topmatchevmultidetail$cds_code1, topmatchevmultidetail$cds_code2), c(repli.invar, sort(setdiff(colnames(topmatchevmultidetail), repli.invar)))]
		topmatchevmultidetail$physical_dist = apply(topmatchevmultidetail[,c('cds_begin1', 'cds_begin2', 'replicon_size')], 1, function(x){
			if (x[1] < x[2]){ a = x[1] ; b = x[2] }else{ b = x[1] ; a = x[2] }
			min((b - a), (x[3] - b + a))
		})
		write.table(topmatchevmultidetail, file=paste(file.path(dirout, prefixout), 'top_associations_same-replicon_gene_annot.tab', sep='.'), sep='\t', row.names=F)
		
		# plot relation between co-evolution scaore of top associated lineages and the physical or locus distance of their representatives on replicons
		pdf(paste(file.path(dirout, prefixout), 'top_associations_same-replicon.jointfreq_vs_physical_dist.pdf', sep='.'), height=10, width=15)
		plot(topmatchevmultidetail$physical_dist, topmatchevmultidetail$freq, xlab='physical distance', ylab='Co-evolution score', main='Co-evolution score and inter-gene distance\n(lineages leading to genes residing in the same genome)')
		lmgenedistfreq = lm(topmatchevmultidetail$freq ~ topmatchevmultidetail$physical_dist)
		abline(lmgenedistfreq, col='red')
		print(summary(lmgenedistfreq))
		dev.off()
		
		save(toplinereprannot, topmatchevdetail, toplinemultiannot, topmatchevmultidetail, lmgenedistfreq, file=nftopmatchevdetail)
	}
}
if (endscript<=4){ quit(save='no') }

nftopassannotgraph = paste(file.path(dirout, prefixout), 'top_association_annotated_network.RData', sep='.')
if ((!is.null(opt$load_annot_graph) | loadup>=5) & file.exists(nftopassannotgraph)){
	load(nftopassannotgraph)
	if (verbose) print(sprintf("loaded annotated version of network of top associated gene lineages from file: '%s'", nftopassannotgraph))
}else{
	### step 5: compute and plot network of top associated lineages along with associated metadata
	if (verbose) print("annotate network of top associated gene lineages from file and plot it in full")
	vertex_pop = as.factor(sapply(as.numeric(names(V(graph_topmatchev))), function(rlocdsid){
		toplinereprannot[toplinereprannot$rlocds_id==rlocdsid,'species_pop']
	}))
	vertex_rep = as.factor(sapply(as.numeric(names(V(graph_topmatchev))), function(rlocdsid){
		toplinereprannot[toplinereprannot$rlocds_id==rlocdsid,'chromorplas']
	}))
	vertex_fam = as.factor(sapply(as.numeric(names(V(graph_topmatchev))), function(rlocdsid){
		toplinereprannot[toplinereprannot$rlocds_id==rlocdsid,'gene_family_id']
	}))
	vertex_famcopy = as.factor(sapply(as.numeric(names(V(graph_topmatchev))), function(rlocdsid){
		toplinereprannot[toplinereprannot$rlocds_id==rlocdsid,'singmulti']
	}))
	vertex_fampan = as.factor(sapply(as.numeric(names(V(graph_topmatchev))), function(rlocdsid){
		toplinereprannot[toplinereprannot$rlocds_id==rlocdsid,'coreoracc']
	}))
	colorpop = rainbow(nlevels(vertex_pop)) ; names(colorpop) = levels(vertex_pop)
	shaperep = c("circle", "square", "csquare") ; names(shaperep) = c("chromosome", "plasmid", "both")
	colorpan = c("black", "red") ; names(colorpan) = c("core", "accessory")
	graph_topmatchev = set_vertex_attr(graph_topmatchev, "spe.pop", value=as.character(vertex_pop))
	graph_topmatchev = set_vertex_attr(graph_topmatchev, "color.pop", value=colorpop[as.character(vertex_pop)])
	graph_topmatchev = set_vertex_attr(graph_topmatchev, "repli.type", value=as.character(vertex_rep))
	graph_topmatchev = set_vertex_attr(graph_topmatchev, "shape.rep", value=shaperep[as.character(vertex_rep)])
	graph_topmatchev = set_vertex_attr(graph_topmatchev, "fam", value=as.character(vertex_fam))
	graph_topmatchev = set_vertex_attr(graph_topmatchev, "pang.fam.status", value=as.character(vertex_fampan))
	graph_topmatchev = set_vertex_attr(graph_topmatchev, "color.pan", value=colorpan[as.character(vertex_fampan)])
	graph_topmatchev = set_vertex_attr(graph_topmatchev, "nb.copy.fam", value=as.character(vertex_famcopy))

	save(graph_topmatchev, colorpop, shaperep, colorpan, file=nftopassannotgraph)
	
	if (plotnw){
		pdf(paste(file.path(dirout, prefixout), 'top_associations_network_population_colours.pdf', sep='.'), height=20, width=30)
		plot(graph_topmatchev, vertex.color=vertex_attr(graph_topmatchev, "color.pop"), vertex.shape=vertex_attr(graph_topmatchev, "shape.rep"), vertex.label=NA, vertex.size=2.5, vertex.frame.color=graph_topmatchev$color.pan)
		legend('topleft', legend=levels(vertex_pop), fill=colorpop, ncol=(nlevels(vertex_pop)%/%100)+1)
		dev.off()
	}
}

if (endscript<=5){ quit(save='no') }

nftopasscomm = paste(file.path(dirout, prefixout), 'top_association_annotated_network_communities.RData', sep='.')
if ((!is.null(opt$load_annot_graph) | loadup>=6) & file.exists(nftopasscomm)){
	load(nftopasscomm)
	if (verbose) print(sprintf("loaded community structures of network of top associated gene lineages from file: '%s'", nftopasscomm))
}else{
	### step 6: community analysis: clustering 
	if (verbose) print("compute community structures of network of top associated gene lineages, and plot/write annotation of the components separately")
	comm_topmatchev = mclapply(comm.algos, computeCommunities, graph=graph_topmatchev, mc.cores=ncores) ; names(comm_topmatchev) = comm.algos
	comm.sizes = lapply(comm_topmatchev, function(comm){ sort(sapply(groups(comm), length)) })
	save(comm_topmatchev, comm.sizes, file=nftopasscomm)
}
if (endscript<=6){ quit(save='no') }

### step 7: community analysis: plot 
for (coma in comm.algos){
	dircoma = paste(file.path(dirout, prefixout), sprintf("top_associations_network-%s_communities", coma), sep='.')
	if ((!is.null(opt$load_annot_graph) | loadup>=7) & file.exists(dircoma)){
		if (verbose) print(sprintf("assume table and pdf output files for step 7 (%s clustering) are already in folder: '%s'; SKIP step 7", coma, dircoma))
	}else{
		print(sprintf("plotting %s communities", coma))
		dir.create(dircoma, showWarnings=F)
		
		if (is.null(dbcon)) dbcon = connectionDB(opt)
		
		if (plotnw){
			pdf(paste(file.path(dirout, prefixout), sprintf("top_associations_network-%s_communities_population_colours.pdf", coma), sep='.'), height=12, width=20)
			devbare = dev.cur()
			pdf(paste(file.path(dirout, prefixout), sprintf("top_associations_network-%s_communities_genefam_colours.pdf", coma), sep='.'), height=12, width=20)
			devlabel = dev.cur()
		}
		comm = comm_topmatchev[[coma]]
		for (icomg in 1:length(groups(comm))){
			comg = groups(comm)[[icomg]]
			sg = induced_subgraph(graph_topmatchev, comg)
			if (plotnw){
				if (length(comg)<=plotlabellim){ vlabs = vertex_attr(sg, "name")
				}else{ vlabs = NA }
				dev.set(devbare)
				prespop = unique(sort(vertex_attr(sg, "spe.pop")))
				plot(sg, vertex.color=vertex_attr(sg, "color.pop"), vertex.shape=vertex_attr(sg, "shape.rep"), vertex.label=vlabs, vertex.label.dist=0.2, vertex.size=2.5, vertex.frame.color=sg$color.pan)
				legend('topleft', legend=prespop, fill=colorpop[prespop], ncol=(length(prespop)%/%20)+1)
				dev.set(devlabel)
				presfam = unique(sort(vertex_attr(sg, "fam")))
				colorfam = rainbow(length(presfam)) ; names(colorfam) = presfam
				graph_topmatchev = set_vertex_attr(graph_topmatchev, "color.fam", value=colorfam[as.character(vertex_fam)])
				plot(sg, vertex.color=colorfam[vertex_attr(sg, "fam")], vertex.shape=vertex_attr(sg, "shape.rep"), vertex.label=vlabs, vertex.label.dist=0.2, vertex.size=2.5, vertex.frame.color=sg$color.pan)
				legend('topleft', legend=presfam, fill=colorfam[presfam], ncol=(length(presfam)%/%20)+1)
			}
			# save tables describing communities
			submultiannot = toplinemultiannot[toplinemultiannot$rlocds_id %in% as.numeric(comg),]
			submultidetail = merge(merge(topmatchev, submultiannot, by.x='rlocdsid1', by.y='rlocds_id'), submultiannot, by.x=c('rlocdsid2', repli.invar), by.y=c('rlocds_id', repli.invar), suffixes=1:2)
			write.table(submultidetail, file=paste(file.path(dircoma, sprintf("top_associations_network-%s_community%d_same-replicon_gene_annot.tab", coma, icomg))), sep='\t', row.names=F)
			subreprannot = toplinereprannot[toplinereprannot$rlocds_id %in% as.numeric(comg),]
			subreprdetail = merge(merge(topmatchev, subreprannot, by.x='rlocdsid1', by.y='rlocds_id'), subreprannot, by.x='rlocdsid2', by.y="rlocds_id", suffixes=1:2)
			subreprdetail = subreprdetail[order(subreprdetail$repr_cds_code1, subreprdetail$repr_cds_code2), outannotcols]
			write.table(subreprdetail, file=paste(file.path(dircoma, sprintf("top_associations_network-%s_community%d_repr_gene_annot.tab", coma, icomg))), sep='\t', row.names=F)
			if (length(comg)<=eventstablim){ 
				# get detail of event
				alleventsincomm = NULL
				for (rlocdsid in comg){
					lineageeventsquery = sprintf("
					select gene_family_id, rlocds_id, event_id, freq, event_type, rec.branch_name as rec_branch, don.branch_name as don_branch
					from phylogeny.gene_lineage_events 
					inner join replacement_label_or_cds_code2gene_families using (replacement_label_or_cds_code) 
					inner join phylogeny.species_tree_events using (event_id) inner join phylogeny.species_tree as rec on rec_branch_id=rec.branch_id 
					left join phylogeny.species_tree as don on don_branch_id=don.branch_id 
					where rlocds_id=%s;
					", rlocdsid)
					if (verbose) cat(lineageeventsquery)
					alleventsonlineage = dbGetQuery(dbcon, lineageeventsquery)
					alleventsonlineage$gene_family_id = as.factor(trimws(as.character(alleventsonlineage$gene_family_id)))
					if (is.null(alleventsincomm)){ alleventsincomm = alleventsonlineage
					}else{ alleventsincomm = rbind(alleventsincomm, alleventsonlineage) }
				}
				write.table(alleventsincomm, file=paste(file.path(dircoma, sprintf("top_associations_network-%s_community%d_all_lineage_events.tab", coma, icomg))), sep='\t', row.names=F, quote=F)
			}
		}
		if (plotnw){
			dev.off(devbare)
			dev.off(devlabel)
		}
}
}

if (endscript<=7){ quit(save='no') }

### step 8: plot heatmap of lineages association intensity as projected on maps of selected replicons
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
