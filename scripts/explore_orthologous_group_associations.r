#!/usr/bin/env Rscript
library(gplots)
library(igraph)
library(getopt)
library(parallel)
library(RPostgreSQL)
library(MASS)

options(width=140)

eventstablim = 50
plotlabellim = 100

aggreg.coev.metrics = c('coev_score_mean','coev_score_min', 'coev_score_Q.25', 'coev_score_Q.50', 'coev_score_Q.75', 'coev_score_max')
repli.invar = c("assembly_id", "genomic_accession", "replicon_type", "replicon_size")
comm.algos =  c("louvain", "fast_greedy", "infomap", "label_prop")

outannotcols = c(aggreg.coev.metrics, sapply(1:2, function(i){ 
	paste( c("gene_family_id", "og_id", "repr_lineage", "repr_cds_code", "repr_product", 
	          "tot_og_freq", "chromorplas", "species_pop", "coreoracc", "go_terms", "pathways"), i, sep='_')
}))

plotCoEvolQuantiles = function(qmatchog, detail.quant, sample.name=''){
	layout(matrix(c(1,1,1,2), 1, 4, byrow=T))
	plot(qmatchog[as.character(0:100/100)], xlab='quantiles', ylab='co-evolution score',
	 main=sprintf("distribution of co-evolution score\nin sample %s", sample.name))
	plot(qmatchog[c(detail.quant, 1.0)], xlab='quantiles')
}

parseOGMatchFile = function(nfogmatch, quant.only=F, nmaxrlocsid=NULL, refrepliordlineages=NULL, aggregate.metric='coev_score_mean', 
                                detail.quant=((0:10)*.0001)+.9990, top.quant.cutoff=.99, top.val.cutoff=NULL,
                                comp.quant=T, plot.quant=F, mc.cores=1, verbose=F){
	# load OG match report
	if (verbose) print(sprintf("loading %s ...", nfogmatch), quote=F)
	matchog = read.table(nfogmatch, sep='\t', header=F, stringsAsFactors=F)
	colnames(matchog) = c('gene_family_id_1', 'og_id_1', 'gene_family_id_2', 'og_id_2', 
	                      'n_links', 'n_besthit_links', aggreg.coev.metrics)
	p = unique(sort(c(0:99/100, detail.quant, top.quant.cutoff, 1.0)))
	ncomp = nrow(matchog)
	
	## empirical distrubution of co-evolution scores
	if (comp.quant){
		# compute
		qmatchog = quantile(matchog[, aggregate.metric], probs=p, names=F) ; names(qmatchog) = as.character(p)
		if (plot.quant){
			plotCoEvolQuantiles(qmatchog, detail.quant, sample.name=basename(nfogmatch))
		}
	}else{
		qmatchog = NULL
	}
	qn = c(qmatchog, ncomp) ; names(qn) = c(names(qmatchog), 'nb.reported.comp')
	if (verbose) print(c(sprintf("%s ... done. quantiles:", nfogmatch), qn), quote=F)
	# can stop here
	if (quant.only){ return( qn ) }
	## get set of unique query lineages and can derive the total count of comparisons that could have been attempted
	urlocdsid1 = unique(matchog$rlocdsid1)
	if (!is.null(nmaxrlocsid)){ totcomp = sum(sapply(urlocdsid1, function(k){ nmaxrlocsid - k }))
	}else{ totcomp = NULL }
	## reduce the table to top hits
	if (!is.null(top.val.cutoff)){
		topmatchog = matchog[matchog[, aggregate.metric]>=top.val.cutoff,]
		topstr = sprintf("scores >= %.2g", top.val.cutoff)
	}else{
		topmatchog = matchog[matchog[, aggregate.metric]>=qmatchog[as.character(top.quant.cutoff)],]
		topstr = sprintf("top %.g", 1-top.quant.cutoff)
	}
	if (verbose){ print(sprintf("%s ... done. Top co-evolution scores (%s) [OG-OG aggregate metric: %s]:", nfogmatch, topstr, aggregate.metric), quote=F) ; print(head(topmatchog)) ; print(sprintf("... (%d lines)", nrow(topmatchog))) }
	## can project results on a reference replicon map (the resulting matrix can be summed with those from other partial results)
	if (!is.null(refrepliordlineages)){
		if (verbose){ print("build matrix of gene-to-gene coevolution score based on selected replicon map", quote=F) }
		matmatchog = data.matrix(apply(refrepliordlineages[,c('gene_family_id', 'og_id')], 1, function(x){
			genefamid1 = x[1] ; ogid1 = x[2]
			i1 = which(matchog$gene_family_id_1==genefamid1 & matchog$og_id_1==ogid1)
			i2 = which(matchog$gene_family_id_2==genefamid1 & matchog$og_id_2==ogid1)
			rowmat = apply(refrepliordlineages[,c('gene_family_id', 'og_id')], 1, function(y){
				genefamid2 = y[1] ; ogid2 = y[2]
				j1 = which(matchog$gene_family_id_1[i2]==genefamid2 & matchog$og_id_1[i2]==ogid2)
				j1 = which(matchog$gene_family_id_2[i1]==genefamid2 & matchog$og_id_2[i1]==ogid2)
				k = c(i1[j2], i2[j1])
				if (length(k)==1){ return(matchog[k, aggregate.metric])
				}else{ if (length(k)==0){ return(NA) 
				}else{ print(matchog[k,]) ; stop(sprintf("multiple co-evolution score values for (fam, og) pair %s, %d, %s, %d", genefamid1, ogid1, genefamid2, ogid2)) }}
			})
			if (verbose) cat(rlocdsid1, " ")
			return(rowmat) 
		}))
		if (verbose) cat("... done.\n")
	}else{ matmatchog = NULL }
	return( list(quantiles=qmatchog, unique.query.rlocdsid=urlocdsid1, nb.reported.comp=ncomp, tot.comp=totcomp, top.matches=topmatchog, ref.repli.proj.mat=matmatchog) )
}

plotHMevents = function(mat, matname, cap=99, excl.neighbour=5, fixed.max=30){
	if (fixed.max){ brk = 0:fixed.max
	}else{ brk = quantile(mat, (0:30)/30, na.rm=T) }
	print(brk)
	mat[is.na(mat)] = 0
	heatmap.2(mat, Rowv=F, Colv=F, dendrogram='none', scale='none', trace='none', breaks=brk, col=rainbow(30), na.col='white', main=matname)
	heatmap.2(mat, scale='none', trace='none', breaks=brk, col=rainbow(30), na.col='white', main=matname)
#~ 	for (i in 1:(ncol(mat)-1)){
#~ 	 for (j in 2:ncol(mat)){
#~ 	  mat[j,i] = mat[i,j]
#~ 	}}
	layout(matrix(1:2, 2,1))
	sme = apply(mat, 1, sum, na.rm=T)
	barplot(sme, las=2, main=sprintf("%s\nsummed matched event freq.", matname))
	barplot(sapply(1:nrow(mat), function(i){
		js = sapply(seq(i-excl.neighbour, i+excl.neighbour), function(k){ ifelse(k>0, ifelse(k<=ncol(mat), k, k-nrow(mat)), k+ncol(mat)) })
		smen = sum(mat[i,-js], na.rm=T) 
		return(smen)
	}), las=2, main=sprintf("%s\nsummed matched event freq. excluding %d next genes", matname, excl.neighbour))
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

#~ dirmatchogents = '/enterobac/PanteroDB_v0.3/07.compare_scenarios/bs_stem70_withinmedian35/replaceCCinGasinS-collapsePOPinSnotinG/ale_collapsed_undat/gene_event_matches'
#~ dirmatchogents = '/enterobac/PanteroDB_v0.3/07.compare_scenarios/bs_stem70_withinmedian35/replaceCCinGasinS-collapsePOPinSnotinG/ale_collapsed_undat/gene_event_matches_10253genes_06062018/plasmid2genefamily.mat_excl-transposases.NC_022651.1-centred_community_size13'
#~ refrepli = 'NC_022651.1'
#~ nfrefrepliordlineages = '/enterobac/PanteroDB_v0.3/07.compare_scenarios/bs_stem70_withinmedian35/replaceCCinGasinS-collapsePOPinSnotinG/ale_collapsed_undat/gene_event_matches_10253genes_06062018/plasmid2genefamily.mat_excl-transposases.NZ_CP008933.1-centred_community_size24/NC_022651.1_ordered_gene_lineage_products.tab'
#~ nftopogreprannot = '/enterobac/PanteroDB_v0.3/07.compare_scenarios/bs_stem70_withinmedian35/replaceCCinGasinS-collapsePOPinSnotinG/ale_collapsed_undat/gene_event_matches_10253genes_06062018/covered_10253genes_lineages.details'
#~ nmaxrlocsid = 640228
#~ opt = list()

### argument (option) parsing
spec = matrix(c(
  'input_match_ogs',    'i', 1, "character", "[mandatory] (preferably absolute) path of file listing associated ortholog groups (OGs) (output of condense_coevolution_network_by_orthologous_groups.py)",
  'ref_replicons',      'r', 2, "character", "comma-separated list of replicons on which to project the LD-like matrix of co-evolution scores",
  'replicon_annot_map', 'm', 2, "character", "path to table file of gene lineage ordered as in reference replicons, with functional annotation of gene products; table must be tab-delimited and contain at least columns named: 'gene_family_id', 'og_id', 'cds_code', 'product'",
  'nb_genomes',         'G', 1, "integer",   "number of genomes on the dataset",
  'output_dir',         'o', 2, "character", "path to folder where to write output files; default to the parent folder of input match event folder",
  'output_prefix',      'p', 2, "character", "prefix for output files; default to the name of input match event folder",
  'aggregate_metric',   'g', 2, "character", paste0("name of aggregate metric used to summarize coevolution scores between single gene lineages of two OGs; possible values are: '", paste0(aggreg.coev.metrics, collapse="', '", "'")),
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
  'load_quantiles',     'Q', 0, "logical",   "will attempt to load quantile information from R object archive file saved in a previous run; file must be named like `dir_match_events` argument value with prefix '.co-evolution_quantiles.RData', and contain a data.frame named `matchogdf` (bypass quantile computation)",
  'load_top_ass',       'A', 0, "logical",   "will attempt to load top lineage association table from R object archive file saved in a previous run; file must be named like `dir_match_events` argument value with prefix '.top_association.RData', and contain a list named `lparsematchog` (bypass quantile computation and top association retrieval)",
  'load_top_network',   'N', 0, "logical",   "will attempt to load top lineage association network from R object archive file saved in a previous run; file must be named like `dir_match_events` argument value with prefix '.top_association_network.RData', and contain a list named `lparsematchog` (bypass quantile computation, top association retrieval and graph computation)",
  'load_up_to_step',    'L', 2, "integer",   "will attempt to load data from R object archive file saved during a previous run, up to the desired step (1:quantile computation; 2:top association retrieval; 3:top association graph computation; 4:top associated lineage annotation retrival from SQL db)",
  'plot_network',       'n', 0, "logical",   "enable plotting of the networks (long)",
  'same_replicon',      'R', 2, "integer",   "enable making and writing annotation table for associated gene pairs from the same replicon; integer argument provide max limit rows per associated OG pair",
  'stop_after_step',    'Z', 2, "integer",   "will run until reaching the desired progress of computation, then quit (steps: quantile computation=1; top association retrieval (+OPTIONALLY: plot heatmaps of lineages association intensity as projected on maps of selected replicons)=2; top association graph computation=3; top associated lineage annotation retrival from SQL db=4; compute association network=5; search communities in network=6; plot top association network with metadata=7 [end])"
), byrow=TRUE, ncol=5)
opt = getopt(spec, opt=commandArgs(trailingOnly=T))

print("options:", quote=F)
print(opt)

if (!is.null(opt$verbose)){ verbose = TRUE }else{ verbose = FALSE }
if (!is.null(opt$plot_network)){ plotnw = TRUE }else{ plotnw = FALSE }

if ( is.null(opt$input_match_ogs) ){ 
	stop("path to input file required! please specify it through '--input_match_ogs' option")  
}else{ 
	nfmatchog = opt$input_match_ogs
	print(sprintf("will load lineage co-evolution scores from file '%s'", nfmatchog), quote=F)
}
if ( is.null(opt$output_dir) ){ 
	dirout = dirname(nfmatchog)  
}else{
	dirout = opt$output_dir
}
if ( is.null(opt$output_prefix) ){ 
	prefixout = basename(nfmatchog)  
}else{
	prefixout = opt$output_prefix
}
qcutoff = scutoff = NULL
if (is.null(opt$score_cutoff)){
	if (is.null(opt$quant_cutoff)){
		qcutoff = .99
	}else{
		qcutoff = opt$quant_cutoff
	}
}else{
	scutoff = opt$score_cutoff
}
if (is.null(opt$aggregate_metric)){
	aggregate.metric = 'coev_score_mean'
}else{
	aggregate.metric = opt$aggregate_metric
}
refrepli = opt$ref_replicons
nmaxrlocsid = opt$max_lineage_id
nfrefrepliordlineages = opt$replicon_annot_map
if (!is.null(opt$replicon_annot_map)){
	refrepliordlineages = read.table(opt$replicon_annot_map, h=T, sep='\t', quote='')
	if ( is.null(refrepli) ){ 
		# try some matching in file name...
		refrepli = gsub("map_(.+)_cds.tab", "\\1", basename(opt$replicon_annot_map))
	}
}else{
	refrepliordlineages = NULL
}
if (is.null(opt$same_replicon)){
	samereplimax = NULL
}else{
	samereplimax = opt$same_replicon
}
if (is.null(opt$threads)){
	ncores = 6
}else{
	ncores = opt$threads
}
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
if (!is.null(opt$nb_genomes)){
	nb_genome = as.numeric(opt$nb_genomes)
}else{
	stop("the number of genomes is a required parameter; specify it with option --nb_genomes")
}

### step 1: parse data to extract co-evolution score quantiles
if (!is.null(qcutoff)){
	nfsavematchog = paste(file.path(dirout, prefixout), "co-evolution_quantiles.RData", sep='.')
	if ((!is.null(opt$load_quantiles) | loadup>=1) & file.exists(nfsavematchog)){
		load(nfsavematchog)
		if (verbose) print(sprintf("loaded quantiles of co-evolution score from pre-exesting file: '%s", nfsavematchog), quote=F)
	}else{
		if (verbose) print("will read OG co-evolution score tables a first time to compute quantiles", quote=F)
		qmatchog = parseOGMatchFile(nfmatchog, quant.only=T, top.quant.cutoff=qcutoff, verbose=verbose, mc.cores=ncores, aggregate.metric=aggregate.metric)

		save(qmatchog, file=nfsavematchog)
		if (verbose) print(sprintf("saved quantiles of OG co-evolution score to file: '%s", nfsavematchog), quote=F)
	}
	q1 = as.character(0:100/100)
	q2 = names(qmatchog)[!(names(qmatchog) %in% c(q1, 'nb.reported.comp'))]
	pdf(paste(file.path(dirout, prefixout), "co-evolution_score_quantiles.pdf", sep='.'), height=6, width=10)
	layout(matrix(c(1,1,1,2), 1, 4, byrow=T))
	plot(y=qmatchog[q1], x=as.numeric(q1), xlab="quantiles", ylab=sprintf("aggregate co-evolution score (%s)", aggregate.metric), main="distribution of OG-OG co-evolution scores")
	plot(y=qmatchog[q2], x=as.numeric(q2), xlab="quantiles", ylab=sprintf("aggregate co-evolution score (%s)", aggregate.metric))
	dev.off()
	
	totrepcomp = qmatchog['nb.reported.comp']
	scutoff = qmatchog[as.character(qcutoff)]

	print(sprintf("found a value of %f for quantile at p=%f of empirical distribution of co-evolution scores between %g reported gene (gene family, OG) comparisons",
				   scutoff, qcutoff, totrepcomp), quote=F)
}
if (endscript<=1){ quit(save='no') }

### step 2: parse data to extract top co-evolution scores
nftopass = paste(file.path(dirout, prefixout), "top_associations.RData", sep='.')
nftopasstab = paste(file.path(dirout, prefixout), "top_associations.tab", sep='.')
if ((!is.null(opt$load_top_ass) | loadup>=2) & file.exists(nftopasstab) & is.null(refrepliordlineages)){
	topmatchog = read.table(nftopasstab, sep='\t', header=T)
	if (verbose) print(sprintf("loaded top association data from file: '%s'", nftopasstab))
}else{
	if ((!is.null(opt$load_top_ass) | loadup>=2) & file.exists(nftopass)){
	load(nftopass)
	if (verbose) print(sprintf("loaded top association data from file: '%s'", nftopass))
	}else{
		if (verbose) print("extract top associations from bulk data")
		parsedmatchog = parseOGMatchFile(nfmatchog, comp.quant=F, top.val.cutoff=scutoff, aggregate.metric=aggregate.metric,
								 refrepliordlineages=refrepliordlineages, verbose=verbose, mc.cores=ncores)
		save(parsedmatchog, file=nftopass)
	}
	topmatchog = parsedmatchog$top.matches
	### OPTIONAL: plot heatmap of lineages association intensity as projected on maps of selected replicons
	if (!is.null(refrepliordlineages)){
		nfmatrepli = paste(file.path(dirout, prefixout), sprintf("co-evol_scores_projection_%s.mat.RData", refrepli), sep='.')
		if ((!is.null(opt$load_top_ass) | loadup>=2) & file.exists(nfmatrepli)){
			load(nfmatrepli)
		}else{
			matmatchog = lparsematchog[[1]]$ref.repli.proj.mat
			if (length(lparsematchog)>1){
				for (i in 2:length(lparsematchog)){
					matmatchog = matmatchog + lparsematchog[[i]]$ref.repli.proj.mat
				}
			}
			colnames(matmatchog) = rownames(matmatchog) = refrepliordlineages$cds_code
			save(matmatchog, file=nfmatrepli)
		}
		nfpdf = paste(file.path(dirout, prefixout), sprintf("co-evol_scores_projection_%s.pdf", refrepli), sep='.')
		pdf(nfpdf, height=20, width=20)
		plotHMevents(matmatchog, sprintf("co-evolution scores projected on %s map", refrepli))
		dev.off()
		print(nfpdf)
	}
	write.table(topmatchog, file=nftopasstab, sep='\t', row.names=F)
}
if (endscript<=2){ quit(save='no') }

### step 3: create network of gene lineage association
nftopassgraph = paste(file.path(dirout, prefixout), "top_association_network.RData", sep='.')
if ((!is.null(opt$load_top_network) | loadup>=3) & file.exists(nftopassgraph)){
	load(nftopassgraph)
	if (verbose) print(sprintf("loaded network of top associated gene lineages from file: '%s'", nftopassgraph))
}else{
	if (verbose) print("compute network of top associated gene lineages")
	topmatchfamog = as.data.frame(cbind(
	 paste(topmatchog[,'gene_family_id_1'], topmatchog[,'og_id_1'], sep='_'),
	 paste(topmatchog[,'gene_family_id_2'], topmatchog[,'og_id_2'], sep='_')
	), stringsAsFactors=F)
	topmatchfamog[[aggregate.metric]] = topmatchog[[aggregate.metric]]
	colnames(topmatchfamog) = c('genefamog_1', 'genefamog_2', aggregate.metric)
	graph_topmatchog = graph_from_data_frame(topmatchfamog, directed=F)
	save(graph_topmatchog, topmatchfamog, file=nftopassgraph)
}
#~ pdf(file=paste(file.path(dirout, prefixout), 'top_association_network.pdf', sep='.'), height=30, width=30)
#~ plot(graph_topmatchog)
#~ dev.off()
if (endscript<=3){ quit(save='no') }

### step 4: collect detailed annotation data on top associated lineages
if (!is.null(opt$db_name)){
	if (verbose) print(sprintf("retrieve annotation of (genes in) top associated OGs from %s database: '%s'", dbtype, opt$db_name))
	
	if (is.null(dbcon)) dbcon = connectionDB(opt)
	
	u1 = unique(topmatchfamog$genefamog_1)
	u2 = unique(topmatchfamog$genefamog_2)
	u = as.data.frame(t(simplify2array(strsplit(union(u1, u2), split='_'))), stringsAsFactors=F)
	colnames(u) = c('gene_family_id', 'og_id')
	u$og_id = as.numeric(u$og_id)
	if (verbose){ print("toporthogroups:"); print(head(u)) }
	
	dbExecute(dbcon, "CREATE TEMP TABLE toporthogroups (gene_family_id varchar(20), og_id smallint); ")
	dbWriteTable(dbcon, "toporthogroups", value=u, append=T, row.names=F)
	
	nftopmatchogdetail = paste(file.path(dirout, prefixout), "top_associated_OG_repr_annot.RData", sep='.')
	if ((!is.null(opt$load_top_annot) | loadup>=4) & file.exists(nftopmatchogdetail)){
		load(nftopmatchogdetail)
	if (verbose) print(sprintf("loaded annotation of (genes in) top associated OGs from file: '%s'", nftopmatchogdetail))
	}else{
		## query annotations - 1 row / OG, with a representative lineage / CDS
		
		# first generate table of enriched gene family_OG annotation
		createfamogstat = "
		  CREATE TEMP TABLE toporthogroupstats AS
		   SELECT gene_family_id, og_id, max(cds_code) AS cds_code,
		    count(cds_code) AS tot_og_freq, count(distinct assembly_id) AS og_pres_freq, 
		    sum(CASE WHEN replicon_type='chromosome' THEN 1 ELSE 0 END) AS count_chromosome, 
		    sum(CASE WHEN replicon_type='plasmid' THEN 1 ELSE 0 END) AS count_plasmid
		     FROM toporthogroups
		     INNER JOIN orthologous_groups using (gene_family_id, og_id)
		     INNER JOIN gene_tree_label2cds_code USING (replacement_label_or_cds_code)
		     INNER JOIN coding_sequences USING (gene_family_id, cds_code)
		     INNER JOIN replicons USING (genomic_accession)
		   GROUP BY gene_family_id, og_id ;
		   "
		if (verbose) cat(createfamogstat)
		dbExecute(dbcon, createfamogstat)
		
		### !! replace  (genbank_nr_protein_id) with (nr_protein_id)
		annotreprquery = "
		SELECT DISTINCT gene_family_id, og_id,
		 replacement_label_or_cds_code AS repr_lineage, cds_code AS repr_cds_code, product as repr_product,
		 tot_og_freq, og_pres_freq, count_chromosome, count_plasmid, go_terms, pathways
		  FROM toporthogroupstats
		  INNER JOIN coding_sequences USING (gene_family_id, cds_code)
		  INNER JOIN gene_tree_label2cds_code USING (cds_code)
		  INNER JOIN proteins USING (genbank_nr_protein_id)
		  LEFT JOIN prot2goterms_oneliner USING (genbank_nr_protein_id)
		  LEFT JOIN prot2pathways_oneliner USING (genbank_nr_protein_id)
		ORDER BY repr_cds_code ;
		"
		if (verbose) cat(annotreprquery)
		topogreprannot = dbGetQuery(dbcon, annotreprquery)
		topogreprannot$gene_family_id = as.factor(trimws(as.character(topogreprannot$gene_family_id)))
		topogreprannot$chromorplas = as.factor(ifelse(topogreprannot$count_chromosome>0, ifelse(topogreprannot$count_plasmid==0, 'chromosome', 'both'), ifelse(topogreprannot$count_plasmid>0, 'plasmid', 'none')))
		topogreprannot$species_pop = as.factor(sapply(strsplit(as.character(topogreprannot[,'repr_lineage']), split='_'), `[`, 1))
		topogreprannot$coreoracc = ifelse(topogreprannot$tot_og_freq >= nb_genome*0.98, "core", "accessory")
		write.table(topogreprannot, file=paste(file.path(dirout, prefixout), "top_OG_repr_annot.tab", sep='.'), sep='\t', row.names=F)

		topmatchogdetail = merge(merge(topmatchog, topogreprannot, by.x=c('gene_family_id_1','og_id_1'), by.y=c('gene_family_id','og_id')),
		                         topogreprannot, by.x=c('gene_family_id_2','og_id_2'), by.y=c('gene_family_id','og_id'), suffixes=paste0('_', 1:2))
		if (verbose){ 
			print("topmatchogdetail (dim, head):")
			print(dim(topmatchogdetail))
			print(head(topmatchogdetail))
		}
		topmatchogdetail = topmatchogdetail[order(topmatchogdetail[['repr_cds_code_1']], topmatchogdetail[['repr_cds_code_2']]), outannotcols]
		write.table(topmatchogdetail, file=paste(file.path(dirout, prefixout), "top_associated_OG_repr_annot.tab", sep='.'), sep='\t', row.names=F)
		save(topogreprannot, topmatchogdetail, file=nftopmatchogdetail)
	}
	
	if (!is.null(samereplimax)){
		nftopmatchogmultidetail = paste(file.path(dirout, prefixout), "top_associated_OG_same-replicon_gene_annot.RData", sep='.')
		if ((!is.null(opt$load_top_annot) | loadup>=4) & file.exists(nftopmatchogmultidetail)){
			load(nftopmatchogmultidetail)
		if (verbose) print(sprintf("loaded annotation of gene pairs from same replicons in top associated OGs from file: '%s'", nftopmatchogmultidetail))
		}else{
			if (verbose) print(sprintf("retrieve annotation of gene pairs from same replicons in top associated OGs from file: '%s'", nftopmatchogmultidetail))
			## query annotations - for all CDS in the lineages in th OGs
			multiannotquery = "
			SELECT DISTINCT gene_family_id, og_id, replacement_label_or_cds_code, cds_code, product,
			  assembly_id, genomic_accession, replicon_type, replicon_size, cds_begin, go_terms, pathways
			  FROM toporthogroups
			  INNER JOIN orthologous_groups USING (gene_family_id, og_id)
			  INNER JOIN gene_tree_label2cds_code USING (replacement_label_or_cds_code)
			  INNER JOIN coding_sequences USING (gene_family_id, cds_code)
			  INNER JOIN proteins USING (genbank_nr_protein_id)
			  INNER JOIN replicons USING (genomic_accession)
			  LEFT JOIN prot2goterms_oneliner USING (genbank_nr_protein_id)
			  LEFT JOIN prot2pathways_oneliner USING (genbank_nr_protein_id);"
			if (verbose) cat(gsub("^\t\t", "", multiannotquery))
			toplinemultiannot = dbGetQuery(dbcon, multiannotquery)
			dbExecute(dbcon, "drop table toporthogroups; ")
			toplinemultiannot$assembly_id = as.factor(trimws(as.character(toplinemultiannot$assembly_id)))
			toplinemultiannot$replicon_type = as.factor(trimws(as.character(toplinemultiannot$replicon_type)))
			toplinemultiannot$gene_family_id = as.factor(trimws(as.character(toplinemultiannot$gene_family_id)))
			if (verbose){ 
				print("toplinemultiannot (dim, head):")
				print(dim(toplinemultiannot))
				print(head(toplinemultiannot))
			}
			# add the details of all genes represented in the lineages
			# but keep only rows where the associated genes are on the same replicon.
			# very heavy merge operation
			gc()
#~ 			topmatchogmultidetail = merge(merge(topmatchog, toplinemultiannot, by.x=c('gene_family_id_1','og_id_1'), by.y=c('gene_family_id','og_id')),
#~ 										  toplinemultiannot, by.x=c('gene_family_id_2','og_id_2', repli.invar), by.y=c('gene_family_id','og_id', repli.invar), suffixes=1:2)
			# try to proceed with divide-and-conquer approach instead, proceeding by replicon
			replicons = unique(toplinemultiannot[['genomic_accession']])
			ltopmatchogmultidetail = mclapply(replicons, function(repli){
				replitoplinemultiannot = toplinemultiannot[toplinemultiannot[['genomic_accession']]==repli,]
				replitopmatchogmultidetail = merge(merge(topmatchog, replitoplinemultiannot, by.x=c('gene_family_id_1','og_id_1'), by.y=c('gene_family_id','og_id')),
			                              toplinemultiannot, by.x=c('gene_family_id_2','og_id_2', repli.invar), by.y=c('gene_family_id','og_id', repli.invar), suffixes=1:2)
				if (verbose) print(repli)
				if (verbose) print(head(replitopmatchogmultidetail, n=2))
				if (nrow(replitopmatchogmultidetail) == 0){
					return(NA)
				}else{ if (nrow(replitopmatchogmultidetail) <= samereplimax){
					return(replitopmatchogmultidetail)
				}else{ 
					return(replitopmatchogmultidetail[sample.int(nrow(replitopmatchogmultidetail), size=samereplimax),])
				}}
			}, mc.cores=ncores)
			ltopmatchogmultidetail = ltopmatchogmultidetail[which(!is.na(ltopmatchogmultidetail))]
			if (length(ltopmatchogmultidetail) > 1){
				topmatchogmultidetail = do.call(rbind, ltopmatchogmultidetail)
			}else{ topmatchogmultidetail = ltopmatchogmultidetail[[1]] }
			rm(ltopmatchogmultidetail)
			gc()
			if (verbose){ 
				print("topmatchogmultidetail (dim, head):")
				print(dim(topmatchogmultidetail))
				print(head(topmatchogmultidetail))
			}
			topmatchogmultidetail = topmatchogmultidetail[order(topmatchogmultidetail[['cds_code1']], topmatchogmultidetail[['cds_code2']]), c(repli.invar, sort(setdiff(colnames(topmatchogmultidetail), repli.invar)))]
			topmatchogmultidetail$physical_dist = apply(topmatchogmultidetail[,c('cds_begin1', 'cds_begin2', 'replicon_size')], 1, function(x){
				if (x[1] < x[2]){ a = x[1] ; b = x[2] }else{ b = x[1] ; a = x[2] }
				return(min((b - a), (x[3] - b + a)))
			})
			write.table(topmatchogmultidetail, file=paste(file.path(dirout, prefixout), "top_associated_OG_same-replicon_gene_annot.tab", sep='.'), sep='\t', row.names=F)
			
			save(toplinemultiannot, topmatchogmultidetail, lmgenedistfreq, file=nftopmatchogmultidetail)
			
			# plot relation between co-evolution scaore of top associated lineages and the physical or locus distance of their representatives on replicons
			pdf(paste(file.path(dirout, prefixout), "top_associated_OG_same-replicon.jointfreq_vs_physical_dist.pdf", sep='.'), height=10, width=15)
			contour(kde2d(topmatchogmultidetail[['physical_dist']], topmatchogmultidetail[[aggregate.metric]]), xlab='physical distance', ylab=sprintf("aggregate co-evolution score (%s)", aggregate.metric), main='Co-evolution score and inter-gene distance\n(lineages leading to genes residing in the same genome)')
			lmgenedistfreq = lm(topmatchogmultidetail[[aggregate.metric]] ~ topmatchogmultidetail[['physical_dist']])
			abline(lmgenedistfreq, col='red')
			print(summary(lmgenedistfreq))
			dev.off()
		}
	}
}
if (endscript<=4){ quit(save='no') }

nftopassannotgraph = paste(file.path(dirout, prefixout), "top_association_annotated_network.RData", sep='.')
if ((!is.null(opt$load_annot_graph) | loadup>=5) & file.exists(nftopassannotgraph)){
	load(nftopassannotgraph)
	if (verbose) print(sprintf("loaded annotated version of network of top associated gene lineages from file: '%s'", nftopassannotgraph))
}else{
	### step 5: compute and plot network of top associated lineages along with associated metadata
	if (verbose) print("annotate network of top associated gene lineages from file and plot it in full")
	vertex_pop = as.factor(sapply(names(V(graph_topmatchog)), function(genefamog){
		famog = strsplit(genefamog, split='_')[[1]]
		topogreprannot[topogreprannot$gene_family_id==famog[1] & topogreprannot$og_id==famog[2],'species_pop']
	}))
	vertex_fam = as.factor(sapply(names(V(graph_topmatchog)), function(genefamog){
		famog = strsplit(genefamog, split='_')[[1]]
		topogreprannot[topogreprannot$gene_family_id==famog[1] & topogreprannot$og_id==famog[2],'gene_family_id']
	}))
	vertex_rep = as.factor(sapply(names(V(graph_topmatchog)), function(genefamog){
		famog = strsplit(genefamog, split='_')[[1]]
		topogreprannot[topogreprannot$gene_family_id==famog[1] & topogreprannot$og_id==famog[2],'chromorplas']
	}))
	vertex_fampan = as.factor(sapply(names(V(graph_topmatchog)), function(genefamog){
		famog = strsplit(genefamog, split='_')[[1]]
		topogreprannot[topogreprannot$gene_family_id==famog[1] & topogreprannot$og_id==famog[2],'coreoracc']
	}))
	colorpop = rainbow(nlevels(vertex_pop)) ; names(colorpop) = levels(vertex_pop)
	shaperep = c("circle", "square", "csquare") ; names(shaperep) = c("chromosome", "plasmid", "both")
	colorrep = c("yellow", "red", "orange") ; names(shaperep) = c("chromosome", "plasmid", "both")
	colorpan = c("black", "red") ; names(colorpan) = c("core", "accessory")
	graph_topmatchog = set_vertex_attr(graph_topmatchog, "spe.pop", value=as.character(vertex_pop))
	graph_topmatchog = set_vertex_attr(graph_topmatchog, "color.pop", value=colorpop[as.character(vertex_pop)])
	graph_topmatchog = set_vertex_attr(graph_topmatchog, "repli.type", value=as.character(vertex_rep))
	graph_topmatchog = set_vertex_attr(graph_topmatchog, "shape.rep", value=shaperep[as.character(vertex_rep)])
	graph_topmatchog = set_vertex_attr(graph_topmatchog, "color.rep", value=colorrep[as.character(vertex_rep)])
	graph_topmatchog = set_vertex_attr(graph_topmatchog, "fam", value=as.character(vertex_fam))
	graph_topmatchog = set_vertex_attr(graph_topmatchog, "pang.fam.status", value=as.character(vertex_fampan))
	graph_topmatchog = set_vertex_attr(graph_topmatchog, "color.pan", value=colorpan[as.character(vertex_fampan)])

	save(graph_topmatchog, colorpop, shaperep, colorpan, file=nftopassannotgraph)
	
	if (plotnw){
		pdf(paste(file.path(dirout, prefixout), "top_associations_network_population_colours.pdf", sep='.'), height=20, width=30)
		plot(graph_topmatchog, vertex.color=vertex_attr(graph_topmatchog, "color.pop"), vertex.shape=vertex_attr(graph_topmatchog, "shape.rep"), vertex.label=NA, vertex.size=2.5, vertex.frame.color=graph_topmatchog$color.pan)
		legend('topleft', legend=levels(vertex_pop), fill=colorpop, ncol=(nlevels(vertex_pop)%/%100)+1)
		dev.off()
		pdf(paste(file.path(dirout, prefixout), "top_associations_network_coreoracc_colours.pdf", sep='.'), height=20, width=30)
		plot(graph_topmatchog, vertex.color=vertex_attr(graph_topmatchog, "color.pan"), vertex.shape=vertex_attr(graph_topmatchog, "shape.rep"), vertex.label=NA, vertex.size=2.5)
		legend('topleft', legend=levels(vertex_fampan), fill=colorpan, ncol=(nlevels(vertex_fampan)%/%100)+1)
		dev.off()
		pdf(paste(file.path(dirout, prefixout), "top_associations_network_chromorplas_colours.pdf", sep='.'), height=20, width=30)
		plot(graph_topmatchog, vertex.color=vertex_attr(graph_topmatchog, "color.rep"), vertex.shape=vertex_attr(graph_topmatchog, "shape.rep"), vertex.label=NA, vertex.size=2.5)
		legend('topleft', legend=levels(vertex_rep), fill=colorrep, ncol=(nlevels(vertex_rep)%/%100)+1)
		dev.off()
	}
}

if (endscript<=5){ quit(save='no') }

nftopasscomm = paste(file.path(dirout, prefixout), "top_association_annotated_network_communities.RData", sep='.')
if ((!is.null(opt$load_annot_graph) | loadup>=6) & file.exists(nftopasscomm)){
	load(nftopasscomm)
	if (verbose) print(sprintf("loaded community structures of network of top associated gene lineages from file: '%s'", nftopasscomm))
}else{
	### step 6: community analysis: clustering 
	if (verbose) print("compute community structures of network of top associated gene lineages, and plot/write annotation of the components separately")
	comm_topmatchog = mclapply(comm.algos, computeCommunities, graph=graph_topmatchog, mc.cores=ncores) ; names(comm_topmatchog) = comm.algos
	comm.sizes = lapply(comm_topmatchog, function(comm){ sort(sapply(groups(comm), length)) })
	save(comm_topmatchog, comm.sizes, file=nftopasscomm)
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
		comm = comm_topmatchog[[coma]]
		for (icomg in 1:length(groups(comm))){
			comg = groups(comm)[[icomg]]
			if (verbose) print(sprintf("community %d: vertices %s ...", icomg, paste(head(comg), collapse=', ')))
			sg = induced_subgraph(graph_topmatchog, comg)
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
				graph_topmatchog = set_vertex_attr(graph_topmatchog, "color.fam", value=colorfam[as.character(vertex_fam)])
				plot(sg, vertex.color=colorfam[vertex_attr(sg, "fam")], vertex.shape=vertex_attr(sg, "shape.rep"), vertex.label=vlabs, vertex.label.dist=0.2, vertex.size=2.5, vertex.frame.color=sg$color.pan)
				legend('topleft', legend=presfam, fill=colorfam[presfam], ncol=(length(presfam)%/%20)+1)
			}
			# save tables describing communities
			comannotrows = paste(as.character(topogreprannot[,'gene_family_id',]), as.character(topogreprannot[,'og_id']), sep='_')
#~ 			if (verbose) print(head(comannotrows))
			comannotrowsi = which(comannotrows %in% comg)
			submultiannot = toplinemultiannot[comannotrowsi,]
			submultidetail = merge(merge(topmatchog, submultiannot, by.x=c('gene_family_id_1','og_id_1'), by.y=c('gene_family_id','og_id')), submultiannot, by.x=c('gene_family_id_2','og_id_2', repli.invar), by.y=c('gene_family_id','og_id', repli.invar), suffixes=paste0('_', 1:2))
			write.table(submultidetail, file=paste(file.path(dircoma, sprintf("top_associations_network-%s_community%d_same-replicon_gene_annot.tab", coma, icomg))), sep='\t', row.names=F)
			subreprannot = topogreprannot[comannotrowsi,]
			subreprdetail = merge(merge(topmatchog, subreprannot, by.x=c('gene_family_id_1','og_id_1'), by.y=c('gene_family_id','og_id')), subreprannot, by.x=c('gene_family_id_2','og_id_2'), by.y=c('gene_family_id','og_id'), suffixes=paste0('_', 1:2))
			subreprdetail = subreprdetail[order(subreprdetail[['repr_cds_code_1']], subreprdetail[['repr_cds_code_2']]), outannotcols]
			write.table(subreprdetail, file=paste(file.path(dircoma, sprintf("top_associations_network-%s_community%d_repr_gene_annot.tab", coma, icomg))), sep='\t', row.names=F)
		}
		if (plotnw){
			dev.off(devbare)
			dev.off(devlabel)
		}
}
}

if (endscript<=7){ quit(save='no') }

print("warnings:")
print(warnings())
