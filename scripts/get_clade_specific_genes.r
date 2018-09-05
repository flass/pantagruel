#!/usr/bin/Rscript --vanilla
library('RSQLite')

cargs = commandArgs(trailing=T)
# env variables relating to the full Pantagruel pangenome db, expected to be set
#~ nflasscode =  file.path(Sys.getenv('database'), 'genome_codes.tab')
#~ nffamgenomemat = file.path(Sys.getenv('protali'), 'full_families_genome_counts-noORFans.mat')
#~ sqldb = Sys.getenv('sqldb')
nffamgenomemat = cargs[1]
sqldb = cargs[2]
ogcolid = as.numeric(cargs[3])
# genome subset-specific input files
filerad = cargs[4]
foutrad = cargs[5]
if (length(cargs)>5){
	nflasscode = cargs[6]
}else{
	nflasscode = NULL
}
if (length(cargs)>6){
	preferredgenomes = strsplit(cargs[7], split=',')[[1]]
}else{
	preferredgenomes = c()
}
nfrestrictlist = sprintf("%s_genome_codes", filerad)
nfcladedef = sprintf("%s_clade_defs", filerad)
# output files
nfabspresmat = sprintf("%s_gene_abspres.mat.RData", foutrad)
nfoutspege = sprintf("%s_specific_genes.tab", foutrad)

#~ # clade definitions
#~ if (!file.exists(nfcladedef)){
#~ 	clade0 = c('RHIZOB27', 'RNT25')
#~ 	outgroup0 = 'RHISP2'
#~ 	clade1 = c(clade0, outgroup0)
#~ 	clade2 = c('RKHAN', 'RHAB21', 'RFYW14')
#~ 	clade12 = c(clade1, clade2)
#~ 	clade3 = c('REQ54', 'REJC140')
#~ 	clade123 = c(clade12, clade3)
#~ 	clade4 = c('PSEPEL1', 'PSEPEL2', 'RHIMAR')
#~ 	clade1234 = c(clade123, clade4)
#~ 	outgroup = 'RHIZOB54'
#~ 	clades = list(clade0=clade0, outgroup0=outgroup0, clade1=clade1, clade2=clade2, clade12=clade12, clade3=clade3, clade123=clade123, clade4=clade4, clade1234=clade1234, outgroup=outgroup)
#~ 	sisterclades = list(clade0=outgroup0, outgroup0=clade0, clade1=clade2, clade2=clade1, clade12=clade3, clade3=clade12, clade123=clade4, clade4=clade123, clade1234=outgroup, outgroup=clade1234)
#~ 	write.table(t(sapply(names(clades), function(cla){ sapply(list(clades, sisterclades), function(x){ paste(x[[cla]], collapse=',') }) })),
#~ 	 sep='\t', quote=F, row.names=names(clades), col.names=c('clade', 'sisterclade'),
#~ 	 file=nfcladedef)
#~ }
cladedefcsv = read.table(nfcladedef, sep='\t', header=T, row.names=1, stringsAsFactors=F)
cladedefs = apply(cladedefcsv, 1, strsplit, split=',')

# load gene presence / absence data
if (file.exists(nfabspresmat)){
	load(nfabspresmat)
}else{
	genocount = data.matrix(read.table(file=nffamgenomemat))
	restrictgenomelist = readLines(nfrestrictlist)
	if (!is.null(nflasscode)){
		lasscode = read.table(nflasscode, row.names=1, stringsAsFactors=F)
		colnames(genocount) = lasscode[colnames(genocount),1]
	}
#~ 	print(setdiff(restrictgenomelist, colnames(genocount)))
	genocount = genocount[,restrictgenomelist]
	gc()
	save(genocount, file=nfabspresmat)
}

# compute gene sets
specifigenes = lapply(cladedefs, function(cladedef){
	allpres = apply(genocount[,cladedef$clade,drop=F], 1, function(x){ all(x>0) })
	allabssis = apply(genocount[,cladedef$sisterclade,drop=F], 1, function(x){ all(x==0) })
	return(which(allpres & allabssis))
})

dbcon = dbConnect(SQLite(), sqldb)

for (i in 1:length(cladedefs)){
	cla = names(cladedefs)[i]
	if (!is.null(cladedefcsv$name)){ claname = cladedefcsv$name[[i]] }else{ claname = "" }
	a = as.logical(i-1)
	write(paste(sprintf("# %s %s", cla, claname), paste(c("# gene families present in all genomes of clade:", "and absent in all genomes of sister clade:"), cladedefcsv[cla,c('clade', 'sisterclade')], sep=' ', collapse='; '), sep='\t\t\t\t\t\t'),
	 file=nfoutspege, append=a)
	dbBegin(dbcon)
	spefamogs = as.data.frame(t(sapply(strsplit(rownames(genocount)[specifigenes[[cla]]], split='-'), function(x){ if (length(x)==2) return(x) else return(c(x, NA)) })), stringsAsFactors=F)
	spefamogs[,2] = as.numeric(spefamogs[,2]) ; colnames(spefamogs) = c("gene_family_id", "og_id")
	dbWriteTable(dbcon, "specific_genes", spefamogs, temporary=T)
	# choose adequate reference genome
	refgenome = NULL
	for (prefgenome in preferredgenomes){
		if (prefgenome %in% cladedefs[[cla]]$clade){
			refgenome = prefgenome
			break
		}
	}
	if (is.null(refgenome)){ refgenome = cladedefs[[cla]]$clade[1] }
	print(sprintf("%s: %s; ref: %s", cla, claname, refgenome))
	if ( ogcolid >= 0 ){
	# use ortholog classification of homologous genes
	spegeneinfo = dbGetQuery(dbcon, paste( c(
	 "SELECT gene_family_id, og_id, genomic_accession, locus_tag, cds_begin, cds_end, product",
	 "FROM coding_sequences",
	 "LEFT JOIN orthologous_groups USING (cds_code, gene_family_id)",
	 "INNER JOIN proteins USING (nr_protein_id)",
	 "INNER JOIN specific_genes USING (gene_family_id, og_id)",
	 "WHERE cds_code LIKE :c AND ( ortholog_col_id = :o OR ortholog_col_id IS NULL) ;"),
	 collapse=" "), params=list(c=sprintf("%s%%", refgenome), o=ogcolid))
	}else{
	# only use the homologous family mapping of genes (coarser homology mapping and stricter clade-specific gene finding)
	spegeneinfo = dbGetQuery(dbcon, paste( c(
	 "SELECT gene_family_id, genomic_accession, locus_tag, cds_begin, cds_end, product",
	 "FROM coding_sequences",
	 "INNER JOIN proteins USING (nr_protein_id)",
	 "INNER JOIN specific_genes USING (gene_family_id)",
	 "WHERE cds_code LIKE :c ESCAPE '@';"),
	 collapse=" "), params=list(c=sprintf("%s@_%%", refgenome)))
	}
	dbExecute(dbcon, "DROP TABLE specific_genes;")
	dbCommit(dbcon)
	write.table(spegeneinfo, file=nfoutspege, sep='\t', quote=F, row.names=F, col.names=T, append=T)
}


