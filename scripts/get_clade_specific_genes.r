#!/usr/bin/Rscript --vanilla
library('RSQLite')

# env variables relating to the full Pantagruel pangenome db, expected to be set
nflasscode =  file.path(Sys.getenv('database'), 'genome_codes.tab')
nffamgenomemat = file.path(Sys.getenv('protali'), 'full_families_genome_counts-noORFans.mat')
sqldb = Sys.getenv('sqldb')
# genome subset-specific input files
cargs = commandArgs(trailing=T)
if (length(cargs)>0)
	filerad = cargs[1]
}else{
	dircore = Sys.getenv('coregenome')
	ana = Sys.getenv('NT26like')
	filerad = file.path(dircore, ana, ana)
}
nfrestrictlist = sprintf("%s_genome_codes", filerad)
nfcladedef = sprintf("%s_clade_defs", filerad)
# output files
nfabspresmat = sprintf("%s_gene_abspres.mat.RData", filerad)
nfoutspege = sprintf("%s_specific_genes.tab", filerad)

# clade definitions
if (!file.exists(nfcladedef)){
#~ 	clade0 = c('RHIZOB16', 'RNT25')
#~ 	clade0 = c('RHIZOB27', 'RHISP1')
	clade0 = c('RHIZOB27', 'RNT25')
#~ 	outgroup0 = 'RTCK'
	outgroup0 = 'RHISP2'
	clade1 = c(clade0, outgroup0)
	clade2 = c('RKHAN', 'RHAB21', 'RFYW14')
#~ 	clade2 = c('RHISP3', 'RHIHAL', 'RHIFLA')
	clade12 = c(clade1, clade2)
	clade3 = c('REQ54', 'REJC140')
#~ 	clade3 = c('RHIEND1', 'RHIEND2')
	clade123 = c(clade12, clade3)
#~ 	clade4 = c('NEOGAL1', 'NEOGAL2', 'RHIZOB17')
	clade4 = c('PSEPEL1', 'PSEPEL2', 'RHIMAR')
	clade1234 = c(clade123, clade4)
#~ 	outgroup = 'RHIZOB33'
	outgroup = 'RHIZOB54'
	clades = list(clade0=clade0, outgroup0=outgroup0, clade1=clade1, clade2=clade2, clade12=clade12, clade3=clade3, clade123=clade123, clade4=clade4, clade1234=clade1234, outgroup=outgroup)
	sisterclades = list(clade0=outgroup0, outgroup0=clade0, clade1=clade2, clade2=clade1, clade12=clade3, clade3=clade12, clade123=clade4, clade4=clade123, clade1234=outgroup, outgroup=clade1234)
	write.table(t(sapply(names(clades), function(cla){ sapply(list(clades, sisterclades), function(x){ paste(x[[cla]], collapse=',') }) })),
	 sep='\t', quote=F, row.names=names(clades), col.names=c('clade', 'sisterclade'),
	 file=nfcladedef)
}
cladedefcsv = read.table(nfcladedef, sep='\t', header=T, row.names=1, stringsAsFactors=F)
cladedefs = apply(cladedefcsv, 1, strsplit, split=',')

# load gene presence / absence data
if (!file.exists(nfabspresmat)){
	load(nfabspresmat)
else{
	genocount = data.matrix(read.table(file=nffamgenomemat))
	lasscode = read.table(nflasscode, row.names=1, stringsAsFactors=F)
	colnames(genocount) = lasscode[colnames(genocount),1]
	restrictgenomelist = readLines(nfrestrictlist)
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
	cla = names(cladedefs)[i] ; print(cla)
	a = as.logical(i-1)
	write(paste(c("# gene families present in all genomes of clade:", "and absent in all genomes of sister clade:"), cladedefcsv[cla,], sep=' ', collapse='; '),
	 file=nfoutspege, append=a)
	dbBegin(dbcon)
	dbWriteTable(dbcon, "specific_genes", data.frame(gene_family_id=rownames(genocount)[specifigenes[[cla]]]), temporary=T)
	spegeneinfo = dbGetQuery(dbcon, paste( c(
	 "SELECT gene_family_id, genomic_accession, locus_tag, cds_begin, cds_end, product",
	 "FROM coding_sequences",
	 "INNER JOIN proteins USING (nr_protein_id)",
	 "INNER JOIN specific_genes USING (gene_family_id)",
	 "WHERE cds_code LIKE :c ;"),
	 collapse=" "), params=list(c=sprintf("%s%%", cladedefs[[cla]]$clade[1])))
	dbExecute(dbcon, "DROP TABLE specific_genes;")
	dbCommit(dbcon)
	write.table(spegeneinfo, file=nfoutspege, sep='\t', quote=F, row.names=F, col.names=T, append=T)
}


