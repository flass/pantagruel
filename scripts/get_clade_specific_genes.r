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
infilerad = cargs[4]
outfilerad = cargs[5]
# genome subset-specific input files
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
if (length(cargs)>7){
	interstfams = strsplit(cargs[8], split=',')[[1]]
}else{
	interstfams = c()
}

nfrestrictlist = sprintf("%s_genome_codes", infilerad)
nfcladedef = sprintf("%s_clade_defs", infilerad)
# output files
nfabspresmat = sprintf("%s_gene_abspres.mat.RData", outfilerad)
nfoutspege = sprintf("%s_specific_genes.tab", outfilerad)
bnoutspege = sub(".tab$", "", basename(nfoutspege))
diroutspegedetail = sprintf("%s_specific_genes.tables_byclade_goterms_pathways", outfilerad)
dir.create(diroutspegedetail, showWarnings=F)

# # clade definition object structure example :
# if (!file.exists(nfcladedef)){
# 	clade0 = c('RHIZOB27', 'RNT25')
# 	outgroup0 = 'RHISP2'
# 	clade1 = c(clade0, outgroup0)
# 	clade2 = c('RKHAN', 'RHAB21', 'RFYW14')
# 	clade12 = c(clade1, clade2)
# 	clade3 = c('REQ54', 'REJC140')
# 	clade123 = c(clade12, clade3)
# 	clade4 = c('PSEPEL1', 'PSEPEL2', 'RHIMAR')
# 	clade1234 = c(clade123, clade4)
# 	outgroup = 'RHIZOB54'
# 	clades = list(clade0=clade0, outgroup0=outgroup0, clade1=clade1, clade2=clade2, clade12=clade12, clade3=clade3, clade123=clade123, clade4=clade4, clade1234=clade1234, outgroup=outgroup)
# 	sisterclades = list(clade0=outgroup0, outgroup0=clade0, clade1=clade2, clade2=clade1, clade12=clade3, clade3=clade12, clade123=clade4, clade4=clade123, clade1234=outgroup, outgroup=clade1234)
# 	write.table(t(sapply(names(clades), function(cla){ sapply(list(clades, sisterclades), function(x){ paste(x[[cla]], collapse=',') }) })),
# 	 sep='\t', quote=F, row.names=names(clades), col.names=c('clade', 'sisterclade'),
# 	 file=nfcladedef)
# }
cladedefcsv = read.table(nfcladedef, sep='\t', header=T, row.names=1, stringsAsFactors=F)
cladedefs = apply(cladedefcsv[,c('clade', 'sisterclade')], 1, strsplit, split=',')

for (i in 1:length(cladedefs)){
	cla = names(cladedefs)[i]
	cladedefs[[cla]]$name = ifelse(!is.null(cladedefcsv$name), cladedefcsv$name[i], "")
	cladedefs[[cla]]$maxabsin = ifelse(!is.null(cladedefcsv$maxabsin), cladedefcsv$maxabsin[i], 0)
	cladedefs[[cla]]$maxpresout = ifelse(!is.null(cladedefcsv$maxpresout), cladedefcsv$maxpresout[i], 0)
}

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
#~ 	allpres = apply(genocount[,cladedef$clade,drop=F], 1, function(x){ all(x>0) })
#~ 	allabssis = apply(genocount[,cladedef$sisterclade,drop=F], 1, function(x){ all(x==0) })
	allpres = apply(genocount[,cladedef$clade,drop=F], 1, function(x){ length(which(x==0))<=cladedef$maxabsin })
	allabssis = apply(genocount[,cladedef$sisterclade,drop=F], 1, function(x){ length(which(x>0))<=cladedef$maxpresout })
	return(which(allpres & allabssis))
})

dbcon = dbConnect(SQLite(), sqldb)

cladedefwritesep = paste(rep('  ...  \t', 6), sep='', collapse='')

for (i in 1:length(cladedefs)){
	cla = names(cladedefs)[i]
	cladedef = cladedefs[[cla]]
	a = as.logical(i-1)
	write(sprintf("# %s %s", cla, cladedef$name), file=nfoutspege, append=a)
	write( paste(sprintf("# gene families present in all genomes but %d of clade:", cladedef$maxabsin), cladedefcsv[cla,'clade'], sep=cladedefwritesep), file=nfoutspege, append=T)
	write( paste(sprintf("# and absent in all but %d genomes of sister clade:", cladedef$maxpresout), cladedefcsv[cla,'sisterclade'], sep=cladedefwritesep), file=nfoutspege, append=T)
	ncsg = length(specifigenes[[cla]])
	if (ncsg==0){
		write("# no specific gene found", file=nfoutspege, append=T)
		print(sprintf("%s: %s; no specific gene found", cla, cladedef$name))
	}else{
		for (interstfam in interstfams){
			if (interstfam %in% rownames(genocount)[specifigenes[[cla]]]){
				print(sprintf("%s family:", interstfam), quote=F)
				miss = which(genocount[interstfam,cladedef$clade,drop=F]==0)
				print(sprintf("  missing in %d clade genomes: %s", length(miss), paste(cladedef$clade[miss], sep=' ', collapse=' ')), quote=F)
				extra = which(genocount[interstfam,cladedef$sisterclade,drop=F]>0)
				print(sprintf("  present in %d sister clade genomes: %s", length(extra), paste(cladedef$sisterclade[extra], sep=' ', collapse=' ')), quote=F)
			}
		}
		spefamogs = as.data.frame(t(sapply(strsplit(rownames(genocount)[specifigenes[[cla]]], split='-'), function(x){ if (length(x)==2) return(x) else return(c(x, NA)) })), stringsAsFactors=F)
		spefamogs[,2] = as.numeric(spefamogs[,2]) ; colnames(spefamogs) = c("gene_family_id", "og_id")
		dbBegin(dbcon)
		dbWriteTable(dbcon, "specific_genes", spefamogs, temporary=T)
		# choose adequate reference genome
		refgenome = NULL
		for (prefgenome in preferredgenomes){
			if (prefgenome %in% cladedefs[[cla]]$clade){
				refgenome = prefgenome
				break
			}
		}
		if (is.null(refgenome)){ refgenome = cladedef$clade[1] }
		print(sprintf("%s: %s; %d clade-specific genes; ref: %s", cla, cladedef$name, ncsg, refgenome), quote=F)
		if ( ogcolid >= 0 ){
		# use ortholog classification of homologous genes
		dbExecute(dbcon, paste( c(
		 "CREATE TABLE spegeneannots AS SELECT gene_family_id, og_id, genomic_accession, locus_tag, cds_begin, cds_end, nr_protein_id, product",
		 "FROM coding_sequences",
		 "LEFT JOIN orthologous_groups USING (cds_code, gene_family_id)",
		 "INNER JOIN proteins USING (nr_protein_id)",
		 "INNER JOIN specific_genes USING (gene_family_id, og_id)",
		 "WHERE cds_code LIKE :c AND ( ortholog_col_id = :o OR ortholog_col_id IS NULL) ;"),
		 collapse=" "), params=list(c=sprintf("%s%%", refgenome), o=ogcolid))
		}else{
		# only use the homologous family mapping of genes (coarser homology mapping and stricter clade-specific gene finding)
		dbExecute(dbcon, paste( c(
		 "CREATE TABLE spegeneannots AS SELECT gene_family_id, genomic_accession, locus_tag, cds_begin, cds_end, product, nr_protein_id, product",
		 "FROM coding_sequences",
		 "INNER JOIN proteins USING (nr_protein_id)",
		 "INNER JOIN specific_genes USING (gene_family_id)",
		 "WHERE cds_code LIKE :c ESCAPE '@';"),
		 collapse=" "), params=list(c=sprintf("%s@_%%", refgenome)))
		}
		spegeneinfo = dbGetQuery(dbcon, "SELECT gene_family_id, og_id, genomic_accession, locus_tag, cds_begin, cds_end, product FROM spegeneannots;")
		spegeneinfoplus = dbGetQuery(dbcon, paste( c(
		 "SELECT distinct gene_family_id, og_id, genomic_accession, locus_tag, cds_begin, cds_end, product, interpro_id, interpro_description, go_terms, pathways",
		 "FROM spegeneannots",
		 "LEFT JOIN functional_annotations USING (nr_protein_id)",
		 "LEFT JOIN interpro_terms USING (interpro_id) ;"), collapse=" "))
		spegallgoterms = dbGetQuery(dbcon, paste( c(
		 "SELECT distinct gene_family_id, og_id, genomic_accession, locus_tag, go_id",
		 "FROM spegeneannots",
		 "LEFT JOIN functional_annotations USING (nr_protein_id)",
		 "LEFT JOIN interpro2GO USING (interpro_id) ;"), collapse=" "))
		spegallpathways = dbGetQuery(dbcon, paste( c(
		 "SELECT distinct gene_family_id, og_id, genomic_accession, locus_tag, pathway_db, pathway_id",
		 "FROM spegeneannots",
		 "LEFT JOIN functional_annotations USING (nr_protein_id)",
		 "LEFT JOIN interpro2pathways USING (interpro_id) ;"), collapse=" "))
		dbExecute(dbcon, "DROP TABLE specific_genes;")
		dbExecute(dbcon, "DROP TABLE spegeneannots;")
		dbCommit(dbcon)
		write.table(spegeneinfo, file=nfoutspege, sep='\t', quote=F, row.names=F, col.names=T, append=T)
		write.table(spegeneinfoplus, file=file.path(diroutspegedetail, paste(bnoutspege, cla, "details.tab", sep='_')), sep='\t', quote=F, row.names=F, col.names=T, append=F)
		write.table(spegallgoterms, file=file.path(diroutspegedetail, paste(bnoutspege, cla, "goterms.tab", sep='_')), sep='\t', quote=F, row.names=F, col.names=T, append=F)
		write.table(spegallpathways, file=file.path(diroutspegedetail, paste(bnoutspege, cla, "pathways.tab", sep='_')), sep='\t', quote=F, row.names=F, col.names=T, append=F)
	}
}
#~ print("Warnings:", quote=F)
#~ print(warnings(), quote=F)
print(sprintf("wrote ouput in file '%s'", nfoutspege), quote=F)
