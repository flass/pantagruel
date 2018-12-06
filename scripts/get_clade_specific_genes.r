#!/usr/bin/Rscript --vanilla
library('RSQLite')

genesetscopes = c("reprseq", "allseq")

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
if ( ogcolid < 0 ){
	print("will only use the homologous family mapping of genes (coarser homology mapping and stricter clade-specific gene finding", quote=F)
}else{
	print("use ortholog classification of homologous genes", quote=F)
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

cladedefwritesep = paste(c(rep('\t', 4), rep('  ...  \t', 3)), sep='', collapse='')

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
		write.table(spefamogs, file=file.path(diroutspegedetail, paste(bnoutspege, cla, "spegene_fams_ogids.tab", sep='_')), sep='\t', quote=F, row.names=F, col.names=T, append=F)
		# choose adequate reference genome
		refgenome = NULL
		for (prefgenome in preferredgenomes){
			if (prefgenome %in% cladedefs[[cla]]$clade){
				refgenome = prefgenome
				break
			}
		}
		if (is.null(refgenome)){ refgenome = min(cladedef$clade) }
		print(sprintf("%s: '%s'; %d clade-specific genes; ref genome: %s", cla, cladedef$name, ncsg, refgenome), quote=F)
		dbExecute(dbcon, "DROP TABLE IF EXISTS spegeneannots;")
		creaspegeneannots = paste( c(
		 "CREATE TABLE spegeneannots AS ",
		 "SELECT gene_family_id, og_id, genomic_accession, locus_tag, cds_code, cds_begin, cds_end, nr_protein_id, product",
		 "FROM (",
		 "  SELECT gene_family_id, og_id, ortholog_col_id, coding_sequences.* FROM specific_genes",
		 "   INNER JOIN orthologous_groups USING (gene_family_id, og_id)",
		 "   INNER JOIN coding_sequences USING (cds_code, gene_family_id)",
		 " UNION",
		 "  SELECT gene_family_id, og_id, ortholog_col_id, coding_sequences.* FROM specific_genes",
		 "   LEFT JOIN orthologous_groups USING (gene_family_id, og_id)",
		 "   INNER JOIN coding_sequences USING (gene_family_id)",
		 "   WHERE og_id IS NULL ) as famog2cds",
		 "INNER JOIN proteins USING (nr_protein_id)",
		 "INNER JOIN replicons USING (genomic_accession)",
		 "INNER JOIN assemblies USING (assembly_id)",
		 sprintf("WHERE code IN ( '%s' )", paste(cladedef$clade, collapse="','", sep='')),
		 "AND (ortholog_col_id = :o OR ortholog_col_id IS NULL) ;"),
		 collapse=" ")
#~ 		print(creaspegeneannots)
		dbExecute(dbcon, creaspegeneannots, params=list(o=ogcolid))
		testtest = dbGetQuery(dbcon, "select * from spegeneannots where gene_family_id='RHIZOC019762';")
		print('test RHIZOC019762')
		print(testtest)
		dbExecute(dbcon, "DROP TABLE specific_genes;")
		genesetclauses = list(sprintf("WHERE cds_code LIKE '%s@_%%' ESCAPE '@' AND", refgenome), "WHERE") ; names(genesetclauses) = genesetscopes
		for (genesetscope in genesetscopes){
			gsc = genesetclauses[[genesetscope]] # = "WHERE [clause AND]"
			spegeneinfo = dbGetQuery(dbcon, paste( c(
			"SELECT gene_family_id, og_id, cds_code, genomic_accession, locus_tag, cds_begin, cds_end, product", 
			"FROM spegeneannots", 
			gsc, "1 ORDER BY locus_tag ;"), collapse=" "))
			spegeneinfoplus = dbGetQuery(dbcon, paste( c(
			 "SELECT distinct gene_family_id, og_id, cds_code, genomic_accession, locus_tag, cds_begin, cds_end, product, interpro_id, interpro_description, go_terms, pathways",
			 "FROM spegeneannots",
			 "LEFT JOIN functional_annotations USING (nr_protein_id)",
			 "LEFT JOIN interpro_terms USING (interpro_id)", 
			 gsc, "1 ORDER BY locus_tag ;"), collapse=" "))
			spegallgoterms = dbGetQuery(dbcon, paste( c(
			 "SELECT distinct gene_family_id, og_id, cds_code, genomic_accession, locus_tag, go_id",
			 "FROM spegeneannots",
			 "LEFT JOIN functional_annotations USING (nr_protein_id)",
			 "LEFT JOIN interpro2GO USING (interpro_id)",
			 gsc, "go_id NOT NULL ORDER BY locus_tag ;"), collapse=" "))
			spegallpathways = dbGetQuery(dbcon, paste( c(
			 "SELECT distinct gene_family_id, og_id, cds_code, genomic_accession, locus_tag, pathway_db, pathway_id",
			 "FROM spegeneannots",
			 "LEFT JOIN functional_annotations USING (nr_protein_id)",
			 "LEFT JOIN interpro2pathways USING (interpro_id)",
			 gsc, "pathway_id NOT NULL ORDER BY locus_tag ;"), collapse=" "))
			if (genesetscope=="reprseq"){ write.table(spegeneinfo, file=nfoutspege, sep='\t', quote=F, row.names=F, col.names=T, append=T) }
			write.table(spegeneinfoplus, file=file.path(diroutspegedetail, paste(bnoutspege, cla, genesetscope, "details.tab", sep='_')), sep='\t', quote=F, row.names=F, col.names=T, append=F)
			write.table(spegallgoterms, file=file.path(diroutspegedetail, paste(bnoutspege, cla, genesetscope, "goterms.tab", sep='_')), sep='\t', quote=F, row.names=F, col.names=T, append=F)
			write.table(spegallpathways, file=file.path(diroutspegedetail, paste(bnoutspege, cla, genesetscope, "pathways.tab", sep='_')), sep='\t', quote=F, row.names=F, col.names=T, append=F)
		}
		dbExecute(dbcon, "DROP TABLE spegeneannots;")
		dbCommit(dbcon)
	}
}
#~ print("Warnings:", quote=F)
#~ print(warnings(), quote=F)
print(sprintf("wrote ouput in file '%s'", nfoutspege), quote=F)
