#!/usr/bin/Rscript --vanilla
library('RSQLite')
library('getopt')

#~ # clade definition file format example :
#~ 	name	maxabsin	maxpresout	clade	sisterclade
#~ cladeA	"P. endolithicum"	0	0	REJC140,REQ54	RFYW14,RHAB21,RKHAN,RHIZOB27,RNT25,RTCK,PSEPEL1,PSEPEL2,RHIMAR
#~ cladeB	"P. banfieldii"	0	0	RHIZOB27,RNT25,RTCK	RFYW14,RHAB21,RKHAN,REJC140,REQ54,PSEPEL1,PSEPEL2,RHIMAR
#~ cladeC	"P. halotolerans"	0	0	RFYW14,RHAB21,RKHAN	RHIZOB27,RNT25,RTCK,REJC140,REQ54,PSEPEL1,PSEPEL2,RHIMAR
#~ cladeD	"P. pelagicum"	0	0	PSEPEL1,PSEPEL2,RHIMAR	RFYW14,RHAB21,RKHAN,RHIZOB27,RNT25,RTCK,REJC140,REQ54
#~ cladeE	"Pban+Phalo"	0	0	RFYW14,RHAB21,RKHAN,RHIZOB27,RNT25,RTCK	REJC140,REQ54,PSEPEL1,PSEPEL2,RHIMAR
#~ cladeF	"Pban+Phalo+Pendo"	0	0	RFYW14,RHAB21,RKHAN,RHIZOB27,RNT25,RTCK,REJC140,REQ54	PSEPEL1,PSEPEL2,RHIMAR
#~ cladeG	"Psedorhizobium"	0	0	RFYW14,RHAB21,RKHAN,RHIZOB27,RNT25,RTCK,REJC140,REQ54,PSEPEL1,PSEPEL2,RHIMAR	

genesetscopes = c("reprseq", "allseq")

cargs = commandArgs(trailing=T)

##### main
spec = matrix(c(
  'gene_count_matrix',    'm', 1, "character", "path of metrix of counts of each gene family (or gene family_ortholog group id) (rows) in each genome (columns)",
  'sqldb',                'd', 1, "character", "path to SQLite database file",
  'clade_defs',           'C', 1, "character", "path to file describing clade composition; this must be a (tab-delimited) table file with named rows and a header with the following collumns: (mandatory:) clade, siterclade, (facultative:) name, maxabsin, maxpresout",
  'outrad',               'o', 1, "character", "path to output dir+file prefix for output files",
  'restrict_to_genomes',  'g', 2, "character", "(optional) path to file listing the genomes (UniProt-like codes) to which the analysis will be restricted",
  'og_col_id',            'c', 2, "integer",   "orthologous group collection id in SQL database; if not provided, will only use the homologous family mapping of genes (coarser homology mapping, meaning stricter clade-specific gene definition)",
  'ass_to_code',          'a', 2, "character", "(optional) path to file providing correspondency between assembly ids and UniProt-like genome codes; only if the input matrix has assembly ids in column names (deprecated)",
  'preferred_genomes',    'p', 2, "character", "(optional) comma-separated list of codes of genomes which CDS info will be reported in reference tables",
  'interesting_families', 'f', 2, "character", "(optional) comma-separated list of gene families for which detail of presence/absence distribution will be printed out"
), byrow=TRUE, ncol=5);
opt = getopt(spec, opt=commandArgs(trailingOnly=T))

nffamgenomemat = opt$gene_count_matrix
sqldb = opt$sqldb
nfcladedef = opt$clade_defs
outfilerad = opt$outrad
ogcolid = opt$og_col_id
nfrestrictlist = opt$restrict_to_genomes
nflasscode = opt$ass_to_code
if (!is.null(opt$preferred_genomes)){
	preferredgenomes = strsplit(opt$preferred_genomes, split=',')[[1]]
}else{
	preferredgenomes = c()
}
if (!is.null(opt$interesting_families)){
	interstfams = strsplit(opt$interesting_families, split=',')[[1]]
}else{
	interstfams = c()
}
if ( is.null(ogcolid) | ogcolid < 0 ){
	print("will only use the homologous family mapping of genes (coarser homology mapping and stricter clade-specific gene finding)", quote=F)
}else{
	print("use ortholog classification of homologous genes", quote=F)
}


# output files
nfabspresmat = sprintf("%s_gene_abspres.mat.RData", outfilerad)
nfoutspege = sprintf("%s_specific_genes.tab", outfilerad)
bnoutspege = sub(".tab$", "", basename(nfoutspege))
diroutspegedetail = sprintf("%s_specific_genes.tables_byclade_goterms_pathways", outfilerad)
dir.create(diroutspegedetail, showWarnings=F)

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
	if (!is.null(nfrestrictlist)){
		restrictgenomelist = readLines(nfrestrictlist)
		if (!is.null(nflasscode)){
			lasscode = read.table(nflasscode, row.names=1, stringsAsFactors=F)
			colnames(genocount) = lasscode[colnames(genocount),1]
		}
	#~ 	print(setdiff(restrictgenomelist, colnames(genocount)))
		genocount = genocount[,restrictgenomelist]
		gc()
	}
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

today = Sys.Date()
write( paste("#", c(format(today, format="%B %d %Y"), "Pantagruel version:", system("cd ${ptgscripts} ; git log | head -n 3", intern=T), "- - - - -")), file=nfoutspege, append=F)

for (i in 1:length(cladedefs)){
	cla = names(cladedefs)[i]
	cladedef = cladedefs[[cla]]
	write(sprintf("# %s %s", cla, cladedef$name), file=nfoutspege, append=T)
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
