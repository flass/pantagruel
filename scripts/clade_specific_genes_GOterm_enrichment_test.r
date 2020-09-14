#!/usr/bin/env Rscript
library(topGO)
library(getopt)

parseGeneTermTable = function(filename){
	gt =  read.table(filename, sep='\t', header=F, stringsAsFactors=F)
	colnames(gt) = c('gene', 'term')
	ugene = unique(gt$gene)
	gene2term = lapply(ugene, function(g){ gt$term[gt$gene==g] })
	names(gene2term) = ugene
	return(gene2term)
}

parseGAF = function(filename, gene.name.field=2, GO.term.field=5, header=F){
	# create list map of gene to annotated GO terms from table, by default inthe GAF 2.1 format http://www.geneontology.org/page/go-annotation-file-gaf-format-21
	gaf = read.table(filename, skip=1, sep='\t', header=header, stringsAsFactors=F)[,c(gene.name.field,GO.term.field)]
	colnames(gaf) = c('gene', 'term')
	ugene = unique(gaf$gene)
	gene2term = lapply(ugene, function(g){ gaf$term[gaf$gene==g] })
	names(gene2term) = ugene
	return(gene2term)
}

parseGeneAnnotFile = function(filename){
	if (grepl("\\.gaf$", filename)){
		gene2GO = parseGAF(filename)
	}else{
		gene2GO = parseGeneTermTable(filename)
	}
	return(gene2GO)
}

combineGeneAnnotLists = function(lgene2GO){
	if (length(lgene2GO)==1){
		return(lgene2GO[[1]])
	}else{
		cgene2GO = do.call(c, lgene2GO) ; uncgene2GO = unique(names(cgene2GO))
		gene2GO = lapply(uncgene2GO, function(x){ cgene2GO[[x]] }) ; names(gene2GO) = uncgene2GO
		return(gene2GO)
	}
}

makeGeneList = function(studyset, populationset){
	geneList = factor(as.integer(populationset %in% studyset))
	names(geneList) = populationset
return(geneList)
}

##### main
spec = matrix(c(
  'study_annots',         's', 1, "character", "path(s) or glob-able pattern of paths of files containing the gene to GO term annotations for the study (i.e. test) set; files need to be table, either with just two columns (for 1) gene ids and 2) go terms, no header), or in the GAF 2.1 format, in which case the file must bear the '.gaf' extension",
  'population_annots',    'p', 1, "character", "path(s) or glob-able pattern of paths of files containing the gene to GO term annotations for the population (i.e. reference) set; format as for study set",
  'out',                  'o', 1, "character", "path to output dir/file (depending on if input study set is multiple/single)",
  'combine_study_files',  'c', 2, "integer",   "should the files specified for the study sets (1) be combined into one set, or (0, default) treated iteratively for separate tests?",
  'ontology',             'g', 2, "character", "indicate which type of ontology is considered among 'MF', 'CC' or 'BP' (default)",
  'algo',                 'a', 2, "character", "algorithm used for GO term enrichment with topGO (default 'weight01; for others see https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf)",
  'stat',                 't', 2, "character", "statistic used for GO term enrichment with topGO (default 'fisher; for others see https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf)"
), byrow=TRUE, ncol=5);
opt = getopt(spec, opt=commandArgs(trailingOnly=T))

print("This is 'clade_specific_genes_GOterm_enrichment_test.r', with options:")
print(opt)
## argument parsing
# combine_study_files
if (!is.null(opt$combine_study_files)){ csf = as.logical(opt$combine_study_files)
}else{ csf = FALSE }
if (!is.null(opt$study_annots)){
	lnfstudyset = Sys.glob(strsplit(opt$study_annots, split=',')[[1]])
	lstudy.gene2GO = lapply(lnfstudyset, parseGeneAnnotFile)
	if (csf){
		study.gene2GO = combineGeneAnnotLists(lstudy.gene2GO)
		lstudy.gene2GO = list( study.gene2GO )
		names(lstudy.gene2GO) = "combined_study_sets"
	}else{
		bn = basename(lnfstudyset)
		if (length(bn)==length(unique(bn))){
			names(lstudy.gene2GO) = bn
		}else{
			names(lstudy.gene2GO) = lnfstudyset
		}
	}
}else{
	stop("It is required to specify one or several file paths for the study set(s) with option --study_annots")
}
# pop_annots
if (!is.null(opt$population_annots)){
	lnfpopset = Sys.glob(strsplit(opt$population_annots, split=',')[[1]])
	lpop.gene2GO = lapply(lnfpopset, parseGeneAnnotFile)
	pop.gene2GO = combineGeneAnnotLists(lpop.gene2GO)
}else{
	stop("It is required to specify one or several file paths for the population set with option --population_annots")
}
# out
if (!is.null(opt$out)){
	outloc = opt$out
	if (length(lstudy.gene2GO)>1){ outd = outloc
	}else{ outd = dirname(outloc) }
	if (!file.exists(outd)){
		dir.create(outd)
	}else{
		if (!file.info(outd)$isdir){
			stop(sprintf("'%s' is not a directory", outd))
		}
	}
}else{
	stop("It is required to specify an output location with --out")
}
# ontology
if (!is.null(opt$ontology)){ ontology = opt$ontology }else{ ontology = "BP" }
# algo
if (!is.null(opt$algo)){ tgalgo = opt$algo }else{ tgalgo = "weight01" }
# stat
if (!is.null(opt$stat)){ tgstat = opt$stat }else{ tgstat = "Fisher" }

# run analysis
for (study.name in names(lstudy.gene2GO)){
	study.gene2GO = lstudy.gene2GO[[study.name]]
	geneList = makeGeneList(names(study.gene2GO), names(pop.gene2GO))
	print("head(geneList):")
	print(head(geneList))
	print("head(pop.gene2GO):")
	print(head(pop.gene2GO))
	print("create topGOdata object")
	GOdata = new("topGOdata", ontology=ontology, allGenes=geneList, annot=annFUN.gene2GO, gene2GO=pop.gene2GO)
	print(sprintf("run term enrichment test (algo: %s; statistic: %s)", tgalgo, tgstat))
	enrich = runTest(GOdata, algorithm=tgalgo, statistic=tgstat)
	print("generate annotation table for enriched terms")
	enriched = GenTable(GOdata, weight = enrich, orderBy = "weight")
	if (length(lstudy.gene2GO)>1){ nfout = file.path(outd, paste(study.name, "enriched", sep='_'))
	}else{ nfout = outloc }
	print(sprintf("write annotation table to '%s'", nfout))
	write.table(enriched, file=nfout, sep='\t', col.names=T, row.names=F)
}
