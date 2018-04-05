#!/usr/bin/Rscript

minfracNgenomeshow=0.75

readInput = function(prompt){
	if (interactive() ){
		val <- readline(prompt=prompt)
	}else{
		cat(prompt);
		val <- readLines("stdin",n=1)
	}
	return(val)
}

selectmingenomes = function(pseudocoremingenomes=NA, silent=F){
	while (is.na(pseudocoremingenomes) || pseudocoremingenomes < 0){
		pseudocoremingenomes = as.numeric(readInput(prompt="Please enter non-zero integer value for minimum of genomes represented in pseudo-core unicopy gene families: "))
	}
	if (!silent){ cat(sprintf("Selected %d as the minimum number of genomes to be represented in pseudo-core unicopy gene families\n", pseudocoremingenomes)) }
	return(pseudocoremingenomes)
}

getpseudocorefams = function(pseudocoremingenomes, countsbyfam){
	pseudocorefams = names(countsbyfam[countsbyfam >= pseudocoremingenomes])
	cat(sprintf("results in a set of %d pseudo-core unicopy gene families\n", length(pseudocorefams)))
	return(pseudocorefams)
}

selectMinGenomes = function(countmatrix, dirout, pseudocoremingenomes=NA, ngenomes=NA, interactive.X=FALSE, plot.PCoA=FALSE){
	countsbyfam = apply(countmatrix, 1, sum)
	if (is.numeric(ngenomes)){   N = ngenomes 
	}else{                       N = ncols(countmatrix) }
	cat(sprintf("number of unicopy gene families present in at least n genomes (out of %d):\n", N))
	print(cumsum(rev(table(countsbyfam)))[as.character(floor(N*minfracNgenomeshow):N)])

	pseudocoremingenomes = selectmingenomes(pseudocoremingenomes)	
	pseudocorefams = NULL
	pseudocore = list()
	pcmg = -1
	if (interactive.X){ X11(width=16, height=10) }
	while (pseudocoremingenomes!=pcmg){
		nftabout = file.path(dirout, sprintf("pseudo-core-%d-unicopy_families.tab", pseudocoremingenomes))
		nfpdfout = file.path(dirout, sprintf("pseudo-core-%d-unicopy_families.pdf", pseudocoremingenomes))
		if (!interactive.X){ pdf(nfpdfout, width=40, height=40) }
		pcmg = pseudocoremingenomes
		pseudocorefams = getpseudocorefams(pseudocoremingenomes, countsbyfam)
		cat("plotting heatmap... ")
		pseudocoremat = countmatrix[pseudocorefams,]
		heatmap(t(pseudocoremat), breaks=c(-0.5, 0.5, 1.5, N), col=c('white', 'black', 'red'), scale='none')
		if (interactive.X){
			cat("please verify that the distribution of markers per species is not too skewed (counts per species, white: 0, black: 1, red: >1)\n")
			message("Press Return To Continue") ; invisible(readLines("stdin", n=1))
		}
		if (plot.PCoA){
			cat("computing and plotting PCoA of gene species based on presence/absence... ")
			count.coa = dudi.coa(pseudocoremat, scannf=F, nf=2)
			if (interactive.X){ message("Press Return To Continue") ; invisible(readLines("stdin", n=1)) }
			s.label(count.coa$c1)
		}
		nmissing = apply(pseudocoremat, 2, function(x){ length(which(!x)) })
		barplot(nmissing, horiz=T, las=2, ylab='Species label', xlab='Nb. missing gene markers')
		if (interactive.X){ message("Press Return To Continue") ; invisible(readLines("stdin", n=1)) }
		barplot(nmissing[order(nmissing, decreasing=T)[1:min(20, N)]], horiz=T, las=2, ylab='Species label', xlab='Nb. missing gene markers')
		if (!interactive.X){ dev.off() }
		write(pseudocorefams, file=nftabout)
		cat(sprintf("Written list of pseudo-core unicopy gene families (with min. genome nb. = %d) and graphical representation of their distribution at:\n%s\n%s\n",
		 pseudocoremingenomes, nfpdfout, nftabout))
		cat("please confirm value for minimum of genomes represented in pseudo-core unicopy gene families [or stop by choosing 0]: \n")
		pseudocore[[pseudocoremingenomes]] = pseudocorefams
		pseudocoremingenomes = selectmingenomes(silent=T)
		if (pseudocoremingenomes==0){ break }
	}
	if (interactive.X){ dev.off() }
	return(pseudocore)
}

#~ ngenomes = as.numeric(Sys.getenv('ngenomes'))
#~ pseudocoremingenomes = as.numeric(Sys.getenv('pseudocoremingenomes'))
#~ protali = Sys.getenv('protali')
#~ nflasscode =  file.path(Sys.getenv('database'), 'genome_codes.tab')
#~ dirout = protali
#~ nffamgenomemat = file.path(protali, 'full_families_genome_counts-noORFans.mat')

cargs = commandArgs(trailingOnly=TRUE)

nffamgenomemat = cargs[1]
nflasscode = cargs[2]
dirout = cargs[3]
if (length(cargs) > 3){
	pseudocoremingenomes = as.numeric(cargs[4])
}else{
	pseudocoremingenomes = NA
}
if (length(cargs) > 4){
	interactive.X = as.logical(cargs[5])
}else{
	interactive.X = FALSE
}
cat("Loading matrix of gene families counts in genomes...\n")
genocount = data.matrix(read.table(file=nffamgenomemat))
lasscode = read.table(nflasscode, row.names=1, stringsAsFactors=F)
colnames(genocount) = lasscode[colnames(genocount),1]
onlyunicopy = apply(genocount, 1, function(x){ max(x)==1 })
genocountunicopy = genocount[onlyunicopy,]

pseudocore = selectMinGenomes(genocountunicopy, dirout, pseudocoremingenomes=pseudocoremingenomes, interactive.X=interactive.X)

pseudocoremingenomes = names(pseudocore)[length(pseudocore)]
cat(sprintf("Final choice of %d pseudo-core unicopy gene families (present in at least %d genomes).\n", length(pseudocore[[pseudocoremingenomes]]), as.numeric(pseudocoremingenomes))
nfdataout = file.path(dirout, "pseudo-core-all.RData")
save(genocountunicopy, pseudocore, file=nfdataout)
cat(sprintf("Saved data in file: '%s'.\n", nfdataout))

quit(status=pseudocoremingenomes, save='no')
